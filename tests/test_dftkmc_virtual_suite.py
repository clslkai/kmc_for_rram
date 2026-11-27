"""Protocol and inference checks for the synthetic virtual-data experiments.

These tests do not demand that random 95% intervals always contain a reference.
The artificial replica dictionaries test reporting logic only; they are not KMC
results or material data. Actual KMC is used only for a small paired-field check.
"""
from collections import Counter
from copy import deepcopy
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from dftkmc.geometry import Lattice
from dftkmc.mockdata import import_mock_bundle, make_mock_bundle
from dftkmc.virtual_suite import (
    _simulate_replica, build_suite_plan, run_virtual_suite, summarize_suite,
)


KB_EV_K = 8.617333262145e-5


def estimator_fixtures(plan, diffusion_coefficient=None):
    """Balanced, correlated tracers with a prescribed unbiased scalar D estimate.

    No event sequence is claimed for these abstract statistical fixtures.
    Neutral-field records deliberately preserve the planned pairing.
    """
    records = []
    for condition in plan["conditions"]:
        n = condition["n_replicas"]
        if n % 2:
            raise ValueError("This balanced fixture requires an even replica count")
        duration = condition["duration_s"]
        if condition["family"] == "diffusion":
            target = (diffusion_coefficient(condition) if diffusion_coefficient else
                      condition["reference_finite_time"]["trace_over_3_m2_s"] * 1e-12)
        else:
            target = 1e-25
        # For n/2 displacements at +/-a along x, the unbiased trace(Cov)/(6t)
        # equals n/(n-1) * a**2 * 1e-20 / (6t).
        amplitude = math.sqrt(target * 6.0 * duration * 1e20 * (n - 1) / n)
        indices = list(range(condition["n_vacancies"]))
        for stream in condition["replica_seeds"]:
            sign = -1 if stream["replica_index"] < n // 2 else 1
            displacement = np.tile([sign * amplitude, 0.0, 0.0], (len(indices), 1))
            records.append(dict(
                condition_id=condition["id"], replica_index=stream["replica_index"],
                initialization_seed=stream["initialization_seed"], seed=stream["trajectory_seed"],
                duration_s=duration, charge=condition["charge"], family=condition["family"],
                status="completed_duration", time_s=duration, target_reached=True, events=7,
                initial_vacancy_indices=indices, final_vacancy_indices=indices,
                displacements_A=displacement.tolist(), channel_counts={"statistical_fixture_only": 7},
                final=dict(energy_eV=0.0, msd_A2=amplitude**2, pair_count=0,
                           largest_cluster=1, events=7, z_displacement_A=0.0),
            ))
    return records


class VirtualSuiteValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = import_mock_bundle(make_mock_bundle())
        cls.standard = build_suite_plan(cls.catalog, profile="standard", base_seed=20261111)
        cls.smoke = build_suite_plan(cls.catalog, profile="smoke", base_seed=20261111)

    def test_standard_design_has_fixed_times_and_equilibrium_weights(self):
        plan = self.standard
        self.assertEqual(plan["n_conditions"], 24)
        self.assertEqual(plan["n_trajectory_executions"], 480)
        self.assertEqual(plan["n_distinct_stream_pairs"], 448)
        self.assertEqual(Counter(c["family"] for c in plan["conditions"]),
                         {"diffusion": 12, "mechanism": 6, "field": 6})
        seen = set()
        for condition in plan["conditions"]:
            self.assertGreater(condition["duration_s"], 0.0)
            if condition["family"] != "diffusion":
                continue
            seen.add((condition["charge"], condition["temperature_K"], condition["target_expected_hops"]))
            self.assertEqual(condition["n_replicas"], 24)
            self.assertEqual(condition["n_vacancies"], 1)
            self.assertEqual(condition["initial_mode"], "equilibrium_single")
            rate = condition["reference_long_time"]["mean_jump_rate_Hz"]
            self.assertAlmostEqual(condition["duration_s"] * rate,
                                   condition["target_expected_hops"], places=10)
            self.assertEqual(condition["reference_finite_time"]["duration_s"], condition["duration_s"])
            self.assertAlmostEqual(condition["reference_finite_time"]["expected_hop_count"],
                                   condition["target_expected_hops"], places=9)
            energy = self.catalog.site_energies(condition["charge"])
            boltzmann = np.exp(-(energy - energy.min()) / (KB_EV_K * condition["temperature_K"]))
            boltzmann /= boltzmann.sum()
            np.testing.assert_allclose(condition["basis_probabilities"], boltzmann, rtol=1e-12)
        self.assertEqual(seen, {(q, t, h) for q in (0, 2) for t in (700., 800., 900.) for h in (50, 200)})
        for family in ("mechanism", "field"):
            group = [c for c in plan["conditions"] if c["family"] == family]
            self.assertEqual(len({c["duration_s"] for c in group}), 1)
            self.assertTrue(all(c["n_replicas"] == 16 for c in group))

    def test_random_stream_reuse_is_only_the_declared_neutral_pairing(self):
        seen = {}
        bootstrap = []
        for condition in self.standard["conditions"]:
            self.assertEqual(len(condition["replica_seeds"]), condition["n_replicas"])
            bootstrap.append(condition["bootstrap_seed"])
            for stream in condition["replica_seeds"]:
                pair = (stream["initialization_seed"], stream["trajectory_seed"])
                self.assertNotEqual(*pair)
                seen.setdefault(pair, []).append((condition["family"], condition["charge"], stream["replica_index"]))
        self.assertEqual(len(bootstrap), len(set(bootstrap)))
        self.assertEqual(len(seen), 448)
        reused = [uses for uses in seen.values() if len(uses) > 1]
        self.assertEqual(len(reused), 16)
        for uses in reused:
            self.assertEqual(len(uses), 3)
            self.assertEqual(len(set(uses)), 1)
            self.assertEqual(uses[0][:2], ("field", 0))
        all_initial = {a for a, _ in seen}
        all_trajectory = {b for _, b in seen}
        self.assertFalse(all_initial & all_trajectory)
        self.assertFalse(set(bootstrap) & (all_initial | all_trajectory))

    def test_event_cap_does_not_choose_the_physical_observation_time(self):
        capped = build_suite_plan(self.catalog, profile="smoke", base_seed=20261111, max_events=1)
        for before, after in zip(self.smoke["conditions"], capped["conditions"]):
            self.assertEqual(before["id"], after["id"])
            self.assertEqual(before["duration_s"], after["duration_s"])
            self.assertEqual(before["replica_seeds"], after["replica_seeds"])
            self.assertEqual(after["max_events"], 1)

    def test_confidence_interval_misses_are_diagnostic_not_a_failure_gate(self):
        result = summarize_suite(self.smoke, estimator_fixtures(self.smoke))
        self.assertTrue(result["all_trajectories_complete"])
        self.assertTrue(result["neutral_field_paired_invariance"])
        self.assertTrue(result["deterministic_checks_passed"])
        comparisons = result["reference_comparisons"]
        self.assertTrue(comparisons)
        self.assertTrue(all(c["ci_contains_finite_time_reference"] is False for c in comparisons))
        self.assertEqual(result["statistical_diagnostics"]["n_ci_cover_reference"], 0)
        for condition in self.smoke["conditions"]:
            if condition["family"] == "diffusion":
                comparison = next(c for c in comparisons if c["condition_id"] == condition["id"])
                self.assertEqual(comparison["reference_finite_time_m2_s"],
                                 condition["reference_finite_time"]["trace_over_3_m2_s"])

    def test_one_censored_replica_suppresses_the_group_inference(self):
        records = estimator_fixtures(self.smoke)
        censored = records[0]
        censored.update(status="event_limit_before_end", target_reached=False,
                        time_s=censored["duration_s"] / 2.0)
        result = summarize_suite(self.smoke, records)
        self.assertFalse(result["all_trajectories_complete"])
        self.assertFalse(result["deterministic_checks_passed"])
        group = next(g for g in result["groups"] if g["condition"]["id"] == censored["condition_id"])
        stats = group["statistics"]
        self.assertEqual(stats["n_incomplete_replicas"], 1)
        self.assertFalse(stats["diffusion"]["available"])
        self.assertTrue(all(value["ci"] is None for value in stats["scalar_estimates"].values()))
        self.assertIsNone(group["reference_comparison"]["estimate_m2_s"])
        self.assertIsNone(group["reference_comparison"]["ci_contains_finite_time_reference"])

    def test_missing_replicas_do_not_silently_become_a_smaller_complete_set(self):
        records = estimator_fixtures(self.smoke)
        with self.assertRaisesRegex(ValueError, "Missing replica"):
            summarize_suite(self.smoke, records[1:])

    def test_actual_neutral_field_paths_preserve_pairing_and_independent_initialization(self):
        conditions = [deepcopy(c) for c in self.smoke["conditions"]
                      if c["family"] == "field" and c["charge"] == 0]
        reference = next(c["reference_long_time"] for c in self.smoke["conditions"]
                         if c["family"] == "diffusion" and c["charge"] == 0)
        # A priori duration gives enough hops to exercise the equality meaningfully.
        for condition in conditions:
            condition["duration_s"] = 100.0 / reference["mean_jump_rate_Hz"]
        records = [_simulate_replica(self.catalog, c, c["replica_seeds"][0]) for c in conditions]
        baseline = records[0]
        self.assertTrue(baseline["target_reached"])
        self.assertGreater(baseline["events"], 0)
        for record in records[1:]:
            for key in ("seed", "initialization_seed", "initial_vacancy_indices", "final_vacancy_indices",
                        "events", "time_s", "displacements_A", "channel_counts"):
                self.assertEqual(record[key], baseline[key])
        lattice = Lattice(repeats=tuple(conditions[0]["repeats"]))
        energy = self.catalog.site_energies(0)[lattice.basis_id]
        probability = np.exp(-(energy - energy.min()) / (KB_EV_K * 800.0))
        probability /= probability.sum()
        initial_rng = np.random.default_rng(baseline["initialization_seed"])
        independently_selected = initial_rng.choice(lattice.n, size=1, replace=False, p=probability)
        np.testing.assert_array_equal(baseline["initial_vacancy_indices"], independently_selected)
        self.assertNotEqual(baseline["initialization_seed"], baseline["seed"])

    def test_apparent_slope_uses_correct_sign_units_and_keeps_reference_separate(self):
        plan = deepcopy(self.standard)
        plan["profile_controls"]["bootstrap_samples"] = 100
        activation = {0: 0.37, 2: 0.41}
        records = estimator_fixtures(plan, lambda c: 1e-10 * math.exp(
            -activation[c["charge"]] / (KB_EV_K * c["temperature_K"])))
        result = summarize_suite(plan, records)
        slopes = result["apparent_arrhenius_slopes"]
        for q in (0, 2):
            protocol = next(s for s in slopes if s["charge"] == q and s["source"] == "finite_time_KMC_protocol")
            self.assertAlmostEqual(protocol["apparent_activation_eV"], activation[q], places=11)
            self.assertEqual(protocol["n_temperatures"], 3)
            self.assertTrue(any(s["charge"] == q and s["source"] == "exact_long_time_reference" for s in slopes))
            self.assertIn("not an individual NEB barrier", protocol["interpretation"])

    def test_virtual_suite_requires_explicit_synthetic_permission(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "must_not_run"
            with self.assertRaises(ValueError):
                run_virtual_suite(self.catalog, output, profile="smoke", workers=1)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
