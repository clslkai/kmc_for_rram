"""Output-level acceptance for the illustrative variable-population device.

The fast synthetic fixture deliberately exercises all three event types. Its
numbers are test inputs, not a material calibration or a forming benchmark.
"""
import contextlib
import csv
from dataclasses import replace
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from dftkmc.catalog import make_catalog
from dftkmc.cli import main
from dftkmc.coupled import CoupledConfig
from dftkmc.coupled_workflow import run_coupled
from dftkmc.exchange import ExchangeConfig, ReservoirConfig
from dftkmc.geometry import Lattice
from dftkmc.transport import ElectronicNetwork


KB_EV_K = 8.617333262145e-5
REPEATS = (2, 2, 3)
PULSES = [(0.25, 3e-10), (-0.20, 3e-10), (0.10, 3e-10)]


def configuration():
    pool = ReservoirConfig(capacity=8, initial_oxygen=3, mu0_eV=0.0,
                           kappa_eV=0.10, barrier_eV=0.65, attempt_Hz=1e13,
                           enabled=True)
    return CoupledConfig(
        repeats=REPEATS, temperature_K=1000.0, seed=31,
        bottom_potential_V=0.2, top_potential_V=0.45,
        exchange=ExchangeConfig(vacancy_formation_offset_eV=0.10, bottom=pool,
                                top=replace(pool, mu0_eV=0.02)),
    )


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_npz(path):
    with np.load(path) as data:
        return {key: data[key].copy() for key in data.files}


class CoupledWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = make_catalog(Lattice(repeats=REPEATS))
        cls.config = configuration()
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        cls.output = Path(temporary.name) / "trajectory"
        run_coupled(cls.catalog, cls.config, pulses=PULSES, n_vacancies=0,
                    max_events=1000, sample_every=1, output=cls.output,
                    plot=False, allow_illustrative=True)
        cls.summary = json.loads((cls.output / "summary.json").read_text())
        cls.events = read_csv(cls.output / "events.csv")
        cls.ledger = read_csv(cls.output / "oxygen_ledger.csv")
        cls.trace = read_csv(cls.output / "trace.csv")
        cls.states = read_npz(cls.output / "states.npz")

    def test_atom_and_vacancy_ledgers_replay_every_saved_state_and_final_pool(self):
        summary = self.summary
        self.assertTrue(summary["termination"]["target_reached"])
        self.assertEqual(summary["initial"]["vacancies"], 0)
        self.assertEqual({row["kind"] for row in self.events}, {"hop", "oxygen_out", "oxygen_in"})
        self.assertEqual(len(self.events), len(self.ledger))
        oxygen = np.asarray(summary["initial_oxygen_atom_ids"], dtype=int)
        vacancies = np.asarray(summary["initial_vacancy_ids"], dtype=int)
        pools = {name: list(ids) for name, ids in summary["initial_reservoir_oxygen_atom_ids"].items()}
        initial_atoms = set(oxygen[oxygen >= 0].tolist() + pools["bottom"] + pools["top"])
        born = set(vacancies[vacancies >= 0].tolist())
        retired = set()
        electrons = {"bottom": 0, "top": 0}
        replay = {0: (oxygen.copy(), vacancies.copy(), {name: len(ids) for name, ids in pools.items()})}

        for index, (event, ledger) in enumerate(zip(self.events, self.ledger), start=1):
            self.assertEqual(int(event["event_index"]), index)
            self.assertEqual(int(ledger["event_index"]), index)
            for key in ("time_s", "kind", "oxygen_id", "oxygen_from", "oxygen_to",
                        "reservoir", "vacancies_before", "vacancies_after",
                        "reservoir_oxygen_before", "reservoir_oxygen_after"):
                self.assertEqual(event[key], ledger[key])
            self.assertEqual(int(event["vacancies_before"]), int(np.count_nonzero(oxygen < 0)))
            self.assertEqual(json.loads(event["reservoir_oxygen_before"]),
                             {name: len(ids) for name, ids in pools.items()})
            kind, source, target = event["kind"], int(event["source"]), int(event["target"])
            oxygen_id, vacancy_id = int(event["oxygen_id"]), int(event["vacancy_id"])
            name = event["reservoir"]
            if kind == "hop":
                self.assertEqual(vacancies[source], vacancy_id)
                self.assertEqual(oxygen[target], oxygen_id)
                self.assertEqual(oxygen[source], -1)
                self.assertEqual(vacancies[target], -1)
                self.assertEqual(event["oxygen_from"], f"oxide:{target}")
                self.assertEqual(event["oxygen_to"], f"oxide:{source}")
                self.assertEqual(event["vacancy_identity_action"], "moved")
                oxygen[source], oxygen[target] = oxygen_id, -1
                vacancies[source], vacancies[target] = -1, vacancy_id
            elif kind == "oxygen_out":
                self.assertEqual(oxygen[source], oxygen_id)
                self.assertEqual(vacancies[source], -1)
                self.assertNotIn(vacancy_id, born)
                self.assertEqual(event["oxygen_from"], f"oxide:{source}")
                self.assertEqual(event["oxygen_to"], f"reservoir:{name}")
                self.assertEqual(event["vacancy_identity_action"], "born")
                born.add(vacancy_id)
                pools[name].append(oxygen_id)
                oxygen[source], vacancies[source] = -1, vacancy_id
                electrons[name] += self.config.charge
            else:
                self.assertEqual(kind, "oxygen_in")
                self.assertEqual(oxygen[source], -1)
                self.assertEqual(vacancies[source], vacancy_id)
                self.assertNotIn(vacancy_id, retired)
                self.assertEqual(pools[name].pop(), oxygen_id)
                self.assertEqual(event["oxygen_from"], f"reservoir:{name}")
                self.assertEqual(event["oxygen_to"], f"oxide:{source}")
                self.assertEqual(event["vacancy_identity_action"], "removed")
                retired.add(vacancy_id)
                oxygen[source], vacancies[source] = oxygen_id, -1
                electrons[name] -= self.config.charge
            all_atoms = oxygen[oxygen >= 0].tolist() + pools["bottom"] + pools["top"]
            self.assertEqual(set(all_atoms), initial_atoms)
            self.assertEqual(len(all_atoms), len(initial_atoms))
            self.assertEqual(int(event["vacancies_after"]), int(np.count_nonzero(oxygen < 0)))
            self.assertEqual(json.loads(event["reservoir_oxygen_after"]),
                             {name: len(ids) for name, ids in pools.items()})
            self.assertEqual(int(event["mass_residual"]), 0)
            self.assertEqual(int(event["reaction_charge_residual_e"]), 0)
            self.assertFalse(retired.intersection(vacancies[vacancies >= 0].tolist()))
            self.assertEqual(self.config.charge * int(np.count_nonzero(oxygen < 0)), sum(electrons.values()))
            replay[index] = (oxygen.copy(), vacancies.copy(), {name: len(ids) for name, ids in pools.items()})

        for state_index, event_count in enumerate(self.states["events"]):
            oxy, vac, inventory = replay[int(event_count)]
            np.testing.assert_array_equal(self.states["oxygen_ids"][state_index], oxy)
            np.testing.assert_array_equal(self.states["vacancy_ids"][state_index], vac)
            np.testing.assert_array_equal(self.states["vacancy"][state_index], oxy < 0)
            for name in ("bottom", "top"):
                self.assertEqual(int(self.trace[state_index][name + "_reservoir_oxygen"]), inventory[name])
        self.assertEqual(pools, summary["final_reservoir_oxygen_atom_ids"])
        for name in ("bottom", "top"):
            self.assertEqual(electrons[name], summary["reservoir_ledger"][name]["reaction_electrons"])
        self.assertEqual(summary["conservation"]["initial_total_oxygen"], len(initial_atoms))
        self.assertEqual(summary["conservation"]["final_total_oxygen"], len(initial_atoms))

    def test_saved_state_context_fields_and_event_before_potentials_are_aligned(self):
        state, summary = self.states, self.summary
        nstates, nsites = state["vacancy"].shape
        self.assertEqual(nstates, len(self.trace))
        self.assertEqual(nsites, summary["geometry"]["oxide_sites"])
        for key in ("oxygen_ids", "vacancy_ids", "potential_V"):
            self.assertEqual(state[key].shape, (nstates, nsites))
        self.assertEqual(state["edge_current_A"].shape, state["edge_field_V_nm"].shape)
        self.assertEqual(state["edge_current_A"].shape[0], nstates)
        np.testing.assert_array_equal(state["pbc"], [True, True, False])
        digest = hashlib.sha256(json.dumps(json.loads((self.output / "catalog.json").read_text()),
                                           sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        self.assertEqual(summary["catalog_sha256"], digest)
        self.assertEqual(summary["config"], json.loads((self.output / "coupled_config.json").read_text()))
        self.assertEqual(summary["calibration_status"], "illustrative")
        self.assertEqual(summary["finite_bias_global_detailed_balance"], "not_assumed")
        self.assertFalse(summary["diffusion"]["available"])
        self.assertTrue(all((self.output / name).exists() for name in summary["files"]))
        self.assertNotIn("coupled_overview.png", summary["files"])

        pairs = state["electronic_pairs"]
        distances = np.linalg.norm(state["electronic_displacements_A"], axis=1)
        for index, row in enumerate(self.trace):
            self.assertEqual(float(row["time_s"]), state["time_s"][index])
            self.assertEqual(int(row["events"]), state["events"][index])
            self.assertEqual(int(row["vacancies"]), int(state["vacancy"][index].sum()))
            phi = state["potential_V"][index]
            voltage_drop = phi[pairs[:, 0]] - phi[pairs[:, 1]]
            expected_field = 10.0 * voltage_drop / distances
            np.testing.assert_allclose(state["edge_field_V_nm"][index], expected_field,
                                       rtol=1e-10, atol=2e-14)
            vv = state["vacancy"][index, pairs[:, 0]] & state["vacancy"][index, pairs[:, 1]]
            conductance = np.where(vv, self.config.transport.vacancy_conductance_S,
                                   self.config.transport.background_conductance_S)
            conductance *= self.config.transport.reference_distance_A / distances
            np.testing.assert_allclose(state["edge_current_A"][index], conductance * voltage_drop,
                                       rtol=1e-9, atol=2e-20)
            self.assertEqual(int(row["mass_residual"]), 0)
            self.assertEqual(int(row["reaction_charge_residual_e"]), 0)

        pulse_indices = np.asarray([int(row["pulse_index"]) for row in self.trace])
        for event in self.events:
            event_index, pulse_index = int(event["event_index"]), int(event["pulse_index"])
            before_indices = np.flatnonzero((state["events"] == event_index - 1)
                                           & (pulse_indices == pulse_index)
                                           & (state["time_s"] <= float(event["time_before_s"])))
            after_indices = np.flatnonzero((state["events"] == event_index)
                                          & (pulse_indices == pulse_index)
                                          & (state["time_s"] == float(event["time_s"])))
            self.assertGreater(len(before_indices), 0)
            self.assertEqual(len(after_indices), 1)
            before, after = int(before_indices[-1]), int(after_indices[0])
            phi_before, phi_after = state["potential_V"][before], state["potential_V"][after]
            source, target = int(event["source"]), int(event["target"])
            self.assertAlmostEqual(float(event["source_potential_V"]), phi_before[source], places=13)
            self.assertAlmostEqual(float(event["applied_voltage_before_V"]),
                                   float(self.trace[before]["applied_voltage_V"]), places=13)
            if event["kind"] == "hop":
                self.assertAlmostEqual(float(event["target_potential_V"]), phi_before[target], places=13)
                work = self.config.charge * (phi_before[target] - phi_before[source])
                self.assertTrue(event["channel_provenance"])
            else:
                electrode = self.config.bottom_potential_V
                if event["reservoir"] == "top":
                    electrode += float(event["applied_voltage_before_V"])
                self.assertAlmostEqual(float(event["electrode_potential_V"]), electrode, places=13)
                work = int(event["direction"]) * self.config.charge * (phi_before[source] - electrode)
            self.assertAlmostEqual(float(event["electrical_work_eV"]), work, places=12)
            self.assertAlmostEqual(float(event["delta_eV"]), float(event["chemical_delta_eV"]) + work, places=12)
            self.assertAlmostEqual(float(event["potential_change_max_V"]),
                                   float(np.max(np.abs(phi_after - phi_before))), places=12)
            expected_rate = float(event["attempt_Hz"]) * math.exp(
                -float(event["barrier_eV"]) / (KB_EV_K * self.config.temperature_K))
            self.assertTrue(math.isclose(float(event["rate_Hz"]), expected_rate, rel_tol=2e-12))
            self.assertEqual(event["catalog_sha256"], digest)
            self.assertEqual(int(event["seed"]), self.config.seed)
            self.assertEqual(float(event["temperature_K"]), self.config.temperature_K)
            self.assertEqual(event["kinetic_field_closure"], summary["kinetic_field_closure"])
            self.assertEqual(event["parameter_status"], "illustrative")
            self.assertIn("frozen-field", event["energy_increment_interpretation"])

        for index, (bias, duration) in enumerate(PULSES):
            rows = np.flatnonzero(pulse_indices == index)
            self.assertEqual(state["time_s"][rows[0]], math.fsum(t for _, t in PULSES[:index]))
            self.assertEqual(state["time_s"][rows[-1]], math.fsum(t for _, t in PULSES[:index + 1]))
            self.assertAlmostEqual(float(self.trace[rows[0]]["requested_voltage_V"]), bias, places=14)
        final = read_npz(self.output / "final_transport.npz")
        np.testing.assert_allclose(final["potential_V"], state["potential_V"][-1], rtol=0, atol=0)
        np.testing.assert_allclose(final["edge_current_A"], state["edge_current_A"][-1], rtol=0, atol=0)
        for container in (summary["initial"], summary["final"], self.trace[0], state):
            self.assertFalse(any("msd" in key.lower() for key in container))
            for key in ("energy_eV", "total_energy_eV", "total_energy_relative_bottom_eV",
                        "electrostatic_energy_relative_bottom_eV", "diffusion_tensor_m2_s"):
                self.assertNotIn(key, container)

    def test_explicit_illustrative_permission_and_timing_validation_precede_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "denied"
            with self.assertRaisesRegex(ValueError, "illustrative"):
                run_coupled(self.catalog, self.config, duration_s=1e-10, output=output,
                            allow_illustrative=False)
            self.assertFalse(output.exists())
            for timing in ({}, {"duration_s": 1e-10, "pulses": PULSES}):
                with self.assertRaisesRegex(ValueError, "exactly one"):
                    run_coupled(self.catalog, self.config, output=output,
                                allow_illustrative=True, **timing)
                self.assertFalse(output.exists())

    def test_event_limit_and_model_domain_failure_save_incomplete_results(self):
        invalid_domain = replace(self.config,
                                 exchange=replace(self.config.exchange, vacancy_formation_offset_eV=4.0))
        with tempfile.TemporaryDirectory() as temporary:
            for name, config, limit, expected, expected_events in (
                    ("cap", self.config, 3, "event_limit_before_waveform_end", 3),
                    ("domain", invalid_domain, 100, "model_domain_error", 0)):
                with self.subTest(case=name):
                    output = Path(temporary) / name
                    result = run_coupled(self.catalog, config, duration_s=1.0, n_vacancies=0,
                                         max_events=limit, sample_every=2, output=output,
                                         plot=False, allow_illustrative=True)
                    saved = json.loads((output / "summary.json").read_text())
                    termination = saved["termination"]
                    self.assertEqual(result["termination"], termination)
                    self.assertEqual(termination["status"], expected)
                    self.assertFalse(termination["target_reached"])
                    self.assertTrue(termination["censored"])
                    self.assertEqual(termination["requested_end_time_s"], 1.0)
                    self.assertLess(termination["actual_time_s"], 1.0)
                    self.assertEqual(termination["event_count"], expected_events)
                    self.assertEqual(len(read_csv(output / "events.csv")), expected_events)
                    self.assertEqual(len(read_csv(output / "oxygen_ledger.csv")), expected_events)
                    state = read_npz(output / "states.npz")
                    self.assertEqual(state["time_s"][-1], termination["actual_time_s"])
                    self.assertEqual(int(state["events"][-1]), expected_events)
                    self.assertEqual(saved["final"]["mass_residual"], 0)
                    self.assertEqual(saved["final"]["reaction_charge_residual_e"], 0)
                    self.assertFalse(saved["diffusion"]["available"])
                    self.assertTrue(all((output / filename).exists() for filename in saved["files"]))
                    if name == "domain":
                        self.assertIn("negative exchange barrier", termination["error"])
                        self.assertEqual(termination["actual_time_s"], 0.0)
                    else:
                        self.assertIsNone(termination["error"])

    def test_cli_template_permission_pulses_and_seed_override_reproduce_python(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path, catalog_path = root / "config.json", root / "catalog.json"
            self.catalog.write(catalog_path)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                rc = main(["coupled-template", "--output", str(config_path), "--repeats", "2", "2", "3",
                           "--bottom-voltage", "0.1", "--top-voltage", "0.3"])
                self.assertEqual(rc, 0)
                config = CoupledConfig.load(config_path)
                self.assertEqual(config.bottom_potential_V, 0.1)
                self.assertEqual(config.top_potential_V, 0.3)
                base = ["device-run", "--catalog", str(catalog_path), "--config", str(config_path),
                        "--nvac", "0", "--seed", "47"]
                denied = root / "denied"
                rc = main(base + ["--duration", "1e-10", "--output", str(denied)])
                self.assertEqual(rc, 2)
                self.assertFalse(denied.exists())
                with self.assertRaises(SystemExit) as rejected:
                    main(base + ["--duration", "1e-10", "--pulse", "0.15", "1e-10",
                                 "--allow-illustrative", "--output", str(denied)])
                self.assertEqual(rejected.exception.code, 2)
                self.assertFalse(denied.exists())
                output = root / "cli"
                rc = main(base + ["--pulse", "0.15", "1e-10", "--pulse", "-0.1", "2e-10",
                                 "--allow-illustrative", "--output", str(output)])
                self.assertEqual(rc, 0)
            cli = json.loads((output / "summary.json").read_text())
            python = run_coupled(self.catalog, replace(config, seed=47), n_vacancies=0,
                                 pulses=[(0.15, 1e-10), (-0.1, 2e-10)], allow_illustrative=True)
            self.assertEqual(cli["config"]["seed"], 47)
            self.assertEqual(cli["config"], python["config"])
            self.assertEqual(cli["termination"], python["termination"])
            self.assertEqual(cli["initial"], python["initial"])
            self.assertEqual(cli["final"], python["final"])
            self.assertEqual(cli["event_counts"], python["event_counts"])

    def test_failure_after_one_accepted_event_exports_only_that_valid_state(self):
        original_solve = ElectronicNetwork.solve
        accepted_solution = []

        def fail_after_first_exchange(network, vacancy, bottom_V=None, top_V=None):
            if accepted_solution:
                raise RuntimeError("injected solve after first accepted event")
            solution = original_solve(network, vacancy, bottom_V=bottom_V, top_V=top_V)
            if np.any(vacancy):  # Initial Nv=0; the first event must remove oxygen.
                accepted_solution.append(solution)
            return solution

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "failed_after_event"
            with patch.object(ElectronicNetwork, "solve", new=fail_after_first_exchange):
                run_coupled(self.catalog, self.config, duration_s=1.0, n_vacancies=0,
                            max_events=100, sample_every=10, output=output,
                            allow_illustrative=True)
            summary = json.loads((output / "summary.json").read_text())
            events = read_csv(output / "events.csv")
            ledger = read_csv(output / "oxygen_ledger.csv")
            trace = read_csv(output / "trace.csv")
            state = read_npz(output / "states.npz")
            term = summary["termination"]
            self.assertEqual(term["status"], "model_domain_error")
            self.assertIn("injected solve after first accepted event", term["error"])
            self.assertFalse(term["target_reached"])
            self.assertTrue(term["censored"])
            self.assertEqual(term["event_count"], 1)
            self.assertEqual(len(events), 1)
            self.assertEqual(len(ledger), 1)
            event = events[0]
            self.assertEqual(event["kind"], "oxygen_out")
            self.assertEqual(term["actual_time_s"], float(event["time_s"]))
            self.assertEqual(float(trace[-1]["time_s"]), term["actual_time_s"])
            self.assertEqual(state["time_s"][-1], term["actual_time_s"])
            np.testing.assert_array_equal(state["events"], [0, 1])
            self.assertEqual(len(trace), len(state["time_s"]))
            oxygen = np.asarray(summary["initial_oxygen_atom_ids"], dtype=int)
            vacancies = np.asarray(summary["initial_vacancy_ids"], dtype=int)
            source, atom = int(event["source"]), int(event["oxygen_id"])
            self.assertEqual(oxygen[source], atom)
            oxygen[source], vacancies[source] = -1, int(event["vacancy_id"])
            np.testing.assert_array_equal(state["oxygen_ids"][-1], oxygen)
            np.testing.assert_array_equal(state["vacancy_ids"][-1], vacancies)
            np.testing.assert_array_equal(state["vacancy"][-1], oxygen < 0)
            pools = {name: list(ids) for name, ids in summary["initial_reservoir_oxygen_atom_ids"].items()}
            pools[event["reservoir"]].append(atom)
            self.assertEqual(pools, summary["final_reservoir_oxygen_atom_ids"])
            self.assertEqual(summary["event_counts"], {"oxygen_out": 1})
            solution = accepted_solution[0]
            np.testing.assert_array_equal(state["potential_V"][-1], solution.potential_V)
            np.testing.assert_array_equal(state["edge_current_A"][-1], solution.edge_current_A)
            np.testing.assert_array_equal(state["edge_field_V_nm"][-1], solution.edge_field_V_nm)
            final = read_npz(output / "final_transport.npz")
            np.testing.assert_array_equal(final["potential_V"], solution.potential_V)
            self.assertEqual(summary["final"]["equivalent_conductance_S"], solution.equivalent_conductance_S)
            self.assertEqual(summary["conservation"]["oxygen_balance_residual"], 0)
            self.assertEqual(summary["conservation"]["reaction_charge_residual_e"], 0)

    def test_failure_on_bias_switch_saves_prior_valid_field_and_completed_segment(self):
        config = replace(self.config, exchange=replace(
            self.config.exchange,
            bottom=replace(self.config.exchange.bottom, enabled=False),
            top=replace(self.config.exchange.top, enabled=False)))
        original_solve = ElectronicNetwork.solve
        good_solutions = []

        def fail_on_reverse_bias(network, vacancy, bottom_V=None, top_V=None):
            if top_V < bottom_V:
                raise RuntimeError("injected bias-switch solve")
            solution = original_solve(network, vacancy, bottom_V=bottom_V, top_V=top_V)
            good_solutions.append(solution)
            return solution

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "failed_bias"
            with patch.object(ElectronicNetwork, "solve", new=fail_on_reverse_bias):
                run_coupled(self.catalog, config, pulses=[(0.25, 2e-8), (-0.4, 3e-8)],
                            n_vacancies=0, max_events=100, sample_every=1,
                            output=output, allow_illustrative=True)
            summary = json.loads((output / "summary.json").read_text())
            trace = read_csv(output / "trace.csv")
            state = read_npz(output / "states.npz")
            term = summary["termination"]
            self.assertEqual(term["status"], "model_domain_error")
            self.assertIn("injected bias-switch solve", term["error"])
            self.assertFalse(term["target_reached"])
            self.assertTrue(term["censored"])
            self.assertEqual(term["actual_time_s"], 2e-8)
            self.assertEqual(term["requested_end_time_s"], math.fsum((2e-8, 3e-8)))
            self.assertEqual(term["event_count"], 0)
            self.assertEqual(read_csv(output / "events.csv"), [])
            self.assertEqual(read_csv(output / "oxygen_ledger.csv"), [])
            self.assertEqual(len(trace), len(state["time_s"]))
            self.assertEqual(state["time_s"][-1], 2e-8)
            self.assertTrue(all(int(row["pulse_index"]) == 0 for row in trace))
            self.assertTrue(all(math.isclose(float(row["requested_voltage_V"]), 0.25) for row in trace))
            self.assertAlmostEqual(summary["final"]["requested_voltage_V"], 0.25, places=14)
            self.assertAlmostEqual(summary["final"]["applied_voltage_V"], 0.25, places=14)
            self.assertEqual(summary["final"]["pulse_index"], 0)
            self.assertEqual(summary["initial_oxygen_atom_ids"], state["oxygen_ids"][-1].tolist())
            self.assertEqual(summary["initial_reservoir_oxygen_atom_ids"], summary["final_reservoir_oxygen_atom_ids"])
            self.assertEqual(summary["event_counts"], {})
            self.assertGreaterEqual(len(good_solutions), 1)
            np.testing.assert_array_equal(state["potential_V"][-1], good_solutions[-1].potential_V)
            np.testing.assert_array_equal(state["edge_current_A"][-1], good_solutions[-1].edge_current_A)
            final = read_npz(output / "final_transport.npz")
            np.testing.assert_array_equal(final["potential_V"], good_solutions[-1].potential_V)
            self.assertEqual(summary["conservation"]["oxygen_balance_residual"], 0)
            self.assertEqual(summary["conservation"]["reaction_charge_residual_e"], 0)


if __name__ == "__main__":
    unittest.main()
