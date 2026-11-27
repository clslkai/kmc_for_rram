"""Independent physical invariants for the DFT-ready fixed-charge KMC engine.

All fixtures below are deliberately synthetic. They check bookkeeping, transition
state consistency, geometry and stochastic kinetics; they do not validate HfO2
parameters, a DFT method, electron transport or experimental device switching.
"""
from copy import deepcopy
import math
import unittest

import numpy as np

from dftkmc.catalog import Catalog, make_catalog
from dftkmc.geometry import Lattice
from dftkmc.model import SimulationConfig, Simulator


KB_EV_K = 8.617333262145e-5
REPEATS = (2, 2, 2)


def synthetic_catalog(*, interaction_cutoff=3.6, site_offset=0.13,
                      pair_energy=-0.04, gamma=0.017, reference=None,
                      migration_cutoff=3.2):
    """Use a shared TS reference, independently chosen site and pair energies."""
    lattice = Lattice(repeats=REPEATS, boundary="bulk",
                      migration_cutoff_A=migration_cutoff,
                      interaction_cutoff_A=interaction_cutoff,
                      reference=reference)
    raw = deepcopy(make_catalog(lattice, status="illustrative", charges=(0, 2)).data)
    kinds = raw["structure"]["oxygen_kinds"]
    for q in ("0", "2"):
        raw["site_energies_eV"][q] = [0.0 if k == "O3" else site_offset for k in kinds]
        model = raw["charge_models"][q]
        model.update(pair_cutoff_A=interaction_cutoff,
                     pair_energies_eV={"O3-O3": pair_energy,
                                       "O3-O4": pair_energy,
                                       "O4-O4": pair_energy},
                     environment_mode="symmetric_pair",
                     environment_saddle_eV_per_neighbor=gamma,
                     provenance="Synthetic validation fixture, not DFT or literature data.")
        for channel in raw["channels"]:
            channel["data"][q].update(saddle_eV=0.65, attempt_Hz=1e9,
                                       ts_fraction=0.37, status="illustrative",
                                       provenance="Synthetic shared-saddle test fixture.")
    catalog = Catalog(raw)
    catalog.validate(allow_illustrative=True)
    return lattice, catalog


def new_sim(catalog, *, q=2, indices=(0,), field=(0.0, 0.0, 0.0), seed=17):
    sim = Simulator(catalog, SimulationConfig(repeats=REPEATS, charge=q,
                    temperature_K=600.0, seed=seed, boundary="bulk",
                    field_V_nm=field))
    sim.initialize(n_vacancies=len(indices), vacancy_indices=list(indices))
    return sim


def reverse_event(sim, event):
    return next(e for e in sim.events()
                if e.source == event.target and e.target == event.source
                and e.channel_id == event.channel_id)


def independent_chemical_energy(lattice, catalog, q, vacancy):
    """Explicit endpoint Hamiltonian, independent of the engine's delta formula."""
    site = np.asarray(catalog.data["site_energies_eV"][str(q)], dtype=float)
    energy = float(site[lattice.basis_id[np.flatnonzero(vacancy)]].sum())
    pairs = catalog.data["charge_models"][str(q)]["pair_energies_eV"]
    for i, j in lattice.interaction_pairs:
        if vacancy[i] and vacancy[j]:
            key = "-".join(sorted((str(lattice.kind[i]), str(lattice.kind[j]))))
            energy += pairs[key]
    return energy


class TransitionStateValidation(unittest.TestCase):
    def test_neb_barrier_does_not_double_count_endpoint_difference(self):
        lattice, catalog = synthetic_catalog(interaction_cutoff=1.0,
                                               site_offset=0.20, pair_energy=0.0, gamma=0.0)
        edge = next(e for e in lattice.migration_edges
                    if lattice.kind[e.source] == "O3" and lattice.kind[e.target] == "O4")
        sim = new_sim(catalog, q=0, indices=(edge.source,))
        event = next(e for e in sim.events()
                     if e.target == edge.target and e.channel_id == edge.channel_id)
        self.assertAlmostEqual(event.delta_eV, 0.20, places=12)
        self.assertAlmostEqual(event.barrier_eV, 0.65, places=12)
        sim.apply(event)
        reverse = reverse_event(sim, event)
        self.assertAlmostEqual(reverse.barrier_eV, 0.45, places=12)
        self.assertAlmostEqual(event.barrier_eV - reverse.barrier_eV,
                               event.delta_eV, places=12)

    def test_event_delta_matches_total_energy_with_interacting_vacancies(self):
        lattice, catalog = synthetic_catalog()
        i, j = (int(x) for x in lattice.interaction_pairs[0])
        sim = new_sim(catalog, indices=(i, j))
        for _ in range(12):
            choices = sim.events()
            self.assertTrue(choices)
            event = choices[(_ * 7) % len(choices)]
            initial = independent_chemical_energy(lattice, catalog, 2, sim.vacancy)
            self.assertAlmostEqual(sim.snapshot()["energy_eV"], initial, places=11)
            sim.apply(event)
            final = independent_chemical_energy(lattice, catalog, 2, sim.vacancy)
            self.assertAlmostEqual(sim.snapshot()["energy_eV"], final, places=11)
            self.assertAlmostEqual(event.delta_eV, final - initial, places=11)
            reverse = reverse_event(sim, event)
            self.assertAlmostEqual(event.delta_eV, -reverse.delta_eV, places=11)
            self.assertAlmostEqual(event.barrier_eV - reverse.barrier_eV,
                                   event.delta_eV, places=11)
            self.assertAlmostEqual(math.log(event.rate_Hz / reverse.rate_Hz),
                                   -event.delta_eV / (KB_EV_K * 600.0), places=10)

    def test_interaction_graph_can_exclude_a_migration_pair(self):
        lattice, catalog = synthetic_catalog(interaction_cutoff=1.0,
                                               pair_energy=-0.12, gamma=0.03)
        self.assertEqual(len(lattice.interaction_pairs), 0)
        edge = lattice.migration_edges[0]
        sim = new_sim(catalog, indices=(edge.source, edge.target))
        snapshot = sim.snapshot()
        self.assertEqual(snapshot["pair_count"], 0)
        self.assertEqual(snapshot["largest_cluster"], 1)
        expected = independent_chemical_energy(lattice, catalog, 2, sim.vacancy)
        self.assertAlmostEqual(snapshot["energy_eV"], expected, places=12)
        site = catalog.site_energies(2)
        for event in sim.events():
            expected_delta = site[lattice.basis_id[event.target]] - site[lattice.basis_id[event.source]]
            self.assertAlmostEqual(event.delta_eV, expected_delta, places=12)
            self.assertEqual(event.env_n3 + event.env_n4, 0)

    def test_nonmigration_pair_can_contribute_to_interaction_energy(self):
        lattice, catalog = synthetic_catalog(interaction_cutoff=3.6,
                                               site_offset=0.0, pair_energy=-0.041, gamma=0.0)
        migration = {frozenset((e.source, e.target)) for e in lattice.migration_edges}
        pair = next((int(i), int(j)) for i, j in lattice.interaction_pairs
                    if frozenset((int(i), int(j))) not in migration)
        sim = new_sim(catalog, q=0, indices=pair)
        self.assertAlmostEqual(sim.snapshot()["energy_eV"], -0.041, places=12)
        self.assertEqual(sim.snapshot()["pair_count"], 1)
        self.assertEqual(sim.snapshot()["largest_cluster"], 2)

    def test_missing_channel_fails_closed(self):
        _, catalog = synthetic_catalog()
        raw = deepcopy(catalog.data)
        self.assertTrue(raw["channels"])
        raw["channels"].pop()
        with self.assertRaises(ValueError):
            incomplete = Catalog(raw)
            incomplete.validate(allow_illustrative=True)
            sim = new_sim(incomplete)
            sim.events()

    def test_dilute_only_rejects_an_unmodelled_nearby_vacancy(self):
        lattice, catalog = synthetic_catalog()
        raw = deepcopy(catalog.data)
        raw["charge_models"]["2"]["environment_mode"] = "dilute_only"
        raw["charge_models"]["2"]["pair_energies_eV"] = dict.fromkeys(("O3-O3", "O3-O4", "O4-O4"), 0.0)
        raw["charge_models"]["2"]["environment_saddle_eV_per_neighbor"] = 0.0
        i, j = (int(x) for x in lattice.interaction_pairs[0])
        with self.assertRaises(ValueError):
            sim = new_sim(Catalog(raw), indices=(i, j))
            sim.events()


class GeometryAndFieldValidation(unittest.TestCase):
    def test_bulk_displacements_use_the_full_monoclinic_cell(self):
        lattice, _ = synthetic_catalog()
        inverse = np.linalg.inv(lattice.cell_A)
        self.assertGreater(abs(lattice.cell_A[2, 0]), 1e-2)
        wrapped_c_crossings = 0
        for edge in lattice.migration_edges:
            displacement = np.asarray(edge.displacement_A)
            raw = lattice.oxygen_A[edge.target] - lattice.oxygen_A[edge.source]
            image = (raw - displacement) @ inverse
            np.testing.assert_allclose(image, np.rint(image), atol=1e-10)
            self.assertLessEqual(np.linalg.norm(displacement), 3.2 + 1e-10)
            if abs(image[2]) > 0.5:
                wrapped_c_crossings += 1
                correction = np.rint(image) @ lattice.cell_A
                np.testing.assert_allclose(raw - correction, displacement, atol=1e-10)
                self.assertGreater(abs(correction[0]), 1e-3)
        self.assertGreater(wrapped_c_crossings, 0)

    def test_a_periodic_hop_preserves_label_and_unwrapped_displacement(self):
        lattice, catalog = synthetic_catalog()
        edge = next(e for e in lattice.migration_edges
                    if np.linalg.norm(lattice.oxygen_A[e.target] - lattice.oxygen_A[e.source]
                                      - np.asarray(e.displacement_A)) > 1.0)
        sim = new_sim(catalog, indices=(edge.source,))
        label = int(sim.labels[edge.source])
        event = next(e for e in sim.events()
                     if e.target == edge.target and e.channel_id == edge.channel_id)
        sim.apply(event)
        self.assertEqual(int(sim.vacancy.sum()), 1)
        self.assertEqual(int(sim.labels[event.target]), label)
        np.testing.assert_allclose(sim.displacements_A[label], event.displacement_A, atol=1e-12)
        reverse = reverse_event(sim, event)
        sim.apply(reverse)
        self.assertEqual(int(sim.labels[edge.source]), label)
        np.testing.assert_allclose(sim.displacements_A[label], 0.0, atol=1e-12)

    def test_two_vacancy_labels_and_mass_survive_many_events(self):
        lattice, catalog = synthetic_catalog()
        sim = new_sim(catalog, indices=(0, lattice.n - 1))
        expected_displacement = np.zeros((2, 3))
        for count in range(60):
            choices = sim.events()
            event = choices[(11 * count + 3) % len(choices)]
            label = int(sim.labels[event.source])
            expected_displacement[label] += np.asarray(event.displacement_A)
            sim.apply(event)
            self.assertEqual(int(sim.vacancy.sum()), 2)
            np.testing.assert_array_equal(np.sort(sim.labels[sim.vacancy]), [0, 1])
            np.testing.assert_allclose(sim.displacements_A, expected_displacement, atol=1e-11)

    def test_positive_charge_field_work_and_neutral_invariance(self):
        lattice, catalog = synthetic_catalog(site_offset=0.0, pair_energy=0.0, gamma=0.0)
        edge = lattice.migration_edges[0]
        displacement = np.asarray(edge.displacement_A)
        field_V_nm = 0.1 * displacement / np.linalg.norm(displacement)
        sim = new_sim(catalog, indices=(edge.source,), field=tuple(field_V_nm))
        event = next(e for e in sim.events()
                     if e.target == edge.target and e.channel_id == edge.channel_id)
        work_eV = 2.0 * float(field_V_nm @ displacement) / 10.0
        self.assertAlmostEqual(event.delta_eV, -work_eV, places=12)
        sim.apply(event)
        reverse = reverse_event(sim, event)
        self.assertGreater(event.rate_Hz, reverse.rate_Hz)
        self.assertAlmostEqual(math.log(event.rate_Hz / reverse.rate_Hz),
                               work_eV / (KB_EV_K * 600.0), places=10)
        neutral = new_sim(catalog, q=0, indices=(edge.source,), field=tuple(field_V_nm))
        neutral_zero = new_sim(catalog, q=0, indices=(edge.source,))
        a = {(e.source, e.target, e.channel_id): e.rate_Hz for e in neutral.events()}
        b = {(e.source, e.target, e.channel_id): e.rate_Hz for e in neutral_zero.events()}
        self.assertEqual(a, b)

    def test_out_of_domain_field_does_not_clip_negative_barriers(self):
        lattice, catalog = synthetic_catalog(site_offset=0.0, pair_energy=0.0, gamma=0.0)
        edge = lattice.migration_edges[0]
        direction = np.asarray(edge.displacement_A)
        field = tuple(100.0 * direction / np.linalg.norm(direction))
        with self.assertRaises(ValueError):
            sim = new_sim(catalog, indices=(edge.source,), field=field)
            sim.events()


class TwoStateKineticsValidation(unittest.TestCase):
    def test_continuous_time_equilibrium_weights_differ_from_jump_counts(self):
        # Eight synthetic oxygen basis sites; only basis 0 and 1 share a migration
        # edge. All other sites and translated dimers are outside the cutoff.
        # This fixture is a mathematical two-state process, not a hafnia crystal.
        reference = {
            "cell_A": [[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]],
            "hf_fractional": [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
                              [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]],
            "oxygen_fractional": [[0.10, 0.10, 0.10], [0.20, 0.10, 0.10],
                                  [0.55, 0.10, 0.10], [0.10, 0.55, 0.10],
                                  [0.10, 0.10, 0.55], [0.55, 0.55, 0.10],
                                  [0.55, 0.10, 0.55], [0.10, 0.55, 0.55]],
            "oxygen_kinds": ["O3", "O4", "O3", "O4", "O3", "O4", "O3", "O4"],
            "provenance": "Synthetic disconnected dimers for an analytic CTMC validation only.",
        }
        lattice, catalog = synthetic_catalog(reference=reference,
                                               interaction_cutoff=1.0,
                                               migration_cutoff=3.2,
                                               site_offset=0.08, pair_energy=0.0, gamma=0.0)
        raw = deepcopy(catalog.data)
        for channel in raw["channels"]:
            for q in ("0", "2"):
                channel["data"][q].update(saddle_eV=0.20, attempt_Hz=1e5)
        catalog = Catalog(raw)
        index_a = int(np.flatnonzero((lattice.basis_id == 0)
                                    & np.all(lattice.cell_index == (0, 0, 0), axis=1))[0])
        index_b = int(np.flatnonzero((lattice.basis_id == 1)
                                    & np.all(lattice.cell_index == (0, 0, 0), axis=1))[0])
        sim = new_sim(catalog, q=0, indices=(index_a,), seed=202603)
        a_to_b = 1e5 * math.exp(-0.20 / (KB_EV_K * 600.0))
        b_to_a = 1e5 * math.exp(-0.12 / (KB_EV_K * 600.0))
        dwell_a, dwell_b = [], []
        for _ in range(4000):
            current = int(np.flatnonzero(sim.vacancy)[0])
            self.assertIn(current, (index_a, index_b))
            self.assertEqual(len(sim.events()), 1)
            before = sim.time_s
            self.assertIsNotNone(sim.step())
            (dwell_a if current == index_a else dwell_b).append(sim.time_s - before)
        # Event counts alternate 1:1, but physical occupancy is residence-time
        # weighted. This catches using the jump chain as an equilibrium sample.
        self.assertEqual(len(dwell_a), len(dwell_b))
        self.assertAlmostEqual(np.mean(dwell_a) * a_to_b, 1.0, delta=0.12)
        self.assertAlmostEqual(np.mean(dwell_b) * b_to_a, 1.0, delta=0.12)
        measured_b = sum(dwell_b) / (sum(dwell_a) + sum(dwell_b))
        equilibrium_b = a_to_b / (a_to_b + b_to_a)
        self.assertAlmostEqual(measured_b, equilibrium_b, delta=0.035)


if __name__ == "__main__":
    unittest.main()
