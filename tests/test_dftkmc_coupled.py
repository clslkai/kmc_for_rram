"""Independent integration acceptance for the illustrative coupled device.

These tests use atom inventories, state-function differences, gauge changes,
and physical stopping times as oracles.  Finite-bias detailed balance is tested
at one explicitly frozen potential, never across two refreshed network fields.
They do not establish real-material forming or electronic transport accuracy.
"""
from copy import deepcopy
from dataclasses import replace
import math
import unittest
from unittest.mock import patch
import warnings

import numpy as np

from dftkmc.catalog import Catalog, make_catalog
from dftkmc.coupled import CoupledConfig, CoupledDeviceSimulator
from dftkmc.exchange import ExchangeConfig, ReservoirConfig
from dftkmc.geometry import Lattice
from dftkmc.transport import TransportConfig


KB_EV_K = 8.617333262145e-5
REPEATS = (2, 2, 3)


def configuration(**changes):
    bottom = ReservoirConfig(capacity=6, initial_oxygen=1, mu0_eV=0.04,
                             kappa_eV=0.10, barrier_eV=0.95, enabled=True)
    top = replace(bottom, mu0_eV=-0.03)
    base = CoupledConfig(
        repeats=REPEATS, charge=2, temperature_K=1000.0, seed=17,
        bottom_potential_V=0.0, top_potential_V=0.35,
        exchange=ExchangeConfig(vacancy_formation_offset_eV=0.32,
                                bottom=bottom, top=top),
    )
    return replace(base, **changes)


def event_key(event):
    return (event.kind, event.source, event.target, event.reservoir,
            event.channel_id)


def reverse_of(event, events):
    if event.kind == "hop":
        return next(candidate for candidate in events
                    if candidate.kind == "hop" and candidate.source == event.target
                    and candidate.target == event.source
                    and candidate.channel_id == event.channel_id)
    reverse_kind = "oxygen_in" if event.kind == "oxygen_out" else "oxygen_out"
    return next(candidate for candidate in events
                if candidate.kind == reverse_kind and candidate.source == event.source
                and candidate.reservoir == event.reservoir)


def transaction_state(sim):
    """Save the accepted state, including identities and stochastic progress."""
    return {
        "vacancy": sim.vacancy.copy(),
        "oxygen_ids": sim.oxygen_ids.copy(),
        "vacancy_ids": sim.vacancy_ids.copy(),
        "pools": {name: (pool.count, pool.out_count, pool.in_count,
                          list(sim.reservoir_oxygen_ids[name]))
                  for name, pool in sim.reservoirs.items()},
        "electrons": dict(sim.electrons_to_electrode),
        "progress": (sim.time_s, sim.event_count, sim.state_epoch,
                     sim.requested_voltage_V, dict(sim.counts)),
        "identity_history": (sim._next_vacancy_id, set(sim._retired_vacancy_ids)),
        "formation_history": (deepcopy(sim.formation_observations),
                              list(sim._criterion_durations), sim._criterion_start_s,
                              sim._last_criterion),
        "rng": deepcopy(sim.rng.bit_generator.state),
        "snapshot": deepcopy(sim.snapshot()),
        "fields_object": sim.fields,
        "fields_data": sim.fields.to_dict(),
    }


class CoupledDeviceIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = make_catalog(Lattice(repeats=REPEATS))

    def simulator(self, config=None, indices=()):
        sim = CoupledDeviceSimulator(self.catalog, config or configuration())
        sim.initialize(vacancy_indices=list(indices))
        return sim

    def assert_inventory(self, sim):
        """Count the physical atom labels independently of snapshot residuals."""
        oxide = sim.oxygen_ids[~sim.vacancy].tolist()
        pools = [atom for name in ("bottom", "top")
                 for atom in sim.reservoir_oxygen_ids[name]]
        self.assertEqual(len(oxide) + len(pools), sim.initial_total_oxygen)
        self.assertEqual(len(set(oxide + pools)), sim.initial_total_oxygen)
        self.assertTrue(np.all(sim.oxygen_ids[sim.vacancy] == -1))
        self.assertTrue(np.all(sim.vacancy_ids[~sim.vacancy] == -1))
        for name in ("bottom", "top"):
            pool = sim.reservoirs[name]
            self.assertEqual(len(sim.reservoir_oxygen_ids[name]), pool.count)
            self.assertEqual(pool.count, pool.initial_count + pool.out_count - pool.in_count)
            self.assertEqual(sim.electrons_to_electrode[name],
                             sim.cfg.charge * (pool.out_count - pool.in_count))
        self.assertEqual(sim.cfg.charge * (int(sim.vacancy.sum()) - sim.initial_vacancies),
                         sum(sim.electrons_to_electrode.values()))
        snapshot = sim.snapshot()
        self.assertEqual(snapshot["mass_residual"], 0)
        self.assertEqual(snapshot["reaction_charge_residual_e"], 0)
        sim.check()

    def assert_detailed_balance(self, forward, reverse, delta, temperature):
        self.assertAlmostEqual(forward.delta_eV, delta, places=11)
        self.assertAlmostEqual(reverse.delta_eV, -delta, places=11)
        self.assertAlmostEqual(forward.barrier_eV - reverse.barrier_eV, delta, places=11)
        self.assertAlmostEqual(math.log(forward.rate_Hz / reverse.rate_Hz),
                               -delta / (KB_EV_K * temperature), places=10)

    def assert_transaction_unchanged(self, sim, before):
        after = transaction_state(sim)
        for key in ("vacancy", "oxygen_ids", "vacancy_ids"):
            np.testing.assert_array_equal(after[key], before[key])
        for key in ("pools", "electrons", "progress", "identity_history",
                    "formation_history", "rng", "snapshot", "fields_data"):
            self.assertEqual(after[key], before[key], key)
        self.assertIs(after["fields_object"], before["fields_object"])

    def test_empty_initial_vacancy_population_is_valid_without_a_tracer_claim(self):
        sim = CoupledDeviceSimulator(self.catalog, configuration())
        with self.assertRaisesRegex(ValueError, "initialize first"):
            sim.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            sim.initialize(n_vacancies=0)
            snapshot = sim.snapshot()
            events = sim.events()
        self.assertEqual(snapshot["vacancies"], 0)
        self.assertEqual(snapshot["oxide_oxygen"], sim.n)
        self.assertFalse(snapshot["vacancy_graph_spanning"])
        self.assertFalse(snapshot["model_formation_criterion"])
        self.assertTrue(math.isfinite(snapshot["chemical_free_energy_eV"]))
        self.assertGreater(snapshot["equivalent_conductance_S"], 0.0)
        self.assertTrue(events)
        self.assertEqual({event.kind for event in events}, {"oxygen_out"})
        self.assertFalse(any("msd" in key.lower() or "diffusion" in key.lower()
                             or "tracer" in key.lower() for key in snapshot))
        self.assert_inventory(sim)

    def test_interface_source_and_sink_move_real_oxygen_and_retire_vacancy_ids(self):
        config = configuration()
        exchange = replace(config.exchange,
                           bottom=replace(config.exchange.bottom, initial_oxygen=0),
                           top=replace(config.exchange.top, initial_oxygen=0))
        for name in ("bottom", "top"):
            with self.subTest(electrode=name):
                sim = self.simulator(replace(config, exchange=exchange))
                event = next(e for e in sim.events()
                             if e.kind == "oxygen_out" and e.reservoir == name)
                site, oxygen_id = event.source, int(sim.oxygen_ids[event.source])
                before_ids = set(sim.oxygen_ids.tolist())
                ledger = sim.apply(event)
                vacancy_id = int(sim.vacancy_ids[site])
                self.assertTrue(sim.vacancy[site])
                self.assertEqual(sim.reservoir_oxygen_ids[name], [oxygen_id])
                self.assertEqual(ledger["oxygen_id"], oxygen_id)
                self.assertEqual(ledger["oxygen_from"], f"oxide:{site}")
                self.assertEqual(ledger["oxygen_to"], f"reservoir:{name}")
                self.assertEqual(ledger["vacancy_identity_action"], "born")
                self.assertEqual(sim.electrons_to_electrode[name], config.charge)
                self.assert_inventory(sim)

                inward = reverse_of(event, sim.events())
                ledger = sim.apply(inward)
                self.assertFalse(sim.vacancy[site])
                self.assertEqual(int(sim.oxygen_ids[site]), oxygen_id)
                self.assertEqual(sim.reservoir_oxygen_ids[name], [])
                self.assertEqual(ledger["oxygen_from"], f"reservoir:{name}")
                self.assertEqual(ledger["oxygen_to"], f"oxide:{site}")
                self.assertEqual(ledger["vacancy_id"], vacancy_id)
                self.assertEqual(ledger["vacancy_identity_action"], "removed")
                self.assertEqual(set(sim.oxygen_ids.tolist()), before_ids)
                self.assertEqual(sim.electrons_to_electrode[name], 0)
                self.assert_inventory(sim)

                rebirth = reverse_of(inward, sim.events())
                sim.apply(rebirth)
                self.assertGreater(int(sim.vacancy_ids[site]), vacancy_id)
                self.assert_inventory(sim)

    def test_hop_transports_an_oxygen_atom_opposite_to_the_same_vacancy_id(self):
        sim = CoupledDeviceSimulator(self.catalog, configuration())
        edge = sim.edges[0]
        sim.initialize(vacancy_indices=[edge.source])
        event = next(e for e in sim.events() if e.kind == "hop"
                     and e.source == edge.source and e.target == edge.target
                     and e.channel_id == edge.channel_id)
        oxygen_id = int(sim.oxygen_ids[event.target])
        vacancy_id = int(sim.vacancy_ids[event.source])
        pools_before = {name: list(ids) for name, ids in sim.reservoir_oxygen_ids.items()}
        ledger = sim.apply(event)
        self.assertEqual(int(sim.oxygen_ids[event.source]), oxygen_id)
        self.assertEqual(int(sim.oxygen_ids[event.target]), -1)
        self.assertEqual(int(sim.vacancy_ids[event.target]), vacancy_id)
        self.assertEqual(int(sim.vacancy_ids[event.source]), -1)
        self.assertEqual(sim.reservoir_oxygen_ids, pools_before)
        self.assertEqual(ledger["oxygen_from"], f"oxide:{event.target}")
        self.assertEqual(ledger["oxygen_to"], f"oxide:{event.source}")
        self.assertEqual(ledger["vacancy_identity_action"], "moved")
        self.assert_inventory(sim)

    def test_finite_reservoir_has_no_out_events_when_full_or_in_events_when_empty(self):
        config = configuration()
        exchange = replace(config.exchange,
                           bottom=replace(config.exchange.bottom, capacity=1, initial_oxygen=0),
                           top=replace(config.exchange.top, capacity=1, initial_oxygen=0))
        sim = self.simulator(replace(config, exchange=exchange))
        self.assertFalse(any(e.kind == "oxygen_in" for e in sim.events()))
        outward = next(e for e in sim.events()
                       if e.kind == "oxygen_out" and e.reservoir == "top")
        sim.apply(outward)
        self.assertEqual(sim.reservoirs["top"].count, 1)
        self.assertFalse(any(e.kind == "oxygen_out" and e.reservoir == "top"
                             for e in sim.events()))
        self.assertTrue(any(e.kind == "oxygen_out" and e.reservoir == "bottom"
                            for e in sim.events()))
        sim.apply(reverse_of(outward, sim.events()))
        self.assertEqual(sim.reservoirs["top"].count, 0)
        self.assertFalse(any(e.kind == "oxygen_in" and e.reservoir == "top"
                             for e in sim.events()))
        self.assert_inventory(sim)

    def test_zero_bias_global_chemical_detailed_balance_for_hops_and_both_reactions(self):
        # A nonzero common voltage also checks that a voltage zero does not
        # enter chemical equilibrium. Pair interactions and pool entropy remain on.
        config = configuration(bottom_potential_V=0.7, top_potential_V=0.7)
        for kind in ("hop", "oxygen_out", "oxygen_in"):
            with self.subTest(kind=kind):
                sim = self.simulator(config, indices=[0, 1])
                event = next(e for e in sim.events() if e.kind == kind)
                before = sim.chemical_free_energy_eV()
                sim.apply(event)
                delta = sim.chemical_free_energy_eV() - before
                reverse = reverse_of(event, sim.events())  # Actual refreshed field is uniform.
                self.assertAlmostEqual(event.electrical_work_eV, 0.0, places=13)
                self.assertAlmostEqual(event.chemical_delta_eV, delta, places=11)
                self.assertAlmostEqual(event.oxide_delta_eV + event.reservoir_delta_eV,
                                       delta, places=11)
                self.assert_detailed_balance(event, reverse, delta, config.temperature_K)
                self.assert_inventory(sim)

    def test_frozen_field_detailed_balance_survives_a_separate_network_refresh(self):
        for kind in ("oxygen_out", "hop"):
            with self.subTest(kind=kind):
                sim = self.simulator(indices=[] if kind == "oxygen_out" else [0])
                actual_before = sim.fields
                phi = actual_before.potential_V.copy()
                electrodes = {"bottom": actual_before.bottom_potential_V,
                              "top": actual_before.top_potential_V}
                with self.assertRaisesRegex(ValueError, "both"):
                    sim.events(potential_V=phi)
                event = next(e for e in sim.events(phi, electrodes) if e.kind == kind)
                self.assertIs(sim.fields, actual_before)
                before = sim.frozen_field_state_function_eV(phi, electrodes)
                if kind == "hop":
                    work = sim.cfg.charge * (phi[event.target] - phi[event.source])
                else:
                    work = sim.cfg.charge * (phi[event.source] - electrodes[event.reservoir])
                self.assertAlmostEqual(event.electrical_work_eV, work, places=12)
                sim.apply(event)
                actual_after = sim.fields
                self.assertGreater(float(np.max(np.abs(actual_after.potential_V - phi))), 1e-6)
                after = sim.frozen_field_state_function_eV(phi, electrodes)
                reverse = reverse_of(event, sim.events(phi, electrodes))
                self.assertIs(sim.fields, actual_after)
                self.assert_detailed_balance(event, reverse, after - before,
                                             sim.cfg.temperature_K)
                self.assert_inventory(sim)

    def test_exchange_uses_actual_compliance_limited_electrode_voltage(self):
        config = configuration(top_potential_V=0.8,
                               transport=TransportConfig(current_compliance_A=1e-14))
        sim = self.simulator(config)
        self.assertTrue(sim.fields.compliance_active)
        self.assertLess(sim.fields.applied_voltage_V, 0.8)
        self.assertLessEqual(abs(sim.fields.terminal_current_A), 1e-14 * (1.0 + 1e-9))
        event = next(e for e in sim.events()
                     if e.kind == "oxygen_out" and e.reservoir == "top")
        work = config.charge * (sim.fields.potential_V[event.source] - sim.fields.top_potential_V)
        requested_work = config.charge * (sim.fields.potential_V[event.source] - config.top_potential_V)
        self.assertAlmostEqual(event.electrical_work_eV, work, places=13)
        self.assertAlmostEqual(event.electrode_potential_V, sim.fields.top_potential_V, places=13)
        self.assertGreater(abs(event.electrical_work_eV - requested_work), 0.1)

    def test_common_voltage_gauge_shift_preserves_the_seeded_event_trace(self):
        config = configuration(seed=23)
        shift = 2.125
        a = self.simulator(config, indices=[0, 1])
        b = self.simulator(replace(config,
                                  bottom_potential_V=config.bottom_potential_V + shift,
                                  top_potential_V=config.top_potential_V + shift), indices=[0, 1])
        for _ in range(12):
            left, right = a.step(), b.step()
            for key in ("kind", "source", "target", "reservoir", "channel_id", "oxygen_id", "vacancy_id"):
                self.assertEqual(left[key], right[key])
            for key in ("rate_Hz", "dt_s", "time_s"):
                self.assertTrue(math.isclose(left[key], right[key], rel_tol=2e-12, abs_tol=1e-18), key)
            for key in ("chemical_delta_eV", "electrical_work_eV", "delta_eV"):
                self.assertAlmostEqual(left[key], right[key], places=12)
            np.testing.assert_array_equal(a.vacancy, b.vacancy)
            np.testing.assert_array_equal(a.oxygen_ids, b.oxygen_ids)
            np.testing.assert_allclose(b.fields.potential_V - a.fields.potential_V,
                                       shift, rtol=0.0, atol=2e-12)
        self.assert_inventory(a)
        self.assert_inventory(b)

    def test_common_site_and_saddle_energy_shift_needs_no_formation_offset_change(self):
        shifted = Catalog(self.catalog.data)
        offset = 3.25
        for branch, energies in shifted.data["site_energies_eV"].items():
            shifted.data["site_energies_eV"][branch] = [value + offset for value in energies]
        for channel in shifted.data["channels"]:
            for data in channel["data"].values():
                data["saddle_eV"] += offset
        config = configuration()
        a = self.simulator(config, indices=[0, 1])
        b = CoupledDeviceSimulator(shifted, config)
        b.initialize(vacancy_indices=[0, 1])
        self.assertEqual(a.cfg.exchange.vacancy_formation_offset_eV,
                         b.cfg.exchange.vacancy_formation_offset_eV)
        for _ in range(2):
            left = {event_key(e): e for e in a.events()}
            right = {event_key(e): e for e in b.events()}
            self.assertEqual(set(left), set(right))
            self.assertEqual({e.kind for e in left.values()}, {"hop", "oxygen_out", "oxygen_in"})
            self.assertAlmostEqual(a.chemical_free_energy_eV(), b.chemical_free_energy_eV(), places=12)
            for key in left:
                self.assertAlmostEqual(left[key].delta_eV, right[key].delta_eV, places=12)
                self.assertAlmostEqual(left[key].barrier_eV, right[key].barrier_eV, places=12)
                self.assertTrue(math.isclose(left[key].rate_Hz, right[key].rate_Hz, rel_tol=2e-12))
            event = next(e for e in left.values() if e.kind == "oxygen_out")
            a.apply(event)
            b.apply(right[event_key(event)])
        self.assert_inventory(a)
        self.assert_inventory(b)

    def test_waveform_stops_exactly_at_each_boundary_with_active_events(self):
        sim = self.simulator(configuration(temperature_K=800.0), indices=[0])
        pulses = [(0.25, 3e-8), (-0.20, 4e-8), (0.10, 5e-8)]
        result = sim.run_waveform(pulses, max_events=2000, sample_every=3)
        self.assertEqual(result["status"], "completed_waveform")
        self.assertIsNone(result["error"])
        self.assertTrue(result["target_reached"])
        expected_end = math.fsum(duration for _, duration in pulses)
        self.assertEqual(sim.time_s, expected_end)
        self.assertEqual(result["requested_end_time_s"], expected_end)
        self.assertGreater(len(result["events"]), 0)
        for index, (bias, duration) in enumerate(pulses):
            start = math.fsum(value[1] for value in pulses[:index])
            end = math.fsum(value[1] for value in pulses[:index + 1])
            rows = [row for row in result["records"] if row["pulse_index"] == index]
            self.assertEqual(rows[0]["time_s"], start)
            self.assertEqual(rows[-1]["time_s"], end)
            self.assertEqual(rows[0]["requested_voltage_V"], bias)
            for event in result["events"]:
                if event["pulse_index"] == index:
                    self.assertGreaterEqual(event["time_before_s"], start)
                    self.assertGreater(event["time_s"], start)
                    self.assertLess(event["time_s"], end)
                    self.assertEqual(event["requested_voltage_V"], bias)
        self.assert_inventory(sim)

    def test_event_cap_is_censoring_while_no_event_intervals_complete(self):
        capped = self.simulator(indices=[0])
        result = capped.run(1.0, max_events=1, sample_every=1)
        self.assertEqual(result["status"], "event_limit_before_waveform_end")
        self.assertFalse(result["target_reached"])
        self.assertIsNone(result["error"])
        self.assertEqual(len(result["events"]), 1)
        self.assertLess(capped.time_s, result["requested_end_time_s"])
        self.assertEqual(result["records"][-1]["time_s"], capped.time_s)
        config = configuration()
        exchange = replace(config.exchange, bottom=replace(config.exchange.bottom, enabled=False),
                           top=replace(config.exchange.top, enabled=False))
        idle = self.simulator(replace(config, exchange=exchange))
        result = idle.run_waveform([(0.35, 2e-6), (-0.25, 3e-6)], max_events=1)
        self.assertEqual(result["status"], "completed_waveform")
        self.assertTrue(result["target_reached"])
        self.assertEqual(result["events"], [])
        self.assertEqual(idle.time_s, math.fsum((2e-6, 3e-6)))
        self.assertEqual(int(idle.vacancy.sum()), 0)
        self.assert_inventory(capped)
        self.assert_inventory(idle)

    def test_stale_events_and_corrupted_identity_or_reaction_ledgers_are_rejected(self):
        sim = self.simulator()
        event = next(e for e in sim.events() if e.kind == "oxygen_out")
        sim.apply(event)
        before_ids, before_vacancy = sim.oxygen_ids.copy(), sim.vacancy.copy()
        count = sim.event_count
        with self.assertRaisesRegex(ValueError, "stale"):
            sim.apply(event)
        self.assertEqual(sim.event_count, count)
        np.testing.assert_array_equal(sim.oxygen_ids, before_ids)
        np.testing.assert_array_equal(sim.vacancy, before_vacancy)
        sites = np.flatnonzero(~sim.vacancy)
        saved = int(sim.oxygen_ids[sites[1]])
        sim.oxygen_ids[sites[1]] = sim.oxygen_ids[sites[0]]
        with self.assertRaisesRegex(RuntimeError, "oxygen atom conservation"):
            sim.check()
        sim.oxygen_ids[sites[1]] = saved
        sim.electrons_to_electrode["bottom"] += 1
        with self.assertRaisesRegex(RuntimeError, "reaction electron ledger"):
            sim.check()
        sim.electrons_to_electrode["bottom"] -= 1
        self.assert_inventory(sim)

    def test_failed_candidate_field_solve_never_commits_any_event_kind(self):
        for kind in ("hop", "oxygen_out", "oxygen_in"):
            with self.subTest(kind=kind):
                sim = self.simulator(indices=[0, 1])
                event = next(e for e in sim.events() if e.kind == kind)
                before = transaction_state(sim)
                with patch.object(sim.network, "solve", side_effect=RuntimeError("injected candidate solve")) as solve:
                    with self.assertRaisesRegex(RuntimeError, "injected candidate solve"):
                        sim.apply(event)
                    self.assertEqual(solve.call_count, 1)
                self.assert_transaction_unchanged(sim, before)
                self.assert_inventory(sim)
                # The unsuccessful attempt did not age the event or consume an ID.
                sim.apply(event)
                self.assertEqual(sim.event_count, before["progress"][1] + 1)
                self.assert_inventory(sim)

    def test_failed_step_restores_time_and_rng_so_retry_matches_an_uninterrupted_run(self):
        sim = self.simulator(indices=[0, 1])
        uninterrupted = self.simulator(indices=[0, 1])
        self.assertEqual(sim.step(), uninterrupted.step())
        self.assertGreater(sim.time_s, 0.0)
        before = transaction_state(sim)
        with patch.object(sim.network, "solve", side_effect=RuntimeError("injected step solve")):
            with self.assertRaisesRegex(RuntimeError, "injected step solve"):
                sim.step()
        self.assert_transaction_unchanged(sim, before)
        self.assert_inventory(sim)
        self.assertEqual(sim.step(), uninterrupted.step())
        np.testing.assert_array_equal(sim.oxygen_ids, uninterrupted.oxygen_ids)
        np.testing.assert_array_equal(sim.vacancy_ids, uninterrupted.vacancy_ids)
        self.assert_inventory(sim)

    def test_failed_bias_solve_keeps_the_previous_voltage_epoch_and_fields(self):
        sim = self.simulator(indices=[0])
        sim.step()
        before = transaction_state(sim)
        with patch.object(sim.network, "solve", side_effect=RuntimeError("injected bias solve")):
            with self.assertRaisesRegex(RuntimeError, "injected bias solve"):
                sim.set_bias(-0.5)
        self.assert_transaction_unchanged(sim, before)
        self.assert_inventory(sim)

    def test_invalid_existing_ledger_prevents_candidate_solve_and_additional_mutation(self):
        sim = self.simulator()
        event = next(e for e in sim.events() if e.kind == "oxygen_out")
        sim.electrons_to_electrode["bottom"] += 1
        before = transaction_state(sim)
        with patch.object(sim.network, "solve", wraps=sim.network.solve) as solve:
            with self.assertRaisesRegex(RuntimeError, "reaction electron ledger"):
                sim.apply(event)
            solve.assert_not_called()
        self.assert_transaction_unchanged(sim, before)
        sim.electrons_to_electrode["bottom"] -= 1
        self.assert_inventory(sim)


if __name__ == "__main__":
    unittest.main()
