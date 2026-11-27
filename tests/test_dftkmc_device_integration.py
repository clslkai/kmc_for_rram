"""Independent physical and workflow acceptance for ideal electrode boundaries.

The oracle is the potential between two ideal parallel planes. These tests
verify boundary bookkeeping and ionic KMC, not electronic transport or forming.
"""
import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
import warnings

import numpy as np

from dftkmc.catalog import make_catalog
from dftkmc.cli import main
from dftkmc.device import DeviceConfig, ElectrodeConfig
from dftkmc.geometry import Lattice
from dftkmc.model import SimulationConfig, Simulator
from dftkmc.workflow import run_simulation, run_study


KB_EV_K = 8.617333262145e-5
REPEATS = (2, 2, 3)


def configuration(*, bottom=0.0, top=1.0, charge=2, seed=23):
    return SimulationConfig(
        repeats=REPEATS, boundary="slab", charge=charge, seed=seed,
        temperature_K=1000.0,
        device=DeviceConfig(
            bottom=ElectrodeConfig(potential_V=bottom, contact_depth_A=2.0),
            top=ElectrodeConfig(potential_V=top, contact_depth_A=2.0),
        ),
    )


def parallel_plane_oracle(lattice, bottom, top):
    cell = np.asarray(lattice.cell_A)
    normal = np.cross(cell[0], cell[1])
    normal /= np.linalg.norm(normal)
    height = float(cell[2] @ normal)
    depth = np.asarray(lattice.oxygen_A) @ normal
    potential = bottom + (top - bottom) * depth / height
    field_V_A = -(top - bottom) * normal / height
    return normal, height, depth, potential, field_V_A


def nontrivial_cycle(edges, depths):
    """Find a graph cycle with a real excursion in the electrode-normal direction."""
    neighbors = {}
    for edge in edges:
        neighbors.setdefault(edge.source, set()).add(edge.target)
        neighbors.setdefault(edge.target, set()).add(edge.source)
    root = next(iter(neighbors))
    parent = {root: None}
    queue = [root]
    for i in queue:
        for j in sorted(neighbors[i]):
            if j not in parent:
                parent[j] = i
                queue.append(j)
    for edge in edges:
        i, j = edge.source, edge.target
        if parent.get(i) == j or parent.get(j) == i:
            continue
        a, b = [i], [j]
        while parent[a[-1]] is not None:
            a.append(parent[a[-1]])
        while parent[b[-1]] is not None:
            b.append(parent[b[-1]])
        common = next(x for x in b if x in set(a))
        cycle = a[:a.index(common) + 1] + b[:b.index(common)][::-1]
        if len(cycle) >= 3 and np.ptp(np.asarray(depths)[cycle]) > 0.1:
            return cycle + cycle[:1]
    raise AssertionError("Validation geometry should contain a nontrivial normal-direction cycle")


class DeviceIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = make_catalog(Lattice(repeats=REPEATS))

    def test_device_snapshot_before_initialization_has_a_clear_precondition(self):
        sim = Simulator(self.catalog, configuration())
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            with self.assertRaisesRegex(ValueError, "initialize first"):
                sim.snapshot()

    def test_device_hop_matches_dirichlet_work_total_energy_and_reverse_rate(self):
        sim = Simulator(self.catalog, configuration())
        normal, height, depth, potential, field = parallel_plane_oracle(sim.lattice, 0.0, 1.0)
        np.testing.assert_allclose(sim.field_V_A, field, rtol=1e-13, atol=1e-15)
        np.testing.assert_allclose(sim.field_V_nm, 10.0 * field, rtol=1e-13, atol=1e-14)
        self.assertLess(float(sim.field_V_A @ normal), 0.0)
        edge = next(e for e in sim.edges if depth[e.target] - depth[e.source] > 0.1)
        sim.initialize(vacancy_indices=[edge.source])
        event = next(e for e in sim.events() if e.target == edge.target and e.channel_id == edge.channel_id)
        electrostatic_delta = 2.0 * (potential[event.target] - potential[event.source])
        self.assertGreater(electrostatic_delta, 0.0)  # Positive charge moving toward the positive top electrode.
        self.assertAlmostEqual(event.field_work_eV, electrostatic_delta, places=12)
        self.assertAlmostEqual(event.delta_eV, event.chemical_delta_eV + electrostatic_delta, places=12)
        before = sim.snapshot()
        sim.apply(event)
        after = sim.snapshot()
        self.assertAlmostEqual(after["total_energy_relative_bottom_eV"] - before["total_energy_relative_bottom_eV"],
                               event.delta_eV, places=11)
        self.assertAlmostEqual(after["electrostatic_energy_relative_bottom_eV"] - before["electrostatic_energy_relative_bottom_eV"],
                               electrostatic_delta, places=11)
        reverse = next(e for e in sim.events() if e.target == event.source and e.channel_id == event.channel_id)
        self.assertAlmostEqual(event.barrier_eV - reverse.barrier_eV, event.delta_eV, places=11)
        self.assertAlmostEqual(np.log(event.rate_Hz / reverse.rate_Hz),
                               -event.delta_eV / (KB_EV_K * 1000.0), places=10)

    def test_closed_slab_cycle_has_zero_electrostatic_work(self):
        sim = Simulator(self.catalog, configuration(bottom=0.3, top=0.9))
        _, _, depth, potential, _ = parallel_plane_oracle(sim.lattice, 0.3, 0.9)
        cycle = nontrivial_cycle(sim.edges, depth)
        sim.initialize(vacancy_indices=[cycle[0]])
        initial = sim.snapshot()
        work, energy_change = 0.0, 0.0
        for source, target in zip(cycle[:-1], cycle[1:]):
            event = next(e for e in sim.events() if e.source == source and e.target == target)
            self.assertAlmostEqual(event.field_work_eV, 2.0 * (potential[target] - potential[source]), places=12)
            work += event.field_work_eV
            energy_change += event.delta_eV
            sim.apply(event)
        self.assertAlmostEqual(work, 0.0, places=12)
        self.assertAlmostEqual(energy_change, 0.0, places=11)
        self.assertEqual(np.flatnonzero(sim.vacancy).tolist(), [cycle[0]])
        self.assertAlmostEqual(sim.snapshot()["total_energy_relative_bottom_eV"],
                               initial["total_energy_relative_bottom_eV"], places=11)

    def test_effective_graphs_keep_lateral_images_and_remove_cross_plane_wrapping(self):
        sim = Simulator(self.catalog, configuration())
        lattice = sim.lattice
        _, _, depth, _, _ = parallel_plane_oracle(lattice, 0.0, 1.0)
        inverse = np.linalg.inv(lattice.cell_A)
        lateral_crossings = 0
        for edge in sim.edges:
            raw = lattice.oxygen_A[edge.target] - lattice.oxygen_A[edge.source]
            image = (raw - np.asarray(edge.displacement_A)) @ inverse
            np.testing.assert_allclose(image, np.rint(image), atol=1e-10)
            self.assertAlmostEqual(image[2], 0.0, places=10)
            lateral_crossings += int(np.any(abs(image[:2]) > 0.5))
        self.assertGreater(lateral_crossings, 0)
        cutoff = float(self.catalog.charge_model(2)["pair_cutoff_A"])
        self.assertTrue(all(abs(depth[j] - depth[i]) <= cutoff + 1e-10
                            for i, j in lattice.interaction_pairs))
        self.assertTrue(any(np.linalg.norm(lattice.oxygen_A[j] - lattice.oxygen_A[i]) > cutoff
                            for i, j in lattice.interaction_pairs))
        bulk = Lattice(repeats=REPEATS, boundary="bulk", interaction_cutoff_A=cutoff)
        bulk_depth = bulk.oxygen_A @ parallel_plane_oracle(bulk, 0.0, 1.0)[0]
        forbidden_pairs = {frozenset((int(i), int(j))) for i, j in bulk.interaction_pairs
                           if abs(bulk_depth[j] - bulk_depth[i]) > cutoff + 1e-10}
        self.assertTrue(forbidden_pairs)
        actual_pairs = {frozenset((int(i), int(j))) for i, j in lattice.interaction_pairs}
        self.assertFalse(forbidden_pairs & actual_pairs)

    def test_common_voltage_offset_preserves_seeded_paths_and_relative_energy(self):
        sims = [Simulator(self.catalog, configuration(bottom=b, top=t)) for b, t in ((0.0, 1.0), (4.0, 5.0))]
        trajectories = []
        for sim in sims:
            sim.initialize(vacancy_indices=[5, 25, 60])
            trajectories.append(sim.run(1e-8, max_events=400, sample_every=7))
        self.assertGreater(sims[0].event_count, 0)
        self.assertEqual(sims[0].time_s, 1e-8)
        self.assertEqual(sims[0].event_count, sims[1].event_count)
        np.testing.assert_array_equal(sims[0].vacancy, sims[1].vacancy)
        np.testing.assert_array_equal(sims[0].labels, sims[1].labels)
        np.testing.assert_allclose(sims[0].displacements_A, sims[1].displacements_A, atol=1e-12)
        for a, b in zip(trajectories[0]["events"], trajectories[1]["events"]):
            self.assertEqual((a["source"], a["target"], a["channel_id"]),
                             (b["source"], b["target"], b["channel_id"]))
            self.assertAlmostEqual(a["rate_Hz"] / b["rate_Hz"], 1.0, places=12)
            self.assertAlmostEqual(a["time_s"], b["time_s"], places=18)
        for key in ("electrostatic_energy_relative_bottom_eV", "total_energy_relative_bottom_eV"):
            self.assertAlmostEqual(sims[0].snapshot()[key], sims[1].snapshot()[key], places=11)

    def test_neutral_vacancy_dynamics_are_bias_invariant(self):
        sims = [Simulator(self.catalog, configuration(top=v, charge=0)) for v in (-1.0, 0.0, 1.0)]
        for sim in sims:
            sim.initialize(vacancy_indices=[5])
            sim.run(1e-8, max_events=400, sample_every=20)
            self.assertEqual(sim.time_s, 1e-8)
            self.assertEqual(sim.snapshot()["electrostatic_energy_relative_bottom_eV"], 0.0)
        self.assertGreater(sims[0].event_count, 0)
        for sim in sims[1:]:
            self.assertEqual(sim.event_count, sims[0].event_count)
            np.testing.assert_array_equal(sim.vacancy, sims[0].vacancy)
            np.testing.assert_array_equal(sim.displacements_A, sims[0].displacements_A)

    def test_contact_counts_partition_conserved_vacancies(self):
        sim = Simulator(self.catalog, configuration())
        _, height, depth, _, _ = parallel_plane_oracle(sim.lattice, 0.0, 1.0)
        bottom = depth < 2.0
        top = height - depth < 2.0
        interior = ~(bottom | top)
        indices = [int(np.flatnonzero(mask)[0]) for mask in (bottom, top, interior)]
        sim.initialize(vacancy_indices=indices)
        for count in range(15):
            snapshot = sim.snapshot()
            expected = [int(np.count_nonzero(sim.vacancy & mask)) for mask in (bottom, top, interior)]
            observed = [snapshot[k] for k in ("bottom_contact_vacancies", "top_contact_vacancies", "interior_vacancies")]
            self.assertEqual(observed, expected)
            self.assertEqual(sum(observed), 3)
            np.testing.assert_array_equal(np.sort(sim.labels[sim.vacancy]), [0, 1, 2])
            events = sim.events()
            sim.apply(events[(count * 7) % len(events)])
        for unsupported in ("current_A", "conductance_S", "resistance_ohm", "spanning"):
            self.assertNotIn(unsupported, sim.snapshot())

    def test_device_export_potentials_events_and_contact_masks_share_the_same_geometry(self):
        config = configuration(bottom=-0.2, top=0.8)
        lattice = Lattice(repeats=REPEATS, boundary="slab")
        _, height, depth, potential, field = parallel_plane_oracle(lattice, -0.2, 0.8)
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "device_run"
            result = run_simulation(self.catalog, config, vacancy_indices=[5, 25, 60],
                                    duration_s=1e-8, max_events=400, sample_every=7,
                                    output=output, allow_illustrative=True)
            self.assertTrue(result["target_reached"])
            self.assertGreater(result["events"], 0)
            summary = json.loads((output / "summary.json").read_text())
            sidecar = json.loads((output / "device_boundary.json").read_text())
            self.assertEqual(summary["device_boundary"], sidecar)
            np.testing.assert_allclose(summary["effective_field_V_nm"], field * 10.0)
            with (output / "device_sites.csv").open() as handle:
                sites = {int(row["site_id"]): row for row in csv.DictReader(handle)}
            self.assertEqual(set(sites), set(range(lattice.n)))
            for i, row in sites.items():
                self.assertAlmostEqual(float(row["depth_A"]), depth[i], places=10)
                self.assertAlmostEqual(float(row["potential_V"]), potential[i], places=12)
            with np.load(output / "states.npz") as states:
                self.assertEqual(states["potential_V"].shape[-1], lattice.n)
                np.testing.assert_allclose(states["potential_V"], np.broadcast_to(potential, states["potential_V"].shape), atol=1e-12)
                np.testing.assert_allclose(states["device_depth_A"], depth, atol=1e-11)
                np.testing.assert_array_equal(states["bottom_contact"], depth < 2.0)
                np.testing.assert_array_equal(states["top_contact"], height - depth < 2.0)
                self.assertTrue(np.all(states["vacancy"].sum(axis=1) == 3))
            with (output / "events.csv").open() as handle:
                events = list(csv.DictReader(handle))
            self.assertEqual(len(events), result["events"])
            for event in events:
                i, j = int(event["source"]), int(event["target"])
                self.assertAlmostEqual(float(event["source_potential_V"]), potential[i], places=12)
                self.assertAlmostEqual(float(event["target_potential_V"]), potential[j], places=12)
                self.assertEqual(event["source_region"], sites[i]["region"])
                self.assertEqual(event["target_region"], sites[j]["region"])
                self.assertAlmostEqual(float(event["field_work_eV"]), 2.0 * (potential[j] - potential[i]), places=11)
                self.assertAlmostEqual(float(event["applied_voltage_V"]), 1.0, places=12)
                np.testing.assert_allclose(json.loads(event["field_V_nm"]), field * 10.0)

    def test_cli_selects_slab_and_rejects_explicit_conflicting_drives(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            root = Path(tmp)
            catalogue = root / "catalog.json"
            device = root / "device.json"
            self.catalog.write(catalogue)
            self.assertEqual(main(["device-template", "--output", str(device),
                                   "--bottom-voltage", "0.25", "--top-voltage", "1.25",
                                   "--bottom-contact-depth", "2", "--top-contact-depth", "2"]), 0)
            common = ["run", "--catalog", str(catalogue), "--allow-illustrative", "--device", str(device),
                      "--repeats", "2", "2", "3", "--duration", "1e-12", "--max-events", "5"]
            output = root / "implicit_slab"
            self.assertEqual(main(common + ["--output", str(output)]), 0)
            summary = json.loads((output / "summary.json").read_text())
            self.assertEqual(summary["config"]["boundary"], "slab")
            self.assertIsNotNone(summary["device_boundary"])
            self.assertTrue(np.any(np.asarray(summary["effective_field_V_nm"]) != 0.0))
            for index, extra in enumerate((["--boundary", "bulk"], ["--field", "0", "0", "0"], ["--field", "0", "0", "1"])):
                rejected = root / f"rejected_{index}"
                self.assertEqual(main(common + extra + ["--output", str(rejected)]), 2)
                self.assertFalse(rejected.exists())

    def test_python_conflicts_and_slab_study_do_not_produce_bulk_diffusion(self):
        device = configuration().device
        for kwargs in ({"boundary": "bulk"}, {"boundary": "slab", "field_V_nm": (0.0, 0.0, 0.1)}):
            with self.assertRaises(ValueError):
                Simulator(self.catalog, SimulationConfig(device=device, **kwargs))
        with tempfile.TemporaryDirectory() as tmp:
            result = run_study(self.catalog, configuration(top=0.0), n_vacancies=2,
                               duration_s=1e-12, replicas=2, charges=(2,), ablations=("noninteracting",),
                               max_events=5, bootstrap_samples=100, output=Path(tmp) / "study",
                               allow_illustrative=True)
            self.assertTrue(result["groups"][0]["statistics"]["all_replicas_complete"])
            self.assertFalse(result["groups"][0]["statistics"]["diffusion"]["available"])


if __name__ == "__main__":
    unittest.main()
