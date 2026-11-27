"""Analytic circuits and geometry/conservation checks for electronic transport.

These checks validate the declared illustrative resistor closure.  They do not
validate a real material conductance, filament threshold or forming voltage.
"""
from dataclasses import FrozenInstanceError, replace
import json
import math
from types import SimpleNamespace
import unittest

import numpy as np

from dftkmc.device import DeviceBoundary, DeviceConfig, ElectrodeConfig
from dftkmc.geometry import Lattice
from dftkmc.transport import ElectronicNetwork, TransportConfig


def geometry(positions, cell=None, contact_depth=1.1):
    positions = np.asarray(positions, dtype=float)
    lattice = SimpleNamespace(
        n=len(positions), oxygen_A=positions,
        cell_A=np.diag([20., 20., 6.]) if cell is None else np.asarray(cell, dtype=float),
        boundary="slab", pbc=(True, True, False),
    )
    device = DeviceBoundary(lattice, DeviceConfig(
        bottom=ElectrodeConfig(0.0, contact_depth), top=ElectrodeConfig(1.0, contact_depth)))
    return lattice, device


def chain(config=None):
    lattice, device = geometry([[0., 0., 1.], [0., 0., 3.], [0., 0., 5.]])
    if config is None:
        config = TransportConfig(electronic_cutoff_A=2.1,
                                 background_conductance_S=2., vacancy_conductance_S=20.)
    return ElectronicNetwork(lattice, device, config)


class ElectronicTransportTests(unittest.TestCase):
    def test_series_resistors_include_finite_contacts_and_match_closed_form(self):
        network = chain()
        solution = network.solve(np.zeros(3, dtype=bool))
        # Each internal edge has r=2 A and g=2/2=1 S.  Both contacts have
        # r_eff=1 A plane separation + 1 A offset, giving another two 1 S
        # resistors.  Four series 1-ohm resistors give G=1/4 S.
        self.assertAlmostEqual(solution.equivalent_conductance_S, 0.25, 14)
        self.assertAlmostEqual(solution.terminal_current_A, 0.25, 14)
        self.assertAlmostEqual(solution.bottom_terminal_current_A, -0.25, 14)
        np.testing.assert_allclose(solution.potential_V, [0.25, 0.5, 0.75], atol=1e-14)
        np.testing.assert_allclose(solution.edge_current_A, [-0.25, -0.25], atol=1e-14)
        np.testing.assert_allclose(solution.edge_field_V_nm, [-1.25, -1.25], atol=1e-14)
        np.testing.assert_allclose(solution.edge_power_W, [0.0625, 0.0625], atol=1e-14)
        np.testing.assert_allclose(solution.bottom_contact_power_W, [0.0625], atol=1e-14)
        np.testing.assert_allclose(solution.top_contact_power_W, [0.0625], atol=1e-14)
        self.assertAlmostEqual(solution.total_power_W, 0.25, 14)
        # Oxide contact sites are not hard-clamped metal nodes.
        self.assertGreater(solution.potential_V[0], solution.bottom_potential_V)
        self.assertLess(solution.potential_V[-1], solution.top_potential_V)

    def test_parallel_branches_add_conductance_and_resolve_a_vacancy_path(self):
        positions = [[x, 0., z] for x in (0., 6.) for z in (1., 3., 5.)]
        lattice, device = geometry(positions)
        config = chain().config
        network = ElectronicNetwork(lattice, device, config)
        # Two spatially separated three-site chains share the same electrodes.
        vacancy = np.array([False, False, False, True, True, True])
        solution = network.solve(vacancy)
        self.assertAlmostEqual(solution.equivalent_conductance_S, 0.25 + 2.5, 13)
        np.testing.assert_allclose(solution.potential_V, [.25, .5, .75, .25, .5, .75], atol=1e-14)
        np.testing.assert_allclose(solution.top_contact_current_A, [.25, 2.5], atol=1e-14)
        self.assertTrue(solution.terminal_network_connected)
        self.assertTrue(solution.vacancy_graph_spanning)
        self.assertEqual(solution.conducting_path_site_ids, (3, 4, 5))
        self.assertAlmostEqual(solution.total_power_W, 2.75, 13)

    def test_pristine_leakage_and_partial_path_bottleneck_are_distinct_from_spanning(self):
        network = chain()
        pristine = network.solve(np.zeros(3, dtype=bool))
        partial = network.solve(np.array([True, True, False]))
        complete = network.solve(np.ones(3, dtype=bool))
        self.assertTrue(pristine.terminal_network_connected)
        self.assertGreater(pristine.terminal_current_A, 0.)
        self.assertFalse(pristine.vacancy_graph_spanning)
        self.assertFalse(partial.vacancy_graph_spanning)
        self.assertEqual(partial.conducting_path_site_ids, ())
        # Partial chain resistances are 0.1, 0.1, 1, 1 ohms in that order.
        self.assertAlmostEqual(partial.equivalent_conductance_S, 1. / 2.2, 14)
        np.testing.assert_allclose(partial.potential_V, np.array([.1, .2, 1.2]) / 2.2, atol=1e-14)
        self.assertAlmostEqual(complete.equivalent_conductance_S, 2.5, 13)
        self.assertTrue(complete.vacancy_graph_spanning)
        self.assertEqual(complete.conducting_path_site_ids, (0, 1, 2))
        self.assertIsNone(complete.summary()['electronic_on_threshold'])

    def test_equal_vacancy_and_background_conductance_removes_electrical_state_dependence(self):
        network = chain(replace(chain().config, vacancy_conductance_S=2.))
        states = [np.zeros(3, dtype=bool), np.array([True, False, True]), np.ones(3, dtype=bool)]
        solutions = [network.solve(state) for state in states]
        for solution in solutions[1:]:
            np.testing.assert_array_equal(solution.potential_V, solutions[0].potential_V)
            np.testing.assert_array_equal(solution.edge_conductance_S, solutions[0].edge_conductance_S)
            self.assertEqual(solution.equivalent_conductance_S, solutions[0].equivalent_conductance_S)
            self.assertEqual(solution.terminal_current_A, solutions[0].terminal_current_A)
        # Geometry still changes although this configured electrical contrast is 1.
        self.assertFalse(solutions[0].vacancy_graph_spanning)
        self.assertTrue(solutions[-1].vacancy_graph_spanning)

    def test_voltage_gauge_polarity_zero_bias_and_maximum_principle(self):
        network = chain()
        vacancy = np.array([True, True, False])
        positive = network.solve(vacancy, bottom_V=0., top_V=1.)
        shifted = network.solve(vacancy, bottom_V=7.5, top_V=8.5)
        negative = network.solve(vacancy, bottom_V=0., top_V=-1.)
        zero = network.solve(vacancy, bottom_V=-2., top_V=-2.)
        np.testing.assert_allclose(shifted.potential_V, positive.potential_V + 7.5, atol=1e-14)
        np.testing.assert_array_equal(shifted.relative_potential_V, positive.relative_potential_V)
        np.testing.assert_array_equal(shifted.edge_current_A, positive.edge_current_A)
        np.testing.assert_array_equal(shifted.edge_power_W, positive.edge_power_W)
        np.testing.assert_allclose(negative.potential_V, -positive.potential_V, atol=1e-14)
        np.testing.assert_allclose(negative.edge_field_V_nm, -positive.edge_field_V_nm, atol=1e-14)
        self.assertEqual(negative.terminal_current_A, -positive.terminal_current_A)
        self.assertEqual(negative.total_power_W, positive.total_power_W)
        self.assertEqual(zero.equivalent_conductance_S, positive.equivalent_conductance_S)
        self.assertEqual(zero.terminal_current_A, 0.)
        self.assertEqual(zero.total_power_W, 0.)
        np.testing.assert_array_equal(zero.potential_V, [-2., -2., -2.])
        np.testing.assert_array_equal(zero.edge_field_V_nm, [0., 0.])
        for solution in (positive, shifted, negative, zero):
            self.assertGreaterEqual(float(solution.potential_V.min()), min(solution.bottom_potential_V, solution.top_potential_V) - 1e-14)
            self.assertLessEqual(float(solution.potential_V.max()), max(solution.bottom_potential_V, solution.top_potential_V) + 1e-14)

    def test_ideal_compliance_returns_actual_voltage_and_frozen_unclamped_read(self):
        network = chain(replace(chain().config, current_compliance_A=.1, read_voltage_V=2.))
        vacancy = np.zeros(3, dtype=bool)
        solution = network.solve(vacancy, bottom_V=3., top_V=8.)
        self.assertEqual(solution.requested_voltage_V, 5.)
        self.assertTrue(solution.compliance_active)
        self.assertAlmostEqual(solution.applied_voltage_V, .4, 14)
        self.assertEqual(solution.bottom_potential_V, 3.)
        self.assertAlmostEqual(solution.top_potential_V, 3.4, 14)
        np.testing.assert_allclose(solution.potential_V, [3.1, 3.2, 3.3], atol=1e-14)
        self.assertAlmostEqual(solution.terminal_current_A, .1, 14)
        self.assertAlmostEqual(solution.total_power_W, .04, 14)
        self.assertAlmostEqual(solution.read_current_A, .5, 14)
        self.assertGreater(solution.read_current_A, solution.current_compliance_A)
        reverse = network.solve(vacancy, bottom_V=3., top_V=-2.)
        self.assertAlmostEqual(reverse.applied_voltage_V, -.4, 14)
        self.assertAlmostEqual(reverse.top_potential_V, 2.6, 14)
        self.assertAlmostEqual(reverse.terminal_current_A, -.1, 14)
        low_bias = network.solve(vacancy, bottom_V=3., top_V=3.1)
        self.assertFalse(low_bias.compliance_active)

    def test_electronic_graph_does_not_read_migration_or_interaction_graphs(self):
        lattice = Lattice(repeats=(2, 2, 2), boundary='slab')
        device = DeviceBoundary(lattice, DeviceConfig())
        first = ElectronicNetwork(lattice, device)
        positions_before = lattice.oxygen_A.copy()
        self.assertTrue(np.any(first.distances_A > lattice.interaction_cutoff_A))

        class ForbiddenGraph:
            def __iter__(self):
                raise AssertionError('ionic graph was read by electronic transport')

            def __array__(self, *args, **kwargs):
                raise AssertionError('interaction graph was read by electronic transport')

        lattice.migration_edges = ForbiddenGraph()
        lattice.interaction_pairs = ForbiddenGraph()
        lattice.interaction_neighbors = ForbiddenGraph()
        second = ElectronicNetwork(lattice, device)
        np.testing.assert_array_equal(first.pairs, second.pairs)
        np.testing.assert_array_equal(first.displacements_A, second.displacements_A)
        np.testing.assert_array_equal(lattice.oxygen_A, positions_before)
        vacancy = np.zeros(lattice.n, dtype=bool)
        np.testing.assert_array_equal(first.solve(vacancy).potential_V, second.solve(vacancy).potential_V)

    def test_skew_c_lateral_images_rigid_rotation_and_no_normal_wrap(self):
        # A large shear in finite C makes full-cell fractional x/y rounding
        # inappropriate without a stricter full-cell bound.  The physical
        # neighbor is obtained by translating A only, never C.
        cell = np.array([[20., 0., 0.], [7., 20., 0.], [78., 0., 8.]])
        positions = np.array([[0., 0., z / 8.] for z in (1., 3., 5., 7.)]) @ cell
        lattice, device = geometry(positions, cell)
        config = replace(chain().config, electronic_cutoff_A=2.2)
        network = ElectronicNetwork(lattice, device, config)
        np.testing.assert_array_equal(network.pairs, [[0, 1], [1, 2], [2, 3]])
        np.testing.assert_allclose(network.displacements_A, np.tile([-.5, 0., 2.], (3, 1)), atol=1e-13)
        # Site 0 and site 3 must not connect through a fictitious wrapped C.
        self.assertNotIn((0, 3), {tuple(pair) for pair in network.pairs})
        expected_g = 2. / (4. + 3. * math.sqrt(4.25))
        solution = network.solve(np.zeros(4, dtype=bool))
        self.assertAlmostEqual(solution.equivalent_conductance_S, expected_g, 13)
        # An arbitrary lateral image choice and rigid spatial rotation leave
        # a circuit unchanged, while Cartesian displacement vectors rotate.
        shifted = positions.copy()
        shifted[1] += 3 * cell[0] - 2 * cell[1]
        angle = .61
        axis = np.array([1., -2., 3.]) / math.sqrt(14.)
        x, y, z = axis
        skew = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
        rotation = math.cos(angle) * np.eye(3) + (1 - math.cos(angle)) * np.outer(axis, axis) + math.sin(angle) * skew
        rotated_lattice, rotated_device = geometry(shifted @ rotation.T, cell @ rotation.T)
        rotated = ElectronicNetwork(rotated_lattice, rotated_device, config)
        np.testing.assert_array_equal(rotated.pairs, network.pairs)
        np.testing.assert_allclose(rotated.displacements_A, network.displacements_A @ rotation.T, atol=1e-13)
        rotated_solution = rotated.solve(np.zeros(4, dtype=bool))
        np.testing.assert_allclose(rotated_solution.potential_V, solution.potential_V, atol=1e-13)
        self.assertAlmostEqual(rotated_solution.equivalent_conductance_S, expected_g, 13)

    def test_disconnected_one_electrode_components_are_valid_zero_current(self):
        lattice, device = geometry([[0., 0., 1.], [0., 0., 5.]])
        config = replace(chain().config, current_compliance_A=.01)
        network = ElectronicNetwork(lattice, device, config)
        solution = network.solve(np.ones(2, dtype=bool), bottom_V=-1., top_V=2.)
        self.assertEqual(network.pairs.shape, (0, 2))
        self.assertFalse(solution.terminal_network_connected)
        self.assertFalse(solution.vacancy_graph_spanning)
        self.assertEqual(solution.conducting_path_site_ids, ())
        np.testing.assert_array_equal(solution.potential_V, [-1., 2.])
        self.assertEqual(solution.equivalent_conductance_S, 0.)
        self.assertEqual(solution.terminal_current_A, 0.)
        self.assertEqual(solution.read_current_A, 0.)
        self.assertEqual(solution.total_power_W, 0.)
        self.assertEqual(solution.applied_voltage_V, 3.)
        self.assertFalse(solution.compliance_active)
        np.testing.assert_array_equal(solution.node_kcl_residual_A, [0., 0.])

    def test_floating_components_and_ambiguous_periodic_images_fail_explicitly(self):
        lattice, device = geometry([[0., 0., 1.], [0., 0., 3.], [0., 0., 5.]])
        with self.assertRaisesRegex(ValueError, 'floating.*potentials undefined'):
            ElectronicNetwork(lattice, device, replace(chain().config, electronic_cutoff_A=1.1))
        with self.assertRaisesRegex(ValueError, 'unique electronic periodic images'):
            ElectronicNetwork(lattice, device, replace(chain().config, electronic_cutoff_A=10.))
        with self.assertRaisesRegex(ValueError, 'boolean occupancy'):
            chain().solve([0, 0, 0])

    def test_conservation_power_and_rayleigh_monotonicity_on_actual_hafnia_geometry(self):
        lattice = Lattice(repeats=(2, 2, 2), boundary='slab')
        device = DeviceBoundary(lattice, DeviceConfig())
        network = ElectronicNetwork(lattice, device)
        vacancy = np.zeros(lattice.n, dtype=bool)
        pristine = network.solve(vacancy, top_V=1.7)
        vacancy[np.arange(0, lattice.n, 3)] = True
        partial = network.solve(vacancy, top_V=1.7)
        complete = network.solve(np.ones(lattice.n, dtype=bool), top_V=1.7)
        self.assertLess(pristine.equivalent_conductance_S, partial.equivalent_conductance_S)
        self.assertLess(partial.equivalent_conductance_S, complete.equivalent_conductance_S)
        for solution in (pristine, partial, complete):
            # Independent invariants of a passive two-terminal resistor circuit.
            self.assertLess(solution.relative_kcl_residual, 1e-8)
            self.assertLess(solution.relative_terminal_current_imbalance, 1e-8)
            self.assertLess(solution.relative_power_balance_residual, 1e-8)
            self.assertTrue(np.all(solution.edge_power_W >= 0))
            self.assertTrue(np.all(solution.bottom_contact_power_W >= 0))
            self.assertTrue(np.all(solution.top_contact_power_W >= 0))
            self.assertLessEqual(abs(solution.bottom_terminal_current_A + solution.terminal_current_A),
                                 abs(solution.terminal_current_A) * 1e-8)
            measured_power = (solution.edge_power_W.sum() + solution.bottom_contact_power_W.sum()
                              + solution.top_contact_power_W.sum())
            self.assertAlmostEqual(measured_power / (solution.applied_voltage_V * solution.terminal_current_A), 1., 8)
        path = complete.conducting_path_site_ids
        self.assertIn(path[0], network.bottom_sites)
        self.assertIn(path[-1], network.top_sites)
        edges = {tuple(sorted(pair)) for pair in network.pairs}
        self.assertTrue(all(tuple(sorted(pair)) in edges for pair in zip(path[:-1], path[1:])))

    def test_complete_configuration_and_serializable_explicit_scope(self):
        config = TransportConfig()
        self.assertEqual(TransportConfig.from_dict(config.to_dict()), config)
        with self.assertRaises(FrozenInstanceError):
            config.contact_offset_A = 0.
        for field in ('background_conductance_S', 'vacancy_conductance_S', 'contact_offset_A',
                      'electronic_cutoff_A', 'reference_distance_A', 'read_voltage_V'):
            for value in (0., -1., float('nan'), float('inf'), True, '1.0'):
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(config, **{field: value})
        with self.assertRaises(ValueError):
            replace(config, vacancy_conductance_S=config.background_conductance_S / 2)
        with self.assertRaises(ValueError):
            replace(config, current_compliance_A=0.)
        with self.assertRaises(ValueError):
            replace(config, provenance='  ')
        with self.assertRaises(ValueError):
            replace(config, calibration_status='dft')
        raw = config.to_dict()
        raw.pop('contact_offset_A')
        with self.assertRaisesRegex(ValueError, 'missing fields'):
            TransportConfig.from_dict(raw)
        raw = config.to_dict()
        raw['temperature_feedback'] = True
        with self.assertRaisesRegex(ValueError, 'unknown fields'):
            TransportConfig.from_dict(raw)
        network = chain()
        solution = network.solve(np.ones(3, dtype=bool))
        json.dumps(network.summary(), allow_nan=False)
        json.dumps(solution.summary(), allow_nan=False)
        json.dumps(solution.to_dict(), allow_nan=False)
        self.assertNotIn('potential_V', solution.to_dict(include_arrays=False))
        self.assertIn('potential_V', solution.to_dict())
        self.assertIs(solution.summary()['temperature_feedback'], False)
        self.assertIs(solution.summary()['quantum_transport'], False)
        self.assertEqual(solution.calibration_status, 'illustrative')
        for array in (network.pairs, network.displacements_A, network.distances_A,
                      network.bottom_sites, network.top_sites, solution.potential_V,
                      solution.edge_current_A, solution.node_kcl_residual_A):
            self.assertFalse(array.flags.writeable)


if __name__ == '__main__':
    unittest.main()
