"""Boundary geometry and potential conventions independent of KMC rates."""
from copy import deepcopy
from dataclasses import FrozenInstanceError
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

from atomistic.model import CELL
from dftkmc.device import DeviceBoundary, DeviceConfig, ElectrodeConfig
from dftkmc.geometry import Lattice


class DeviceGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lattice = Lattice(repeats=(3, 3, 3), boundary='slab')

    def setUp(self):
        self.boundary = DeviceBoundary(self.lattice, DeviceConfig())

    def test_config_is_frozen_and_roundtrips_complete_strict_json(self):
        config = DeviceConfig(bottom=ElectrodeConfig(-0.25, 1.5), top=ElectrodeConfig(0.75, 2.0))
        self.assertEqual(DeviceConfig.from_dict(config.to_dict()), config)
        with self.assertRaises(FrozenInstanceError):
            config.top.potential_V = 3.0
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'device.json'
            config.write(path)
            self.assertEqual(DeviceConfig.load(path), config)
            self.assertEqual(json.loads(path.read_text()), config.to_dict())

    def test_config_rejects_missing_unknown_and_unsupported_physics(self):
        modifications = [
            (lambda d: d.pop('top'), 'missing fields'),
            (lambda d: d.__setitem__('unknown', 1), 'unknown fields'),
            (lambda d: d['top'].pop('potential_V'), 'missing fields'),
            (lambda d: d['bottom'].__setitem__('exchange_rate', 0), 'unknown fields'),
            (lambda d: d['top'].__setitem__('oxygen_boundary', 'reservoir'), 'blocking'),
            (lambda d: d.__setitem__('potential_model', 'poisson'), 'linear_dielectric'),
            (lambda d: d.__setitem__('schema_version', 2), 'schema_version'),
            (lambda d: d.__setitem__('schema_version', True), 'schema_version'),
            (lambda d: d.__setitem__('schema_version', 1.0), 'schema_version'),
        ]
        for modify, message in modifications:
            with self.subTest(message=message):
                raw = DeviceConfig().to_dict()
                modify(raw)
                with self.assertRaisesRegex(ValueError, message):
                    DeviceConfig.from_dict(raw)
        with self.assertRaisesRegex(ValueError, 'ElectrodeConfig'):
            DeviceConfig(bottom={})

    def test_config_rejects_nonfinite_boolean_and_invalid_contact_depths(self):
        for potential in (float('nan'), float('inf'), True, '1.0'):
            with self.subTest(potential=potential):
                with self.assertRaises(ValueError):
                    ElectrodeConfig(potential_V=potential)
        for depth in (0.0, -1.0, float('nan'), float('inf'), True, '2.0'):
            with self.subTest(depth=depth):
                with self.assertRaises(ValueError):
                    ElectrodeConfig(contact_depth_A=depth)
        with self.assertRaisesRegex(ValueError, 'potential difference'):
            DeviceConfig(bottom=ElectrodeConfig(-1e308), top=ElectrodeConfig(1e308))

    def test_duplicate_json_keys_and_nonfinite_json_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'device.json'
            for text in ('{"schema_version":1,"schema_version":2}', '{"top":NaN}'):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    DeviceConfig.load(path)

    def test_device_requires_finite_normal_direction_and_two_lateral_periodic_axes(self):
        with self.assertRaisesRegex(ValueError, 'requires a slab'):
            DeviceBoundary(Lattice(repeats=(2, 2, 2), boundary='bulk'), DeviceConfig())
        fake = SimpleNamespace(boundary='slab', pbc=(True, False, False))
        with self.assertRaisesRegex(ValueError, 'requires a slab'):
            DeviceBoundary(fake, DeviceConfig())
        invalid = SimpleNamespace(boundary='slab', pbc=(True, True, False),
                                  cell_A=np.array([[5., 0., 0.], [0., 5., 0.], [0., 0., -15.]]))
        with self.assertRaisesRegex(ValueError, 'C dot electrode normal'):
            DeviceBoundary(invalid, DeviceConfig())

    def test_monoclinic_thickness_uses_plane_separation_not_c_vector_length(self):
        a, b, c, beta = CELL
        expected = 3.0 * c * math.sin(math.radians(beta))
        np.testing.assert_allclose(self.boundary.normal, (0., 0., 1.), atol=1e-14)
        self.assertAlmostEqual(self.boundary.thickness_A, expected, 12)
        self.assertNotAlmostEqual(self.boundary.thickness_A, 3.0 * c, places=3)
        self.assertAlmostEqual(self.boundary.area_A2, (3 * a) * (3 * b), 10)
        self.assertAlmostEqual(self.boundary.field_V_nm[2], -10.0 / expected, 12)
        # The occupied atom span is also not the distance between cell-face electrodes.
        self.assertLess(float(np.ptp(self.boundary.depth_A)), self.boundary.thickness_A)

    def test_rigid_rotation_preserves_depth_and_rotates_the_physical_field(self):
        axis = np.array([1., 2., -1.])
        axis /= np.linalg.norm(axis)
        x, y, z = axis
        skew = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
        angle = 0.73
        rotation = math.cos(angle) * np.eye(3) + (1 - math.cos(angle)) * np.outer(axis, axis) + math.sin(angle) * skew
        reference = deepcopy(self.lattice.reference)
        reference['cell_A'] = (np.asarray(reference['cell_A']) @ rotation.T).tolist()
        reference['provenance'] += '; rigid rotation for a geometry test'
        rotated_lattice = Lattice(repeats=(3, 3, 3), boundary='slab', reference=reference)
        rotated = DeviceBoundary(rotated_lattice, DeviceConfig())
        np.testing.assert_allclose(rotated.normal, rotation @ np.array([0., 0., 1.]), atol=1e-14)
        np.testing.assert_allclose(rotated.field_V_nm, rotation @ self.boundary.field_V_nm, atol=1e-14)
        np.testing.assert_allclose(rotated.depth_A, self.boundary.depth_A, atol=1e-13)
        np.testing.assert_allclose(rotated.site_potential_V, self.boundary.site_potential_V, atol=1e-14)
        self.assertAlmostEqual(rotated.thickness_A, self.boundary.thickness_A, 12)
        np.testing.assert_array_equal(rotated.region, self.boundary.region)

    def test_face_potentials_gradient_and_lateral_periodic_images(self):
        boundary = DeviceBoundary(self.lattice, DeviceConfig(bottom=ElectrodeConfig(2.0), top=ElectrodeConfig(3.0)))
        a, b, c = self.lattice.cell_A
        np.testing.assert_allclose(boundary.potential_V(np.array([np.zeros(3), a, b, a+b])), 2.0, atol=1e-14)
        np.testing.assert_allclose(boundary.potential_V(np.array([c, c+a, c+b, c+a+b])), 3.0, atol=1e-14)
        point = 0.37 * c + 0.23 * a + 0.41 * b
        self.assertAlmostEqual(boundary.potential_V(point), 2.37, 12)
        self.assertAlmostEqual(boundary.potential_V(point + 5 * a - 6 * b), 2.37, 12)
        h = 1e-4
        numerical_field = np.array([-(boundary.potential_V(point+h*axis) - boundary.potential_V(point-h*axis))/(2*h)
                                    for axis in np.eye(3)])
        np.testing.assert_allclose(numerical_field, boundary.field_V_A, rtol=1e-10, atol=1e-11)
        self.assertAlmostEqual(float(boundary.field_V_A @ c), -1.0, 12)
        np.testing.assert_allclose(boundary.field_V_nm, 10.0 * boundary.field_V_A)

    def test_contact_tags_are_disjoint_oxide_regions_and_do_not_flatten_potential(self):
        d = self.boundary
        self.assertEqual(int(d.bottom_mask.sum()), 27)
        self.assertEqual(int(d.top_mask.sum()), 27)
        self.assertFalse(np.any(d.bottom_mask & d.top_mask))
        self.assertTrue(np.all(d.region[d.bottom_mask] == 'bottom_contact'))
        self.assertTrue(np.all(d.region[d.top_mask] == 'top_contact'))
        self.assertTrue(np.all(d.region[~d.bottom_mask & ~d.top_mask] == 'interior'))
        self.assertTrue(np.all(d.site_potential_V[d.bottom_mask] > d.config.bottom.potential_V))
        self.assertTrue(np.all(d.site_potential_V[d.top_mask] < d.config.top.potential_V))
        self.assertGreater(float(np.ptp(d.site_potential_V[d.bottom_mask])), 0)
        threshold = float(np.unique(d.depth_A)[1])
        exact = DeviceBoundary(self.lattice, DeviceConfig(bottom=ElectrodeConfig(contact_depth_A=threshold)))
        self.assertTrue(np.all(exact.bottom_mask[np.isclose(exact.depth_A, threshold, atol=1e-12)]))

    def test_empty_or_overlapping_contact_layers_are_not_silently_repaired(self):
        length = self.boundary.thickness_A
        for config, message in [
            (DeviceConfig(bottom=ElectrodeConfig(contact_depth_A=length / 2),
                          top=ElectrodeConfig(1.0, length / 2)), 'sum to less'),
            (DeviceConfig(bottom=ElectrodeConfig(contact_depth_A=0.001)), 'empty layers'),
            (DeviceConfig(top=ElectrodeConfig(1.0, 0.001)), 'empty layers'),
        ]:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    DeviceBoundary(self.lattice, config)

    def test_potential_rejects_outside_points_and_invalid_coordinate_arrays(self):
        for point in (-1e-3 * self.boundary.normal, self.lattice.cell_A[2] + 1e-3 * self.boundary.normal):
            with self.assertRaisesRegex(ValueError, 'outside the finite slab'):
                self.boundary.potential_V(point)
        for positions in (1., [0., 1.], [[0., 1., 2., 3.]], [0., True, 1.], [0., 0., float('nan')]):
            with self.subTest(positions=positions):
                with self.assertRaises(ValueError):
                    self.boundary.potential_V(positions)
        self.assertEqual(self.boundary.potential_V(np.empty((0, 3))).shape, (0,))

    def test_relative_energy_is_gauge_invariant_with_no_half_factor(self):
        shifted = DeviceBoundary(self.lattice, DeviceConfig(bottom=ElectrodeConfig(7.0), top=ElectrodeConfig(8.0)))
        np.testing.assert_allclose(shifted.relative_potential_V, self.boundary.relative_potential_V, atol=1e-14)
        np.testing.assert_allclose(shifted.site_potential_V, self.boundary.site_potential_V + 7.0, atol=1e-14)
        np.testing.assert_allclose(shifted.field_V_A, self.boundary.field_V_A, atol=1e-14)
        vacancies = np.zeros(self.lattice.n, dtype=bool)
        chosen = [int(np.argmin(self.boundary.depth_A)), int(np.argmax(self.boundary.depth_A))]
        vacancies[chosen] = True
        # Symmetric near-bottom/top depths have potentials summing to 1 V.
        self.assertAlmostEqual(self.boundary.electrostatic_energy_eV(vacancies, 2), 2.0, 12)
        self.assertAlmostEqual(shifted.electrostatic_energy_eV(vacancies, 2), 2.0, 12)
        self.assertEqual(self.boundary.electrostatic_energy_eV(vacancies, 0), 0.0)
        self.assertAlmostEqual(self.boundary.electrostatic_energy_eV(vacancies, -2), -2.0, 12)
        with self.assertRaises(ValueError):
            self.boundary.electrostatic_energy_eV(vacancies, 2.0)
        with self.assertRaises(ValueError):
            self.boundary.electrostatic_energy_eV(np.zeros(self.lattice.n), 2)

    def test_device_leaves_lattice_graphs_and_positions_unchanged(self):
        edges = list(self.lattice.migration_edges)
        pairs = self.lattice.interaction_pairs.copy()
        positions = self.lattice.oxygen_A.copy()
        adjacency = deepcopy(self.lattice.interaction_neighbors)
        device = DeviceBoundary(self.lattice, DeviceConfig())
        self.assertEqual(self.lattice.migration_edges, edges)
        np.testing.assert_array_equal(self.lattice.interaction_pairs, pairs)
        np.testing.assert_array_equal(self.lattice.oxygen_A, positions)
        self.assertEqual(self.lattice.interaction_neighbors, adjacency)
        for attribute in ('normal', 'depth_A', 'bottom_mask', 'top_mask', 'region',
                          'relative_potential_V', 'site_potential_V', 'field_V_A', 'field_V_nm'):
            self.assertFalse(getattr(device, attribute).flags.writeable)

    def test_summary_is_serializable_and_states_boundary_scope(self):
        summary = self.boundary.summary()
        json.dumps(summary, allow_nan=False)
        self.assertEqual(summary['counts'], {'oxide_sites': 216, 'bottom_contact': 27, 'interior': 162, 'top_contact': 27})
        self.assertEqual(summary['applied_voltage_V'], 1.0)
        self.assertEqual(summary['pbc'], [True, True, False])
        self.assertIs(summary['explicit_metal_atoms'], False)
        self.assertIs(summary['oxygen_exchange_enabled'], False)
        self.assertIsNone(summary['electronic_current_model'])
        self.assertIsNone(summary['poisson_solver'])
        self.assertIsNone(summary['interface_chemical_calibration'])
        self.assertIn('not equipotential metal nodes', summary['contacts']['role'])


if __name__ == '__main__':
    unittest.main()
