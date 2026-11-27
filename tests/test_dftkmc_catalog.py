"""Scientific-input validation gates, not material-parameter validation."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from dftkmc.catalog import Catalog, make_catalog, PAIR_KEYS
from dftkmc.geometry import Lattice, channel_id


class CatalogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lattice = Lattice(repeats=(2, 2, 2))
        cls.original = make_catalog(cls.lattice).data

    def setUp(self):
        self.catalog = Catalog(self.original)

    def assert_invalid(self, catalog, message):
        report = catalog.audit()
        self.assertFalse(report['ready'])
        self.assertIn(message, '\n'.join(report['errors']))
        with self.assertRaises(ValueError):
            catalog.validate(allow_illustrative=True)

    def test_illustrative_policy_and_classical_conventions_are_explicit(self):
        report = self.catalog.audit()
        self.assertTrue(report['ready'], report['errors'])
        self.assertTrue(report['requires_allow_illustrative'])
        self.assertFalse(report['declared_dft_calibrated'])
        self.assertIn('Synthetic', ' '.join(report['warnings']))
        self.assertIn('No vibrational entropy', ' '.join(report['warnings']))
        self.assertIn('not environment-resolved NEB', ' '.join(report['warnings']))
        with self.assertRaisesRegex(ValueError, 'allow_illustrative'):
            self.catalog.validate()
        self.catalog.validate(allow_illustrative=True)
        self.assertEqual(self.catalog.charges, (0, 2))

    def test_pending_worklist_has_nulls_and_cannot_run(self):
        pending = make_catalog(self.lattice, status='pending')
        for q in pending.charges:
            self.assertTrue(all(value is None for value in pending.data['site_energies_eV'][str(q)]))
            model = pending.charge_model(q)
            self.assertEqual(model['environment_mode'], 'dilute_only')
            self.assertTrue(all(value is None for value in model['pair_energies_eV'].values()))
            for template in pending.data['channels']:
                data = pending.channel(template['id'], q)
                self.assertIsNone(data['saddle_eV'])
                self.assertIsNone(data['attempt_Hz'])
                self.assertEqual(data['status'], 'pending')
                self.assertIn('placeholder', data['provenance'])
        report = pending.audit()
        self.assertGreater(len(report['missing']), 20)
        self.assert_invalid(pending, 'cannot run')
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            pending.site_energies(0)
        with self.assertRaisesRegex(ValueError, 'manufacture DFT'):
            make_catalog(self.lattice, status='dft')

    def test_json_roundtrip_and_pending_inspection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'catalog.json'
            self.catalog.write(path)
            with self.assertRaisesRegex(ValueError, 'allow_illustrative'):
                Catalog.load(path)
            restored = Catalog.load(path, allow_illustrative=True)
            self.assertEqual(restored.data, self.catalog.data)
            pending = make_catalog(self.lattice, status='pending')
            pending.write(path)
            self.assertIn('null', path.read_text())
            with self.assertRaises(ValueError):
                Catalog.load(path, allow_illustrative=True)
            self.assertEqual(Catalog.load(path, validate=False).data, pending.data)

    def test_duplicate_json_keys_and_nonfinite_json_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'catalog.json'
            for contents, message in [('{"schema_version":1,"schema_version":2}', 'duplicate'),
                                      ('{"value":NaN}', 'non-finite'),
                                      ('{"value":Infinity}', 'non-finite')]:
                with self.subTest(contents=contents):
                    path.write_text(contents)
                    with self.assertRaisesRegex(ValueError, message):
                        Catalog.load(path, validate=False)
            path.write_text('original file')
            self.catalog.data['metadata']['invalid_number'] = float('nan')
            with self.assertRaises(ValueError):
                self.catalog.write(path)
            self.assertEqual(path.read_text(), 'original file')

    def test_missing_channels_or_charge_branches_never_use_fallback(self):
        first = self.catalog.data['channels'][0]
        with self.assertRaisesRegex(ValueError, 'no fallback'):
            self.catalog.channel('unknown-channel', 0)
        with self.assertRaisesRegex(ValueError, 'no fallback'):
            self.catalog.channel(first['id'], 7)
        with self.assertRaisesRegex(ValueError, 'cover charge'):
            self.catalog.site_energies(7)
        with self.assertRaisesRegex(ValueError, 'cover charge'):
            self.catalog.charge_model(7)
        del first['data']['2']
        self.assert_invalid(self.catalog, 'charge coverage')
        self.catalog = Catalog(self.original)
        self.catalog.data['channels'].pop()
        self.assert_invalid(self.catalog, 'uncovered geometry channels')

    def test_explicit_channel_truncation_is_distinct_from_missing_data(self):
        removed = self.catalog.data['channels'].pop()
        self.catalog.data['metadata']['channel_exclusions'] = {
            removed['id']: 'Declared reduced-path model while awaiting this NEB calculation.'
        }
        self.catalog.validate(allow_illustrative=True)
        report = self.catalog.audit()
        self.assertEqual(report['n_excluded_channels'], 1)
        self.assertIn('truncated channel graph', ' '.join(report['warnings']))
        self.assertIn('cannot establish full diffusion', ' '.join(report['warnings']))
        self.catalog.data['metadata']['channel_exclusions'][removed['id']] = ' '
        self.assert_invalid(self.catalog, 'truncation reason')
        self.catalog.data['metadata']['channel_exclusions'] = {'unknown': 'intentional exclusion'}
        self.assert_invalid(self.catalog, 'unknown excluded geometry channels')
        self.catalog = Catalog(self.original)
        self.catalog.data['metadata']['channel_exclusions'] = {
            self.catalog.data['channels'][0]['id']: 'This overlaps an actual path.'
        }
        self.assert_invalid(self.catalog, 'overlap actual channels')

    def test_structure_fingerprint_and_displacement_are_checked(self):
        cases = [
            ('hash', lambda c: c['metadata'].__setitem__('structure_id', '0' * 64), 'structure_id'),
            ('basis', lambda c: c['structure']['oxygen_fractional'][0].__setitem__(0, 0.111), 'displacement_A'),
            ('displacement', lambda c: c['channels'][0]['displacement_A'].__setitem__(0, 0.111), 'displacement_A'),
            ('length', lambda c: c['channels'][0].__setitem__('length_A', 0.111), 'length_A'),
            ('kind', lambda c: c['channels'][0].__setitem__('from_kind', 'wrong-kind'), 'O3/O4'),
            ('channel ID', lambda c: c['channels'][0].__setitem__('id', 'unknown'), 'canonical geometry'),
        ]
        for name, mutate, message in cases:
            with self.subTest(name=name):
                catalog = Catalog(self.original)
                mutate(catalog.data)
                self.assert_invalid(catalog, message)

    def test_reverse_duplicates_and_noncanonical_reverse_are_rejected(self):
        reverse = deepcopy(self.catalog.data['channels'][0])
        reverse['from_basis'], reverse['to_basis'] = reverse['to_basis'], reverse['from_basis']
        reverse['from_kind'], reverse['to_kind'] = reverse['to_kind'], reverse['from_kind']
        reverse['cell_shift'] = [-v for v in reverse['cell_shift']]
        reverse['displacement_A'] = [-v for v in reverse['displacement_A']]
        reverse['id'] = channel_id(reverse['from_basis'], reverse['to_basis'], reverse['cell_shift'])
        self.catalog.data['channels'].append(reverse)
        self.assert_invalid(self.catalog, 'or its reverse')
        self.catalog = Catalog(self.original)
        self.catalog.data['channels'][0] = reverse
        self.assert_invalid(self.catalog, 'not a canonical geometry channel')

    def test_duplicate_ids_are_rejected_by_both_audit_and_lookup(self):
        self.catalog.data['channels'].append(deepcopy(self.catalog.data['channels'][0]))
        self.assert_invalid(self.catalog, 'duplicate channel ID')
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            self.catalog.template(self.catalog.data['channels'][0]['id'])

    def test_endpoint_and_translation_indices_must_be_integral(self):
        for key, value in [('from_basis', -1), ('to_basis', 8), ('from_basis', 0.5),
                           ('to_basis', True), ('cell_shift', [0., 0, 0]),
                           ('cell_shift', [0, 0]), ('cell_shift', [False, 0, 0])]:
            with self.subTest(key=key, value=value):
                catalog = Catalog(self.original)
                catalog.data['channels'][0][key] = value
                self.assert_invalid(catalog, 'integers')

    def test_charge_keys_are_canonical_and_charge_requests_are_integral(self):
        for key in ('2.0', '02', '+2', '-0', 'neutral', 2):
            with self.subTest(key=key):
                catalog = Catalog(self.original)
                catalog.data['site_energies_eV'][key] = catalog.data['site_energies_eV'].pop('2')
                self.assert_invalid(catalog, 'charge key')
        for value in (2.0, 0.5, True, '2'):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    self.catalog.site_energies(value)
        for charges in ((), (0, 0), (2.5,), (True,)):
            with self.subTest(charges=charges):
                with self.assertRaises(ValueError):
                    make_catalog(self.lattice, charges=charges)
        self.assertEqual(make_catalog(self.lattice, charges=(2, -1, 0)).charges, (-1, 0, 2))

    def test_nonfinite_or_boolean_scientific_numbers_are_rejected(self):
        mutations = [
            lambda c, value: c['site_energies_eV']['0'].__setitem__(0, value),
            lambda c, value: c['charge_models']['0']['pair_energies_eV'].__setitem__('O3-O3', value),
            lambda c, value: c['charge_models']['0'].__setitem__('pair_cutoff_A', value),
            lambda c, value: c['charge_models']['0'].__setitem__('environment_saddle_eV_per_neighbor', value),
            lambda c, value: c['channels'][0]['data']['0'].__setitem__('saddle_eV', value),
            lambda c, value: c['channels'][0]['data']['0'].__setitem__('attempt_Hz', value),
            lambda c, value: c['channels'][0]['data']['0'].__setitem__('ts_fraction', value),
            lambda c, value: c['structure']['cell_A'][0].__setitem__(0, value),
        ]
        for index, mutate in enumerate(mutations):
            for value in (float('nan'), float('inf'), True):
                with self.subTest(index=index, value=value):
                    catalog = Catalog(self.original)
                    mutate(catalog.data, value)
                    report = catalog.audit()
                    self.assertFalse(report['ready'])
                    self.assertIn('finite', '\n'.join(report['errors']))

    def test_shared_saddle_is_above_both_endpoints_and_gauge_invariant(self):
        first = self.catalog.data['channels'][0]
        sites = self.catalog.site_energies(0)
        a, b = first['from_basis'], first['to_basis']
        first['data']['0']['saddle_eV'] = max(sites[a], sites[b]) - 0.01
        self.assert_invalid(self.catalog, 'below an endpoint')
        shifted = Catalog(self.original)
        for q in shifted.charges:
            key = str(q)
            shifted.data['site_energies_eV'][key] = [value + 7.25 for value in shifted.data['site_energies_eV'][key]]
            for template in shifted.data['channels']:
                template['data'][key]['saddle_eV'] += 7.25
        shifted.validate(allow_illustrative=True)
        for original, moved in zip(self.original['channels'], shifted.data['channels']):
            a, b = original['from_basis'], original['to_basis']
            saddle = moved['data']['0']['saddle_eV']
            new_sites = shifted.site_energies(0)
            forward, reverse = saddle - new_sites[a], saddle - new_sites[b]
            self.assertAlmostEqual(forward - reverse, new_sites[b] - new_sites[a])
            self.assertAlmostEqual(forward, original['data']['0']['saddle_eV'] - self.original['site_energies_eV']['0'][a])

    def test_attempt_frequency_and_saddle_fraction_bounds(self):
        for key, value, message in [('attempt_Hz', 0, 'positive'), ('attempt_Hz', -1, 'positive'),
                                    ('ts_fraction', -0.1, '[0, 1]'), ('ts_fraction', 1.1, '[0, 1]')]:
            with self.subTest(key=key, value=value):
                catalog = Catalog(self.original)
                catalog.data['channels'][0]['data']['0'][key] = value
                self.assert_invalid(catalog, message)

    def test_asymmetric_prefactors_and_vibrational_entropy_cannot_be_ignored(self):
        self.catalog.data['channels'][0]['data']['0']['attempt_reverse_Hz'] = 2e13
        self.assert_invalid(self.catalog, 'one shared attempt_Hz')
        self.catalog = Catalog(self.original)
        self.catalog.data['metadata']['vibrational_entropy_included'] = True
        self.assert_invalid(self.catalog, 'vibrational entropy')
        self.catalog = Catalog(self.original)
        self.catalog.data['metadata']['rate_convention'] = 'unequal_htst_prefactors'
        self.assert_invalid(self.catalog, 'equal forward/reverse')

    def test_units_and_provenance_are_required(self):
        mutations = [
            (lambda c: c['metadata'].__setitem__('energy_unit', 'J'), 'energy_unit'),
            (lambda c: c['metadata'].__setitem__('length_unit', 'nm'), 'length_unit'),
            (lambda c: c['metadata'].__setitem__('provenance', ''), 'metadata.provenance'),
            (lambda c: c['structure'].__setitem__('provenance', 'unspecified reference'), 'structure.provenance'),
            (lambda c: c['charge_models']['0'].__setitem__('provenance', None), 'charge_models[0].provenance'),
            (lambda c: c['channels'][0]['data']['0'].__setitem__('provenance', ' '), 'provenance'),
        ]
        for mutate, message in mutations:
            with self.subTest(message=message):
                catalog = Catalog(self.original)
                mutate(catalog.data)
                self.assert_invalid(catalog, message)

    def test_status_cannot_silently_launder_synthetic_data(self):
        self.catalog.data['metadata']['calibration_status'] = 'dft'
        self.assert_invalid(self.catalog, 'synthetic_data cannot be promoted')
        self.assertIn('status disagrees', '\n'.join(self.catalog.audit()['errors']))
        self.catalog = Catalog(self.original)
        self.catalog.data['channels'][0]['data']['0']['status'] = 'dft'
        self.assert_invalid(self.catalog, 'status disagrees')

    def test_declared_dilute_mode_has_no_environment_saddle_correction(self):
        model = self.catalog.charge_model(0)
        model['environment_mode'] = 'dilute_only'
        self.assert_invalid(self.catalog, 'dilute_only requires zero')
        model['environment_saddle_eV_per_neighbor'] = 0.0
        self.catalog.validate(allow_illustrative=True)
        del model['pair_energies_eV'][PAIR_KEYS[0]]
        self.assert_invalid(self.catalog, 'exactly')

    def test_input_is_copied_and_site_arrays_are_independent(self):
        source = deepcopy(self.original)
        catalog = Catalog(source)
        source['metadata']['structure_id'] = 'changed after construction'
        self.assertEqual(catalog.structure_id, self.lattice.structure_id)
        sites = catalog.site_energies(0)
        sites[0] = 99
        self.assertEqual(catalog.site_energies(0)[0], self.original['site_energies_eV']['0'][0])

    def test_audit_handles_malformed_shapes_without_guessing(self):
        for key, value in [('metadata', []), ('structure', []), ('site_energies_eV', []),
                           ('charge_models', []), ('channels', []), ('channels', [None])]:
            with self.subTest(key=key, value=value):
                catalog = Catalog(self.original)
                catalog.data[key] = value
                self.assertFalse(catalog.audit()['ready'])

    def test_catalogue_identity_does_not_depend_on_simulation_repeats(self):
        larger = Lattice(repeats=(3, 2, 2))
        other = make_catalog(larger)
        self.assertEqual(other.structure_id, self.catalog.structure_id)
        self.assertEqual(other.data['channels'], self.catalog.data['channels'])

    def test_valid_longer_cutoff_does_not_get_rejected_by_audit_supercell_size(self):
        # This pair cutoff needs more than two repeats for unique image geometry.
        larger = Lattice(repeats=(3, 3, 3), interaction_cutoff_A=5.1)
        catalog = make_catalog(larger)
        catalog.validate(allow_illustrative=True)


if __name__ == '__main__':
    unittest.main()
