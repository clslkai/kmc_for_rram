"""Independent invariants for the explicitly synthetic raw-data import path."""
from copy import deepcopy
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from dftkmc.catalog import Catalog
from dftkmc.mockdata import (
    SCHEMA_VERSION, import_mock_bundle, load_mock_bundle, make_mock_bundle,
    source_digest, write_mock_bundle,
)
from dftkmc.model import KB, SimulationConfig, Simulator


class MockDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = make_mock_bundle()

    def setUp(self):
        self.bundle = deepcopy(self.source)

    def mixed_path(self):
        return next(path for path in self.bundle['neb_paths']
                    if path['from_kind'] == 'O3' and path['to_kind'] == 'O4')

    def test_default_bundle_covers_full_geometry_and_remains_illustrative(self):
        catalog = import_mock_bundle(self.bundle)
        self.assertEqual(self.bundle['schema_version'], SCHEMA_VERSION)
        self.assertEqual(len(self.bundle['neb_paths']), 34)
        self.assertEqual(catalog.charges, (0, 2))
        self.assertEqual(catalog.data['metadata']['calibration_status'], 'illustrative')
        self.assertIs(catalog.data['metadata']['synthetic_data'], True)
        self.assertEqual(catalog.data['metadata']['import_consistency']['checked_channel_charge_records'], 68)
        self.assertTrue(catalog.audit()['ready'])
        self.assertFalse(catalog.audit()['declared_dft_calibrated'])
        with self.assertRaisesRegex(ValueError, 'allow_illustrative'):
            catalog.validate()
        self.assertEqual(self.bundle, make_mock_bundle())

    def test_hand_chosen_barriers_produce_a_common_saddle_and_correct_rates(self):
        # Independent fixture: eps_A=0, eps_B=0.1, B_A->B=0.9, B_B->A=0.8.
        # The absolute saddle must be 0.9, and the uphill penalty must not be
        # added again by either import or rate evaluation.
        path = self.mixed_path()
        path['data']['0']['barrier_forward_eV'] = 0.9
        path['data']['0']['barrier_reverse_eV'] = 0.8
        catalog = import_mock_bundle(self.bundle)
        self.assertAlmostEqual(catalog.channel(path['id'], 0)['saddle_eV'], 0.9)
        simulation = Simulator(catalog, SimulationConfig(repeats=(2, 2, 2), charge=0, temperature_K=600.0))
        edge = next(edge for edge in simulation.edges if edge.channel_id == path['id'])
        simulation.initialize(vacancy_indices=[edge.source])
        forward = next(event for event in simulation.events()
                       if event.channel_id == edge.channel_id and event.source == edge.source
                       and event.target == edge.target)
        self.assertAlmostEqual(forward.barrier_eV, 0.9)
        self.assertAlmostEqual(forward.rate_Hz / (1e12 * math.exp(-0.9 / (KB * 600.0))), 1.0)
        simulation.apply(forward)
        reverse = next(event for event in simulation.events()
                       if event.channel_id == edge.channel_id and event.source == edge.target
                       and event.target == edge.source)
        self.assertAlmostEqual(reverse.barrier_eV, 0.8)
        self.assertAlmostEqual(math.log(forward.rate_Hz / reverse.rate_Hz), -0.1 / (KB * 600.0))

    def test_inconsistent_reverse_barrier_is_rejected_independently(self):
        self.mixed_path()['data']['0']['barrier_reverse_eV'] += 0.025
        with self.assertRaisesRegex(ValueError, 'inconsistent forward/reverse barriers'):
            import_mock_bundle(self.bundle)

    def test_reverse_barrier_cannot_be_inferred_when_missing(self):
        del self.mixed_path()['data']['0']['barrier_reverse_eV']
        with self.assertRaisesRegex(ValueError, 'missing required fields.*barrier_reverse_eV'):
            import_mock_bundle(self.bundle)

    def test_prefactor_asymmetry_and_unsupported_energy_fields_are_rejected(self):
        self.mixed_path()['data']['0']['attempt_reverse_Hz'] = 2e12
        with self.assertRaisesRegex(ValueError, 'unequal directional prefactors'):
            import_mock_bundle(self.bundle)
        self.bundle = deepcopy(self.source)
        self.mixed_path()['data']['0']['saddle_eV'] = 0.9
        with self.assertRaisesRegex(ValueError, 'unsupported fields.*saddle_eV'):
            import_mock_bundle(self.bundle)

    def test_canonical_site_identity_survives_reordering(self):
        reference = import_mock_bundle(self.bundle)
        for records in self.bundle['site_energies'].values():
            records.reverse()
        restored = import_mock_bundle(self.bundle)
        for q in restored.charges:
            np.testing.assert_array_equal(restored.site_energies(q), reference.site_energies(q))
        self.assertNotEqual(source_digest(self.bundle), source_digest(self.source))
        self.bundle['site_energies']['0'][0]['basis_id'] = self.bundle['site_energies']['0'][1]['basis_id']
        with self.assertRaisesRegex(ValueError, 'duplicate site-energy basis ID'):
            import_mock_bundle(self.bundle)

    def test_missing_site_charge_and_channel_data_fail_closed(self):
        mutations = [
            (lambda b: b['site_energies']['0'].pop(), 'eight labelled records'),
            (lambda b: b['site_energies'].pop('2'), 'declared canonical charge keys'),
            (lambda b: b['charge_models'].pop('2'), 'declared canonical charge keys'),
            (lambda b: b['neb_paths'][0]['data'].pop('2'), 'declared canonical charge keys'),
            (lambda b: b['neb_paths'].pop(), 'uncovered geometry channels'),
        ]
        for mutate, message in mutations:
            with self.subTest(message=message):
                bundle = deepcopy(self.source)
                mutate(bundle)
                with self.assertRaisesRegex(ValueError, message):
                    import_mock_bundle(bundle)

    def test_charge_declaration_cannot_alias_or_silently_remove_a_branch(self):
        for charges in ([False, 2], [0, 2.0], [0, 0], []):
            with self.subTest(charges=charges):
                bundle = deepcopy(self.source)
                bundle['metadata']['charges'] = charges
                with self.assertRaisesRegex(ValueError, 'distinct integral charges'):
                    import_mock_bundle(bundle)
        self.bundle['site_energies']['02'] = self.bundle['site_energies'].pop('2')
        with self.assertRaisesRegex(ValueError, 'canonical charge keys'):
            import_mock_bundle(self.bundle)

    def test_reference_fingerprint_displacement_and_labels_are_checked(self):
        mutations = [
            (lambda b: b['metadata'].__setitem__('structure_id', '0' * 64), 'structure_id'),
            (lambda b: b['neb_paths'][0]['displacement_A'].__setitem__(0, 99.0), 'displacement_A'),
            (lambda b: b['site_energies']['0'][0].__setitem__('kind', 'O4'), 'kind disagrees'),
            (lambda b: b['neb_paths'][0].__setitem__('id', 'invented'), 'canonical geometry'),
            (lambda b: b['neb_paths'][0].__setitem__('from_basis', True), 'integer basis index'),
            (lambda b: b['neb_paths'].append(deepcopy(b['neb_paths'][0])), 'duplicate channel ID'),
        ]
        for mutate, message in mutations:
            with self.subTest(message=message):
                bundle = deepcopy(self.source)
                mutate(bundle)
                with self.assertRaisesRegex(ValueError, message):
                    import_mock_bundle(bundle)

    def test_nonfinite_null_negative_and_boolean_values_are_not_replaced(self):
        mutations = [
            (lambda b: b['site_energies']['0'][0].__setitem__('energy_eV', None), 'null'),
            (lambda b: b['site_energies']['0'][0].__setitem__('energy_eV', True), 'finite number'),
            (lambda b: b['neb_paths'][0]['data']['0'].__setitem__('barrier_forward_eV', float('nan')), 'finite'),
            (lambda b: b['neb_paths'][0]['data']['0'].__setitem__('barrier_reverse_eV', None), 'null'),
            (lambda b: b['neb_paths'][0]['data']['0'].__setitem__('barrier_forward_eV', -0.1), 'nonnegative'),
            (lambda b: b['neb_paths'][0]['data']['0'].__setitem__('attempt_forward_Hz', 0.0), 'positive'),
            (lambda b: b['neb_paths'][0]['data']['0'].__setitem__('ts_fraction', 1.1), r'\[0, 1\]'),
            (lambda b: b['charge_models']['0']['pair_energies_eV'].__setitem__('O3-O3', None), 'null'),
            (lambda b: b['charge_models']['0'].__setitem__('environment_saddle_eV_per_neighbor', float('inf')), 'finite'),
        ]
        for mutate, message in mutations:
            with self.subTest(message=message):
                bundle = deepcopy(self.source)
                mutate(bundle)
                with self.assertRaisesRegex(ValueError, message):
                    import_mock_bundle(bundle)

    def test_units_and_classical_energy_assumptions_are_not_guessed(self):
        for field, value in [('energy_unit', 'J'), ('length_unit', 'nm'), ('frequency_unit', 'THz')]:
            with self.subTest(field=field):
                bundle = deepcopy(self.source)
                bundle['metadata'][field] = value
                with self.assertRaisesRegex(ValueError, 'no implicit unit conversion'):
                    import_mock_bundle(bundle)
        self.bundle['metadata']['vibrational_entropy_included'] = True
        with self.assertRaisesRegex(ValueError, 'vibrational entropy'):
            import_mock_bundle(self.bundle)

    def test_synthetic_data_cannot_be_promoted_or_have_its_marker_removed(self):
        mutations = [
            lambda b: b['metadata'].__setitem__('calibration_status', 'dft'),
            lambda b: b['metadata'].__setitem__('synthetic_data', False),
            lambda b: b['site_energies']['0'][0].__setitem__('status', 'dft'),
            lambda b: b['charge_models']['2'].__setitem__('status', 'dft'),
            lambda b: b['neb_paths'][0]['data']['0'].__setitem__('status', 'dft'),
            lambda b: b['neb_paths'][0]['data']['0'].__setitem__('synthetic_data', False),
        ]
        for index, mutate in enumerate(mutations):
            with self.subTest(index=index):
                bundle = deepcopy(self.source)
                mutate(bundle)
                with self.assertRaisesRegex(ValueError, 'DFT status is forbidden'):
                    import_mock_bundle(bundle)
        self.bundle['metadata']['status'] = 'dft'
        with self.assertRaisesRegex(ValueError, 'unsupported fields.*status'):
            import_mock_bundle(self.bundle)

    def test_all_input_provenance_is_explicit_and_preserved(self):
        catalog = import_mock_bundle(self.bundle)
        self.assertEqual(catalog.data['metadata']['mock_source']['provenance'], self.bundle['metadata']['provenance'])
        self.assertEqual(catalog.data['metadata']['site_energy_provenance'], self.bundle['site_energies'])
        self.assertEqual(catalog.structure['provenance'], self.bundle['structure']['provenance'])
        for template in self.bundle['neb_paths']:
            for q in catalog.charges:
                self.assertEqual(catalog.channel(template['id'], q)['provenance'], template['data'][str(q)]['provenance'])
        for locate in [lambda b: b['metadata'], lambda b: b['site_energies']['0'][0],
                       lambda b: b['charge_models']['0'], lambda b: b['neb_paths'][0]['data']['0']]:
            bundle = deepcopy(self.source)
            locate(bundle)['provenance'] = 'completed calculation'
            with self.assertRaisesRegex(ValueError, 'SYNTHETIC provenance label'):
                import_mock_bundle(bundle)

    def test_digest_records_complete_source_and_detects_consistent_data_changes(self):
        original = deepcopy(self.bundle)
        digest = source_digest(self.bundle)
        catalog = import_mock_bundle(self.bundle, expected_sha256=digest)
        self.assertEqual(self.bundle, original)
        self.assertEqual(catalog.data['metadata']['mock_source']['sha256'], digest)
        altered = self.mixed_path()['data']['0']
        altered['barrier_forward_eV'] += 0.03
        altered['barrier_reverse_eV'] += 0.03
        # The changed pair remains thermodynamically consistent but it is not the
        # retained source.  A checksum check must catch this independent risk.
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            import_mock_bundle(self.bundle, expected_sha256=digest)
        self.assertNotEqual(import_mock_bundle(self.bundle).data['metadata']['mock_source']['sha256'], digest)
        catalog.data['metadata']['mock_source']['source_metadata']['provenance'] = 'modified catalog copy'
        self.assertEqual(original['metadata']['provenance'], self.source['metadata']['provenance'])

    def test_energy_zero_shift_preserves_imported_activation_barriers(self):
        original = import_mock_bundle(self.bundle)
        # Shifting every state energy in one branch shifts its reconstructed
        # saddles but cannot change a measured activation barrier in that branch.
        for record in self.bundle['site_energies']['0']:
            record['energy_eV'] += 3.75
        shifted = import_mock_bundle(self.bundle)
        for template in self.bundle['neb_paths']:
            a = template['from_basis']
            old = original.channel(template['id'], 0)['saddle_eV'] - original.site_energies(0)[a]
            new = shifted.channel(template['id'], 0)['saddle_eV'] - shifted.site_energies(0)[a]
            self.assertAlmostEqual(new, old)
        np.testing.assert_array_equal(original.site_energies(2), shifted.site_energies(2))

    def test_json_roundtrip_and_catalog_opt_in(self):
        with tempfile.TemporaryDirectory() as directory:
            raw, ready = Path(directory) / 'source.json', Path(directory) / 'catalog.json'
            write_mock_bundle(self.bundle, raw)
            restored = load_mock_bundle(raw, expected_sha256=source_digest(self.bundle))
            self.assertEqual(restored, self.bundle)
            self.assertEqual(source_digest(restored), source_digest(self.bundle))
            catalog = import_mock_bundle(restored)
            catalog.write(ready)
            with self.assertRaisesRegex(ValueError, 'allow_illustrative'):
                Catalog.load(ready)
            self.assertEqual(Catalog.load(ready, allow_illustrative=True).data, catalog.data)

    def test_invalid_write_preserves_existing_file_and_strict_parser_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'source.json'
            path.write_text('existing contents')
            self.bundle['metadata']['calibration_status'] = 'dft'
            with self.assertRaises(ValueError):
                write_mock_bundle(self.bundle, path)
            self.assertEqual(path.read_text(), 'existing contents')
            for contents in ('{"schema_version":1,"schema_version":2}', '{"number":NaN}'):
                path.write_text(contents)
                with self.assertRaises(ValueError):
                    load_mock_bundle(path)

    def test_module_cli_generates_and_imports_explicitly_synthetic_files(self):
        with tempfile.TemporaryDirectory() as directory:
            raw, ready = Path(directory) / 'source.json', Path(directory) / 'catalog.json'
            generate = subprocess.run([sys.executable, '-m', 'dftkmc.mockdata', 'generate', '--output', str(raw)],
                                      capture_output=True, text=True, check=True)
            summary = json.loads(generate.stdout)
            self.assertEqual(summary['status'], 'illustrative')
            self.assertTrue(summary['synthetic_data'])
            self.assertEqual(summary['paths'], 34)
            subprocess.run([sys.executable, '-m', 'dftkmc.mockdata', 'import', '--input', str(raw),
                            '--output', str(ready), '--expected-sha256', summary['source_sha256']],
                           capture_output=True, text=True, check=True)
            catalog = Catalog.load(ready, allow_illustrative=True)
            self.assertEqual(catalog.data['metadata']['mock_source']['sha256'], summary['source_sha256'])


if __name__ == '__main__':
    unittest.main()
