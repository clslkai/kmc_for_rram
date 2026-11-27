import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from dftkmc.analysis import finite_time_tracer_estimate, summarize_replicas
from dftkmc.catalog import make_catalog
from dftkmc.cli import main
from dftkmc.geometry import Lattice
from dftkmc.model import SimulationConfig
from dftkmc.workflow import export_neb, run_simulation, run_study


class ReplicaStatisticsTests(unittest.TestCase):
    def test_centering_and_correlated_tracers_do_not_shrink_uncertainty(self):
        values = np.array([[-2., 1., 0.], [1., -1., 2.], [3., 2., -1.], [-1., -2., 0.]])
        single = [x[None, :] for x in values]
        duplicated = [np.tile(x, (10, 1)) for x in values]
        shifted = [x + np.array([100., -50., 37.]) for x in single]
        args = dict(duration_s=2., bootstrap_samples=200, bootstrap_seed=99)
        a = finite_time_tracer_estimate(single, **args)
        b = finite_time_tracer_estimate(duplicated, **args)
        c = finite_time_tracer_estimate(shifted, **args)
        expected = np.cov(values.T, ddof=1) / 4. * 1e-20
        np.testing.assert_allclose(a['tensor_m2_s'], expected, atol=1e-34)
        np.testing.assert_allclose(a['tensor_m2_s'], b['tensor_m2_s'], atol=1e-34)
        np.testing.assert_allclose(a['tensor_ci_high_m2_s'], b['tensor_ci_high_m2_s'], atol=1e-34)
        np.testing.assert_allclose(a['tensor_m2_s'], c['tensor_m2_s'], rtol=1e-10, atol=1e-31)
        self.assertEqual(a['n_independent_replicas'], 4)

    def test_identical_zero_motion_does_not_claim_zero_uncertainty(self):
        report = finite_time_tracer_estimate([np.zeros((1,3)) for _ in range(8)], duration_s=1., bootstrap_samples=100)
        self.assertEqual(report["trace_over_3_m2_s"],0.)
        self.assertIsNone(report["trace_over_3_ci_m2_s"])
        self.assertFalse(report["uncertainty_resolved"])

    def test_censoring_disables_inference_not_just_warning(self):
        records = [dict(seed=1, status='completed_duration', time_s=1., final={'msd_A2': 4.}, displacements_A=np.array([[2.,0,0]])),
                   dict(seed=2, status='event_limit_before_end', time_s=.5, final={'msd_A2': 3.}, displacements_A=np.array([[1.,1,1.]]))]
        report = summarize_replicas(records, duration_s=1., bootstrap_samples=100)
        self.assertEqual(report['n_complete_replicas'], 1)
        self.assertFalse(report['diffusion']['available'])
        self.assertIsNone(report['scalar_estimates']['msd_A2']['ci'])
        self.assertEqual(report['scalar_estimates']['msd_A2']['mean'], 4.)

    def test_only_equal_time_bulk_zero_field_has_diffusion(self):
        records = [dict(seed=i, status='completed_duration', time_s=1e-12,
                        final={'msd_A2': 1.}, displacements_A=np.array([[float(i),0,0]])) for i in range(2)]
        for conditions in ({'boundary':'slab'}, {'field_V_nm':(0.,0.,.001)}):
            report = summarize_replicas(records, duration_s=1e-12, bootstrap_samples=100, **conditions)
            self.assertFalse(report['diffusion']['available'])
        # An overshoot is not an observation at the requested physical time.
        records[-1]['time_s'] = 2e-12
        report = summarize_replicas(records, duration_s=1e-12, bootstrap_samples=100)
        self.assertFalse(report['all_replicas_complete'])


class WorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lattice = Lattice()
        cls.catalog = make_catalog(cls.lattice)

    def test_run_trace_replays_identities_and_unwrapped_displacements(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'trajectory'
            result = run_simulation(self.catalog, SimulationConfig(temperature_K=1000., seed=4),
                                    vacancy_indices=[100, 2, 179], duration_s=1e-8, max_events=300,
                                    sample_every=3, output=output, allow_illustrative=True)
            self.assertTrue(result['target_reached'])
            self.assertGreater(result['events'], 0)
            summary = json.loads((output/'summary.json').read_text())
            self.assertEqual(summary['initial_sites_by_vacancy_id'], [100,2,179])
            states = np.load(output/'states.npz')
            self.assertEqual(states['time_s'][-1], 1e-8)
            self.assertEqual(states['vacancy'].shape[0], len(states['time_s']))
            labels = states['initial_labels'].copy()
            displacement = np.zeros((3,3))
            with (output/'events.csv').open() as handle:
                rows = list(csv.DictReader(handle))
            for row in rows:
                source,target,identity = int(row['source']),int(row['target']),int(row['vacancy_id'])
                self.assertEqual(labels[source], identity)
                self.assertEqual(labels[target], -1)
                labels[source],labels[target] = -1,identity
                displacement[identity] += json.loads(row['displacement_A'])
                self.assertGreater(float(row['dt_s']), 0)
                self.assertEqual(row['catalog_sha256'], summary['catalog_sha256'])
                self.assertEqual(row['channel_status'], 'illustrative')
            np.testing.assert_array_equal(labels, states['final_labels'])
            np.testing.assert_allclose(displacement, states['displacements_A'][-1])

    def test_run_event_cap_is_explicitly_censored(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = run_simulation(self.catalog, SimulationConfig(), duration_s=1., max_events=1,
                                    output=Path(temporary)/'capped', allow_illustrative=True)
            self.assertFalse(result['target_reached'])
            self.assertTrue(result['censored'])
            self.assertEqual(result['status'], 'event_limit_before_end')
            self.assertEqual(result['events'], 1)

    def test_neb_cross_boundary_uses_short_displacement_and_same_atom_ids(self):
        channel = next(t for t in self.lattice.channel_templates() if any(x < 0 for x in t['cell_shift']))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)/'neb'
            result = export_neb(self.lattice, output=output, charges=(2,), channel_ids=[channel['id']])
            task = output/channel['id']/'q2'
            mapping = json.loads((task/'mapping.json').read_text())
            first_lines = (task/'POSCAR.initial').read_text().splitlines()
            last_lines = (task/'POSCAR.final').read_text().splitlines()
            self.assertEqual(first_lines[5:8], last_lines[5:8])
            first = np.loadtxt(first_lines[8:])
            last = np.loadtxt(last_lines[8:])
            moving = mapping['moving_oxygen_poscar_index_1based']-1
            changed = np.flatnonzero(np.linalg.norm(last-first,axis=1)>1e-9)
            np.testing.assert_array_equal(changed,[moving])
            np.testing.assert_allclose(last[moving]-first[moving], -np.asarray(channel['displacement_A']), atol=1e-12)
            source = np.array(mapping['vacancy_source_wrapped_A'])
            image_offset = (last[moving]-source) @ np.linalg.inv(self.lattice.cell_A)
            np.testing.assert_allclose(image_offset, np.round(image_offset), atol=1e-12)
            self.assertLess(np.linalg.norm(last[moving]-first[moving]), 3.2+1e-9)
            self.assertTrue(mapping['relaxation_required'])
            self.assertEqual(mapping['charge'],2)
            self.assertFalse(list(output.rglob('POTCAR')))
            pending = json.loads((output/'catalog.pending.json').read_text())
            self.assertIsNone(pending['channels'][0]['data']['2']['saddle_eV'])
            self.assertTrue(result['no_jobs_submitted'])

    def test_study_independent_seeds_and_no_cross_charge_energy_contrast(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)/'study'
            result = run_study(self.catalog, SimulationConfig(seed=981), duration_s=1e-12, replicas=2,
                               charges=(0,2), ablations=('noninteracting','pair-only'),
                               bootstrap_samples=100, output=output, allow_illustrative=True)
            raw = result['replica_outcomes']
            self.assertEqual(len(raw),8)
            self.assertEqual(len({r['seed'] for r in raw}),8)
            self.assertTrue(all(r['time_s']==1e-12 for r in raw))
            for contrast in result['contrasts']:
                if 'energy_comparison' in contrast:
                    self.assertNotIn('energy_eV',contrast['statistics']['differences'])
                else:
                    self.assertIn('energy_eV',contrast['statistics']['differences'])
            self.assertTrue((output/'replicas.csv').exists())
            self.assertEqual(len(result['groups']),4)

    def test_cli_pending_catalogue_blocks_run_and_audit_reports_it(self):
        with tempfile.TemporaryDirectory() as temporary, contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            path = Path(temporary)/'pending.json'
            self.assertEqual(main(['template','--output',str(path)]),0)
            self.assertEqual(main(['audit',str(path)]),1)
            self.assertEqual(main(['run','--catalog',str(path),'--duration','1e-10','--output',str(Path(temporary)/'run')]),2)
            self.assertFalse((Path(temporary)/'run').exists())

    def test_declared_exclusions_propagate_and_neb_keeps_them_only_as_provenance(self):
        truncated = make_catalog(self.lattice)
        excluded = truncated.data['channels'].pop(0)['id']
        exclusions = {excluded: 'Deliberate reduced-graph sensitivity test; this channel is omitted.'}
        truncated.data['metadata']['channel_exclusions'] = exclusions
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # Two effectively zero-event replicas exercise serialized output without a costly trajectory.
            study = run_study(truncated, SimulationConfig(), duration_s=1e-15,
                              replicas=2, charges=(2,), ablations=('noninteracting',),
                              max_events=1, bootstrap_samples=100,
                              output=root/'study', allow_illustrative=True)
            first_run = next((root/'study').glob('q*/replica_*/summary.json'))
            summary = json.loads(first_run.read_text())
            for report in (summary, study, json.loads((root/'study'/'study.json').read_text())):
                self.assertEqual(report['migration_cutoff_A'], 3.2)
                self.assertEqual(report['channel_exclusions'], exclusions)
                self.assertTrue(report['conditional_on_truncated_graph'])
                self.assertIn('truncated migration graph', report['coverage_warning'])
                self.assertIn(report['coverage_warning'], report['warnings'])
            diffusion = study['groups'][0]['statistics']['diffusion']
            self.assertIn('conditional on the declared truncated graph', diffusion['label'])
            self.assertEqual(diffusion['channel_exclusions'], exclusions)
            with (root/'study'/'groups.csv').open() as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(json.loads(row['channel_exclusions']), exclusions)
            self.assertEqual(float(row['migration_cutoff_A']), 3.2)
            self.assertEqual(row['conditional_on_truncated_graph'], 'True')
            self.assertEqual(row['diffusion_label'], diffusion['label'])
            source = root/'truncated.json'
            truncated.write(source)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(['export-neb', '--catalog', str(source), '--charges', '2',
                                       '--channels', excluded, '--output', str(root/'neb')]), 0)
            worklist = json.loads((root/'neb'/'worklist.json').read_text())
            pending = json.loads((root/'neb'/'catalog.pending.json').read_text())
            for report in (worklist, pending['metadata']):
                self.assertEqual(report['source_channel_exclusions'], exclusions)
                self.assertEqual(report['source_migration_cutoff_A'], 3.2)
                self.assertFalse(report['source_channel_exclusions_applied'])
                self.assertIn('NOT automatically applied', report['source_channel_exclusions_handling'])
                self.assertEqual(len(report['source_catalog_sha256']), 64)
            self.assertEqual(pending['metadata']['channel_exclusions'], {})
            self.assertEqual(len(pending['channels']), len(self.lattice.channel_templates()))
            self.assertIn(excluded, {channel['id'] for channel in pending['channels']})
            mapping = json.loads((root/'neb'/excluded/'q2'/'mapping.json').read_text())
            self.assertEqual(mapping['source_channel_exclusion_reason'], exclusions[excluded])
            self.assertFalse(mapping['source_channel_exclusions_applied'])


if __name__ == '__main__':
    unittest.main()
