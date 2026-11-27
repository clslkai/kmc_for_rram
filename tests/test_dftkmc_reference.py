"""Independent analytic checks, using solvable synthetic periodic walks."""
import copy
import math
import unittest

import numpy as np

from dftkmc.catalog import make_catalog
from dftkmc.geometry import Lattice
from dftkmc.model import KB
from dftkmc.reference import finite_time_tracer_reference, periodic_tracer_reference


def walk_catalog(*, transverse=4.0, rotation=None, axis_rates=(4., 4., 4.)):
    """Eight equally spaced sites per cell; nearest hops form a known walk.

    transverse=4 gives independent one-dimensional chains. transverse=1 gives
    an ordinary simple-cubic walk represented with eight internal x positions;
    y and z hops preserve the basis ID and therefore test self-basis events.
    """
    rotation = np.eye(3) if rotation is None else np.asarray(rotation)
    reference = {
        'cell_A': (np.diag([8., transverse, transverse]) @ rotation).tolist(),
        'oxygen_fractional': [[i/8., 0., 0.] for i in range(8)],
        'hf_fractional': [[(2*i+.5)/8., .5, .5] for i in range(4)],
        'oxygen_kinds': ['O3']*4+['O4']*4,
        'provenance': 'Synthetic solvable periodic-walk test geometry, not a physical HfO2 structure.',
    }
    lattice = Lattice(reference=reference, migration_cutoff_A=1.01, interaction_cutoff_A=1.01)
    catalog = make_catalog(lattice, charges=(2,))
    catalog.data['metadata']['material'] = 'synthetic solvable periodic random walk'
    catalog.data['site_energies_eV']['2'] = [0.]*8
    model = catalog.data['charge_models']['2']
    model['pair_energies_eV'] = dict.fromkeys(model['pair_energies_eV'], 0.)
    model['environment_saddle_eV_per_neighbor'] = 0.
    model['environment_mode'] = 'dilute_only'
    for channel in catalog.data['channels']:
        original_direction = np.asarray(channel['displacement_A']) @ rotation.T
        axis = int(np.argmax(np.abs(original_direction)))
        channel['data']['2'].update(saddle_eV=0., attempt_Hz=float(axis_rates[axis]),
                                    provenance='Exact synthetic test rates; zero barrier, chosen Poisson frequency.')
    return catalog


class PeriodicReferenceTests(unittest.TestCase):
    def test_uniform_chain_has_known_poisson_diffusion_and_jump_rate(self):
        catalog = walk_catalog(axis_rates=(4.,4.,4.))
        result = periodic_tracer_reference(catalog, 2, 600.)
        np.testing.assert_allclose(result['tensor_m2_s'], np.diag([4.,0.,0.])*1e-20, atol=1e-33)
        np.testing.assert_allclose(result['basis_probabilities'], np.ones(8)/8., atol=1e-15)
        self.assertAlmostEqual(result['mean_jump_rate_Hz'], 8.)
        np.testing.assert_allclose(result['drift_A_s'], 0., atol=1e-14)
        self.assertLess(result['residual_diagnostics']['corrector_equation_relative_residual'], 1e-12)
        self.assertEqual(result['calibration_status'], 'illustrative')

    def test_rotated_cubic_walk_includes_self_basis_jumps_and_covariance(self):
        rotation = np.array([[1/math.sqrt(2),1/math.sqrt(2),0.],
                             [-1/math.sqrt(6),1/math.sqrt(6),2/math.sqrt(6)],
                             [1/math.sqrt(3),-1/math.sqrt(3),1/math.sqrt(3)]])
        rates = (2.,3.,5.)
        catalog = walk_catalog(transverse=1., rotation=rotation, axis_rates=rates)
        expected_D_A2_s = rotation.T @ np.diag(rates) @ rotation
        result = periodic_tracer_reference(catalog, 2, 600.)
        np.testing.assert_allclose(result['tensor_m2_s'], expected_D_A2_s*1e-20, rtol=1e-12, atol=1e-32)
        self.assertAlmostEqual(result['mean_jump_rate_Hz'], 2*sum(rates))
        self.assertEqual(result['residual_diagnostics']['self_basis_directed_event_count'], 32)
        for duration in (.003, .7):
            finite = finite_time_tracer_reference(catalog, 2, 600., duration)
            # Difference of independent +/- Poisson jumps: Var(dx_a)=2*k_a*t.
            np.testing.assert_allclose(finite['covariance_A2'], 2*duration*expected_D_A2_s, rtol=1e-11, atol=1e-12)
            np.testing.assert_allclose(finite['tensor_m2_s'], expected_D_A2_s*1e-20, rtol=1e-11, atol=1e-32)
            np.testing.assert_allclose(finite['mean_displacement_A'], 0., atol=1e-12)
            self.assertAlmostEqual(finite['expected_hop_count'], 2*sum(rates)*duration)
            self.assertEqual(finite['residual_diagnostics']['moment_system_dimension'], 80)

    def test_alternating_bond_chain_has_harmonic_long_time_and_known_finite_moments(self):
        catalog = walk_catalog()
        fast, slow = 9., 1.
        for channel in catalog.data['channels']:
            dx = channel['displacement_A'][0]
            left = channel['from_basis'] if dx > 0 else channel['to_basis']
            channel['data']['2']['attempt_Hz'] = fast if left % 2 == 0 else slow
        result = periodic_tracer_reference(catalog, 2, 500.)
        # Resistors in a 1-D periodic chain add in series, giving the harmonic
        # mean bond rate, not the arithmetic mean instantaneous jump variance.
        expected_D = 2*fast*slow/(fast+slow)
        self.assertAlmostEqual(result['tensor_m2_s'][0][0]/1e-20, expected_D, places=12)
        self.assertGreater(np.linalg.norm(result['corrector_A']), .1)
        self.assertAlmostEqual(result['mean_jump_rate_Hz'], fast+slow)
        duration = .07
        total, contrast = fast+slow, fast-slow
        # Closed two-parity solution: the displacement-bias correlation decays
        # as exp(-2*(fast+slow)*t). This is independent of the 80-state solver.
        expected_variance = (4*fast*slow/total*duration
                             + contrast**2/(2*total**2)*(1-math.exp(-2*total*duration)))
        finite = finite_time_tracer_reference(catalog, 2, 500., duration)
        self.assertAlmostEqual(finite['covariance_A2'][0][0], expected_variance, places=12)
        self.assertAlmostEqual(finite['tensor_m2_s'][0][0]/1e-20, expected_variance/(2*duration), places=12)
        self.assertGreater(finite['tensor_m2_s'][0][0], result['tensor_m2_s'][0][0])

    def test_boltzmann_basis_and_branch_energy_zero(self):
        catalog = walk_catalog()
        temperature = 700.
        energy = KB*temperature*math.log(3.)
        catalog.data['site_energies_eV']['2'] = [0. if i%2 == 0 else energy for i in range(8)]
        for channel in catalog.data['channels']:
            channel['data']['2'].update(saddle_eV=energy, attempt_Hz=6.)
        result = periodic_tracer_reference(catalog, 2, temperature)
        # Outgoing rate per direction is 2 on low-energy sites and 6 on
        # high-energy sites; total equilibrium mass is 3/4 and 1/4 respectively.
        np.testing.assert_allclose(result['basis_probabilities'], [.75/4 if i%2 == 0 else .25/4 for i in range(8)], atol=1e-14)
        self.assertAlmostEqual(result['mean_jump_rate_Hz'], 6., places=12)
        self.assertAlmostEqual(result['tensor_m2_s'][0][0]/1e-20, 3., places=12)
        self.assertLess(result['residual_diagnostics']['detailed_balance_max_abs_log_flux_ratio'], 1e-12)
        shifted = copy.deepcopy(catalog)
        shifted.data['site_energies_eV']['2'] = [value+2. for value in catalog.data['site_energies_eV']['2']]
        for channel in shifted.data['channels']:
            channel['data']['2']['saddle_eV'] += 2.
        other = periodic_tracer_reference(shifted, 2, temperature)
        np.testing.assert_allclose(other['tensor_m2_s'], result['tensor_m2_s'], rtol=1e-12, atol=1e-33)
        np.testing.assert_allclose(other['basis_probabilities'], result['basis_probabilities'], rtol=1e-12)

    def test_disconnected_basis_is_rejected_in_both_references(self):
        catalog = walk_catalog()
        removed = [channel for channel in catalog.data['channels']
                   if {channel['from_basis'],channel['to_basis']} in ({0,7},{3,4})]
        removed_ids = {channel['id'] for channel in removed}
        catalog.data['channels'] = [channel for channel in catalog.data['channels'] if channel['id'] not in removed_ids]
        catalog.data['metadata']['channel_exclusions'] = dict.fromkeys(removed_ids, 'Synthetic disconnected-component test.')
        for function, arguments in ((periodic_tracer_reference, (catalog,2,600.)),
                                    (finite_time_tracer_reference, (catalog,2,600.,1.))):
            with self.assertRaisesRegex(ValueError, 'Disconnected basis transition graph'):
                function(*arguments)

    def test_connected_truncated_graph_is_conditioned_and_can_have_zero_long_time_D(self):
        catalog = walk_catalog()
        removed = next(channel for channel in catalog.data['channels'] if any(channel['cell_shift']))
        catalog.data['channels'].remove(removed)
        reason = 'Synthetic test: remove the intercell edge, leaving disconnected finite cages in real space.'
        catalog.data['metadata']['channel_exclusions'] = {removed['id']:reason}
        result = periodic_tracer_reference(catalog, 2, 600.)
        np.testing.assert_allclose(result['tensor_m2_s'], np.zeros((3,3)), atol=1e-45)
        self.assertTrue(result['conditional_on_truncated_graph'])
        self.assertEqual(result['channel_exclusions'], {removed['id']:reason})
        self.assertIn('conditional on the declared truncated graph', result['label'])
        finite = finite_time_tracer_reference(catalog, 2, 600., 100.)
        # At long time the starting and finishing sites are independent uniform
        # positions within an eight-site cage: Var(Xt-X0)=2 Var(uniform(0,...,7)).
        self.assertAlmostEqual(finite['covariance_A2'][0][0], 2*(8**2-1)/12., places=8)
        self.assertTrue(finite['conditional_on_truncated_graph'])
        self.assertGreater(finite['tensor_m2_s'][0][0], 0.)

    def test_pending_and_invalid_duration_are_rejected(self):
        catalog = walk_catalog()
        with self.assertRaisesRegex(ValueError, 'duration_s'):
            finite_time_tracer_reference(catalog, 2, 600., 0.)
        with self.assertRaisesRegex(ValueError, 'temperature_K'):
            periodic_tracer_reference(catalog, 2, -1.)
        catalog.data['metadata']['calibration_status'] = 'pending'
        with self.assertRaisesRegex(ValueError, 'pending'):
            periodic_tracer_reference(catalog, 2, 600.)


if __name__ == '__main__':
    unittest.main()
