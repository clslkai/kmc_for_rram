"""Finite-reservoir conservation and independent frozen-field thermodynamics."""
from copy import deepcopy
from dataclasses import FrozenInstanceError
import json
import math
import unittest

import numpy as np

from dftkmc.exchange import ExchangeConfig, FiniteReservoir, ReservoirConfig, exchange_pair
from dftkmc.model import KB, ModelDomainError


class ExchangeConfigTests(unittest.TestCase):
    def test_frozen_complete_json_roundtrip_preserves_declared_approximation(self):
        config = ExchangeConfig(reference_basis=3, vacancy_formation_offset_eV=0.35,
                                top=ReservoirConfig(capacity=9, initial_oxygen=4, mu0_eV=-0.1))
        serialized = json.loads(json.dumps(config.to_dict(), allow_nan=False))
        self.assertEqual(ExchangeConfig.from_dict(serialized), config)
        self.assertEqual(ReservoirConfig.from_dict(config.top.to_dict()), config.top)
        self.assertEqual(config.calibration_status, 'illustrative')
        self.assertIn('SYNTHETIC', config.provenance)
        with self.assertRaises(FrozenInstanceError):
            config.top.capacity = 100
        self.assertFalse(ExchangeConfig().bottom.enabled)
        self.assertTrue(ExchangeConfig().top.enabled)

    def test_missing_unknown_or_false_calibration_is_rejected(self):
        alterations = [
            (lambda d: d.pop('reference_basis'), 'missing fields'),
            (lambda d: d.__setitem__('extra', 2), 'unknown fields'),
            (lambda d: d['top'].pop('mixing_entropy'), 'missing fields'),
            (lambda d: d['top'].__setitem__('rate_multiplier', 2), 'unknown fields'),
            (lambda d: d.__setitem__('calibration_status', 'dft'), 'only illustrative'),
            (lambda d: d.__setitem__('provenance', 'completed calculation'), 'SYNTHETIC'),
            (lambda d: d.__setitem__('provenance', ''), 'SYNTHETIC'),
        ]
        for alter, message in alterations:
            with self.subTest(message=message):
                data = ExchangeConfig().to_dict()
                alter(data)
                with self.assertRaisesRegex(ValueError, message):
                    ExchangeConfig.from_dict(data)

    def test_configuration_domains_are_explicit(self):
        bad = [dict(capacity=0), dict(capacity=2.0), dict(capacity=True),
               dict(capacity=3, initial_oxygen=4), dict(initial_oxygen=-1), dict(initial_oxygen=False),
               dict(mu0_eV=float('nan')), dict(kappa_eV=float('inf')), dict(barrier_eV=-0.1),
               dict(attempt_Hz=0), dict(attempt_Hz=True), dict(enabled=1), dict(mixing_entropy='true')]
        for arguments in bad:
            with self.subTest(arguments=arguments):
                with self.assertRaises(ValueError):
                    ReservoirConfig(**arguments)
        for basis in (-1, 8, 0.0, True):
            with self.assertRaises(ValueError):
                ExchangeConfig(reference_basis=basis)
        with self.assertRaises(ValueError):
            ExchangeConfig(vacancy_formation_offset_eV=float('nan'))
        with self.assertRaises(ValueError):
            ExchangeConfig(bottom={})


class FiniteReservoirTests(unittest.TestCase):
    def test_count_and_per_reservoir_transfer_ledger_are_exact(self):
        reservoir = FiniteReservoir(ReservoirConfig(capacity=5, initial_oxygen=2))
        for direction in (1, 1, -1, 1, -1, -1):
            reservoir.apply(direction)
            self.assertEqual(reservoir.count, reservoir.initial_count + reservoir.out_count - reservoir.in_count)
            self.assertGreaterEqual(reservoir.count, 0)
            self.assertLessEqual(reservoir.count, reservoir.config.capacity)
        self.assertEqual((reservoir.initial_count, reservoir.count, reservoir.out_count, reservoir.in_count), (2, 2, 3, 3))
        with self.assertRaises(AttributeError):
            reservoir.count = 9

    def test_capacity_walls_disabled_exchange_and_invalid_direction_do_not_mutate_state(self):
        reservoir = FiniteReservoir(ReservoirConfig(capacity=1))
        self.assertTrue(reservoir.available(1))
        self.assertFalse(reservoir.available(-1))
        with self.assertRaisesRegex(ValueError, 'empty'):
            reservoir.apply(-1)
        reservoir.apply(1)
        self.assertFalse(reservoir.available(1))
        self.assertTrue(reservoir.available(-1))
        with self.assertRaisesRegex(ValueError, 'full'):
            reservoir.apply(1)
        for direction in (0, 2, 1.0, True, 'out'):
            before = (reservoir.count, reservoir.out_count, reservoir.in_count)
            with self.assertRaises(ValueError):
                reservoir.apply(direction)
            self.assertEqual((reservoir.count, reservoir.out_count, reservoir.in_count), before)
        disabled = FiniteReservoir(ReservoirConfig(capacity=2, initial_oxygen=1, enabled=False))
        for direction in (-1, 1):
            self.assertFalse(disabled.available(direction))
            with self.assertRaisesRegex(ValueError, 'disabled'):
                disabled.apply(direction)
        self.assertEqual(disabled.count, 1)

    def test_ideal_mixing_entropy_matches_an_independent_binomial_count(self):
        temperature = 700.0
        reservoir = FiniteReservoir(ReservoirConfig(capacity=4, mu0_eV=0, kappa_eV=0))
        self.assertEqual(reservoir.free_energy_eV(temperature, count=0), 0.0)
        self.assertEqual(reservoir.free_energy_eV(temperature, count=4), 0.0)
        self.assertAlmostEqual(reservoir.free_energy_eV(temperature, count=2), -KB * temperature * math.log(6), 13)
        self.assertAlmostEqual(reservoir.creation_delta_eV(temperature, count_before=0), -KB * temperature * math.log(4), 13)
        self.assertAlmostEqual(reservoir.creation_delta_eV(temperature, count_before=3), KB * temperature * math.log(4), 13)

    def test_analytic_increment_matches_full_free_energy_at_every_capacity_state(self):
        for mixing in (True, False):
            reservoir = FiniteReservoir(ReservoirConfig(capacity=7, initial_oxygen=3, mu0_eV=-0.08,
                                                       kappa_eV=0.31, mixing_entropy=mixing))
            for m in range(7):
                difference = reservoir.free_energy_eV(650, count=m+1) - reservoir.free_energy_eV(650, count=m)
                self.assertAlmostEqual(reservoir.creation_delta_eV(650, count_before=m), difference, 13)
            self.assertAlmostEqual(reservoir.creation_delta_eV(650), reservoir.creation_delta_eV(650, 3))
        no_mixing = FiniteReservoir(ReservoirConfig(capacity=4, mu0_eV=0.15, kappa_eV=0.2, mixing_entropy=False))
        self.assertAlmostEqual(no_mixing.creation_delta_eV(300, 2), 0.25)
        self.assertEqual(no_mixing.free_energy_eV(300, 3), no_mixing.free_energy_eV(900, 3))

    def test_thermodynamic_queries_reject_bad_temperature_and_outside_counts(self):
        reservoir = FiniteReservoir(ReservoirConfig(capacity=4))
        for temperature in (0, -10, float('nan'), float('inf'), True, '300'):
            with self.subTest(temperature=temperature):
                with self.assertRaises(ValueError):
                    reservoir.free_energy_eV(temperature)
                with self.assertRaises(ValueError):
                    reservoir.creation_delta_eV(temperature)
        for count in (-1, 5, 1.0, True):
            with self.assertRaises(ValueError):
                reservoir.free_energy_eV(600, count=count)
        with self.assertRaises(ValueError):
            reservoir.creation_delta_eV(600, count_before=4)


class FrozenFieldExchangeTests(unittest.TestCase):
    def test_hand_chosen_electrochemical_work_and_shared_saddle(self):
        config = ReservoirConfig(capacity=4, mu0_eV=0.1, kappa_eV=0.0, mixing_entropy=False,
                                 barrier_eV=0.8, attempt_Hz=2e11)
        # Oxide +0.3, storage +0.1, electron/electrode term 2*(.6-.8)=-.4 eV.
        pair = exchange_pair(0.3, config, 1, 2, 0.6, 0.8, 600)
        self.assertAlmostEqual(pair['chemical_delta_eV'], 0.4, 13)
        self.assertAlmostEqual(pair['electrical_work_eV'], -0.4, 13)
        self.assertAlmostEqual(pair['delta_eV'], 0.0, 13)
        self.assertAlmostEqual(pair['out_barrier_eV'], 0.8, 13)
        self.assertAlmostEqual(pair['in_barrier_eV'], 0.8, 13)
        expected_rate = 2e11 * math.exp(-0.8 / (KB * 600))
        self.assertAlmostEqual(pair['out_rate_Hz'] / expected_rate, 1.0, 12)
        self.assertAlmostEqual(pair['in_rate_Hz'] / expected_rate, 1.0, 12)

    def test_gauge_shift_and_neutral_limit_are_exact_within_roundoff(self):
        config = ReservoirConfig(capacity=5)
        original = exchange_pair(0.4, config, 2, 2, 0.32, 0.75, 600)
        shifted = exchange_pair(0.4, config, 2, 2, 7.32, 7.75, 600)
        for field in original:
            np.testing.assert_allclose(original[field], shifted[field], rtol=1e-12, atol=1e-12)
        neutral = exchange_pair(0.4, config, 2, 0, 0.32, 0.75, 600)
        neutral_other = exchange_pair(0.4, config, 2, 0, -3.5, 2.2, 600)
        self.assertEqual(neutral, neutral_other)
        self.assertEqual(neutral['electrical_work_eV'], 0.0)

    def test_rate_ratio_contains_finite_capacity_entropy_once(self):
        config = ReservoirConfig(capacity=6, mu0_eV=0, kappa_eV=0, barrier_eV=0.7)
        pair = exchange_pair(0.0, config, 0, 0, 0, 0, 800)
        # Six empty storage slots versus one occupied slot after creation.
        self.assertAlmostEqual(pair['out_rate_Hz'] / pair['in_rate_Hz'], 6.0, 12)
        self.assertAlmostEqual(pair['out_barrier_eV'] - pair['in_barrier_eV'], pair['delta_eV'], 13)
        last = exchange_pair(0.0, config, 5, 0, 0, 0, 800)
        self.assertAlmostEqual(last['out_rate_Hz'] / last['in_rate_Hz'], 1.0 / 6.0, 13)

    def test_complete_small_state_space_has_the_independent_frozen_field_equilibrium(self):
        temperature, q = 900.0, 2
        config = ReservoirConfig(capacity=2, mu0_eV=0.03, kappa_eV=0.12, barrier_eV=1.0, attempt_Hz=1e7)
        site_creation = np.array([0.4, 0.55])
        pair_energy = -0.07
        phi, electrode = np.array([0.20, 0.45]), 0.6
        states = [np.array([bool(s & 1), bool(s & 2)]) for s in range(4)]
        generator = np.zeros((4, 4))
        log_weights = []
        for index, vacancy in enumerate(states):
            m = int(vacancy.sum())  # Total oxygen conservation fixes m=N_v here.
            internal_energy = float(site_creation @ vacancy) + pair_energy * int(vacancy.all())
            storage_energy = config.mu0_eV * m + config.kappa_eV * m * (m-1) / (2 * config.capacity)
            electrical_energy = q * (float(phi @ vacancy) - m * electrode)
            log_weights.append(math.log(math.comb(config.capacity, m)) - (internal_energy + storage_energy + electrical_energy)/(KB*temperature))
            for site in range(2):
                other = 1 - site
                oxide_delta = site_creation[site] + pair_energy * int(vacancy[other])
                canonical_count = m - 1 if vacancy[site] else m
                pair = exchange_pair(oxide_delta, config, canonical_count, q, phi[site], electrode, temperature)
                target = index ^ (1 << site)
                generator[index, target] = pair['in_rate_Hz'] if vacancy[site] else pair['out_rate_Hz']
        generator[np.diag_indices(4)] = -generator.sum(axis=1)
        weights = np.exp(np.asarray(log_weights) - max(log_weights))
        weights /= weights.sum()
        flux = weights[:, None] * generator
        np.testing.assert_allclose(flux, flux.T, rtol=5e-13, atol=1e-12)
        np.testing.assert_allclose(weights @ generator, np.zeros(4), atol=1e-11)
        # Nontrivial loop: create A, create B, fill A, fill B.
        loop = [0, 1, 3, 2, 0]
        affinity = sum(math.log(generator[a, b] / generator[b, a]) for a, b in zip(loop[:-1], loop[1:]))
        self.assertAlmostEqual(affinity, 0.0, 12)

    def test_negative_barriers_are_rejected_in_both_directions_without_clipping(self):
        config = ReservoirConfig(capacity=2, barrier_eV=0.1, mu0_eV=0, kappa_eV=0, mixing_entropy=False)
        for oxide_delta in (-0.3, 0.3):
            with self.subTest(oxide_delta=oxide_delta):
                with self.assertRaisesRegex(ModelDomainError, 'negative exchange barrier.*no clipping'):
                    exchange_pair(oxide_delta, config, 0, 2, 0, 0, 600)
        disabled = ReservoirConfig(enabled=False)
        with self.assertRaisesRegex(ValueError, 'disabled'):
            exchange_pair(0.4, disabled, 0, 2, 0, 0, 600)

    def test_bad_pair_input_and_full_canonical_count_fail_closed(self):
        base = dict(oxide_creation_energy_eV=0.4, reservoir_config=ReservoirConfig(capacity=2), count_before=0,
                    charge=2, site_potential_V=0.0, electrode_potential_V=0.0, temperature_K=600)
        for key, value in [('oxide_creation_energy_eV', None), ('oxide_creation_energy_eV', float('nan')),
                           ('charge', True), ('charge', 2.0), ('site_potential_V', float('inf')),
                           ('electrode_potential_V', False), ('count_before', 2), ('temperature_K', 0)]:
            with self.subTest(key=key, value=value):
                arguments = dict(base, **{key: value})
                with self.assertRaises(ValueError):
                    exchange_pair(**arguments)


if __name__ == '__main__':
    unittest.main()
