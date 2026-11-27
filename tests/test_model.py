import math
import unittest
import numpy as np
from vcm3d.model import Config, Simulator, KB

class ModelTests(unittest.TestCase):
    def sim(self, **kw):
        return Simulator(Config(shape=(3, 3, 4), **kw))

    def test_uniform_network_analytic(self):
        s = self.sim(); f = s.solve_fields(0.8)
        a=s.cfg.spacing_m; nx,ny,nz=s.cfg.shape
        expected=s.cfg.sigma_off_S_m*(nx*ny*a*a)/(nz*a)
        self.assertAlmostEqual(f["conductance_S"]/expected,1,places=10)
        z=np.arange(nz)+.5
        np.testing.assert_allclose(s.phi.reshape(s.cfg.shape)[1,1],.8*z/nz,rtol=1e-10)
        self.assertLess(abs(s.energy_residual_W),s.power_W*1e-10)
        self.assertLess(s.kcl_residual_A,abs(s.current_A)*1e-10)
        self.assertLess(abs(s.thermal_residual_W),s.power_W*1e-10)

    def test_compliance_and_power(self):
        s=self.sim(compliance_A=1e-8)
        s.vacancy[:]=True; s.interstitial[:]=True
        s.solve_fields(3)
        self.assertLess(abs(s.current_A),1.000000001e-8)
        self.assertLess(s.device_voltage_V,3)
        self.assertAlmostEqual(s.power_W/(s.current_A*s.device_voltage_V),1,places=10)
        self.assertTrue(np.all(s.temperature>=s.cfg.ambient_K-1e-9))

    def test_polarity_and_six_neighbors(self):
        s=self.sim()
        self.assertEqual(len(s.neighbors[0]),3)
        self.assertEqual(max(map(len,s.neighbors)),6)
        s.solve_fields(1); positive=s.phi.copy(); temp=s.temperature.copy()
        s.solve_fields(-1)
        np.testing.assert_allclose(s.phi,-positive)
        np.testing.assert_allclose(s.temperature,temp)
        self.assertLess(s.current_A,0)

    def test_pair_detailed_balance_and_conservation(self):
        s=self.sim(); s.solve_fields(1)
        e=next(e for e in s.events() if e.kind=="generate" and e.source != e.target)
        s.apply(e)
        reverse=next(r for r in s.events() if r.kind=="recombine" and (r.source,r.target)==(e.source,e.target))
        temp=(s.temperature[e.source]+s.temperature[e.target])/2
        self.assertAlmostEqual(math.log(e.rate_Hz/reverse.rate_Hz),-e.delta_eV/(KB*temp),places=10)
        self.assertAlmostEqual(reverse.delta_eV,-e.delta_eV)
        s.apply(reverse); self.assertEqual(s.vacancy.sum(),0); self.assertEqual(s.interstitial.sum(),0)

    def test_vacancy_hop_detailed_balance(self):
        s=self.sim(); s.initialize("random",.25); s.solve_fields(-.5)
        e=next(e for e in s.events() if e.kind=="vacancy_hop")
        s.apply(e)
        r=next(r for r in s.events() if r.kind==e.kind and r.source==e.target and r.target==e.source)
        self.assertAlmostEqual(r.delta_eV,-e.delta_eV)
        temp=(s.temperature[e.source]+s.temperature[e.target])/2
        self.assertAlmostEqual(math.log(e.rate_Hz/r.rate_Hz),-e.delta_eV/(KB*temp),places=10)

    def test_interstitial_hop_detailed_balance(self):
        s=self.sim(); s.initialize("random",.25); s.solve_fields(.5)
        e=next(e for e in s.events() if e.kind=="interstitial_hop")
        s.apply(e)
        r=next(r for r in s.events() if r.kind==e.kind and r.source==e.target and r.target==e.source)
        self.assertAlmostEqual(r.delta_eV,-e.delta_eV)

    def test_seed_and_bookkeeping(self):
        a=self.sim(seed=11); b=self.sim(seed=11)
        for s in [a,b]:
            s.initialize("random",.2)
            for _ in range(50): s.step(1)
            s.check_invariants()
        self.assertEqual(a.time_s,b.time_s)
        self.assertEqual(a.counts,b.counts)
        np.testing.assert_array_equal(a.vacancy,b.vacancy)
        self.assertEqual(a.event_count,50)

    def test_no_event_crosses_pulse_boundary(self):
        s=self.sim(); old=s.vacancy.copy()
        e=s.step(1,1e-30)
        self.assertIsNone(e); self.assertEqual(s.time_s,1e-30)
        self.assertEqual(s.event_count,0); np.testing.assert_array_equal(old,s.vacancy)

    def test_true_3d_connectivity(self):
        s=self.sim(); self.assertFalse(s.spanning())
        s.initialize("filament"); self.assertTrue(s.spanning())
        s.vacancy.reshape(s.cfg.shape)[1,1,2]=False
        self.assertFalse(s.spanning())

    def test_zero_bias_no_heating(self):
        s=self.sim(); s.initialize("random"); s.solve_fields(0)
        np.testing.assert_array_equal(s.phi,0)
        np.testing.assert_allclose(s.temperature,s.cfg.ambient_K)
        self.assertEqual(s.power_W,0)

    def test_exact_waiting_time_and_rate_weighting(self):
        from unittest.mock import patch
        from vcm3d.model import Event
        s = self.sim()
        catalog = [Event("generate", 0, 0, 0, 1), Event("generate", 1, 1, 0, 3)]
        class ControlledRng:
            def __init__(self): self.values = iter([math.exp(-2), 0.6])
            def random(self): return next(self.values)
        s.rng = ControlledRng()
        with patch.object(s, "solve_fields"), patch.object(s, "events", return_value=catalog):
            chosen = s.step(1)
        self.assertEqual(chosen.source, 1)
        self.assertAlmostEqual(s.time_s, 0.5)  # -log(exp(-2)) / (1+3)
        self.assertEqual(s.vacancy.sum(), 1)

    def test_frozen_read_conductance_independent_of_write_polarity(self):
        s=self.sim(); s.initialize("filament")
        a=s.snapshot(2); b=s.snapshot(-2)
        self.assertEqual(a["conductance_S"], b["conductance_S"])
        self.assertAlmostEqual(a["resistance_ohm"] * a["conductance_S"], 1)
        self.assertEqual(a["frozen_read_0p1V_current_A"], b["frozen_read_0p1V_current_A"])

    def test_invalid_config(self):
        for kw in ({"shape":(1,2,3)}, {"ambient_K":0}, {"sigma_off_S_m":-1}, {"pair_energy_eV":float("nan")}):
            with self.assertRaises(ValueError): Config(**kw)

if __name__=="__main__": unittest.main()
