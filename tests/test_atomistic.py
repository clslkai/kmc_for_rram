import unittest,math,tempfile,json,subprocess,sys,shutil
from pathlib import Path
import numpy as np
from atomistic.model import AtomConfig,Crystal,AtomSimulator,KB

class AtomisticTests(unittest.TestCase):
    def setUp(self): self.cfg=AtomConfig(repeats=(2,2,3))
    def test_crystal(self):
        cr=Crystal(self.cfg)
        self.assertEqual(len(cr.hf),4*12);self.assertEqual(cr.n,8*12)
        self.assertEqual(cr.bulk_coordination,[3,4]);self.assertEqual(len(np.unique(np.round(cr.oxygen,8),axis=0)),cr.n)
        self.assertTrue(np.all(cr.distance<=self.cfg.hop_cutoff_A));self.assertGreater(cr.distance.min(),2.)
        # Candidate graph has real nonuniform bond lengths, not cubic nearest neighbors.
        self.assertGreater(len(np.unique(np.round(cr.distance,3))),3)
        for i,j in cr.edges:self.assertIn(i,cr.neighbors[j]);self.assertIn(j,cr.neighbors[i])
    def test_periodic_lateral_not_vertical(self):
        cr=Crystal(self.cfg);i,j=cr.edges.T
        raw=cr.oxygen[j]-cr.oxygen[i]
        self.assertTrue(np.any(np.linalg.norm(raw,axis=1)>self.cfg.hop_cutoff_A))
        np.testing.assert_allclose(raw[:,2],cr.displacements[:,2])
    def test_mass_and_boundary(self):
        s=AtomSimulator(self.cfg);s.initialize(.1)
        for _ in range(30):
            s.step(1,1e-3);s.check()
        self.assertEqual(s.n-int(s.vacancy.sum())+s.reservoir,s.n)
        for e in s.events():
            if e.kind!='hop':self.assertIn(e.source,s.crystal.top)
    def test_hop_reverse_detailed_balance(self):
        s=AtomSimulator(self.cfg);s.initialize(.1);s.solve(.2)
        e=next(x for x in s.events() if x.kind=='hop');s.apply(e)
        rev=next(x for x in s.events() if x.kind=='hop' and x.source==e.target and x.target==e.source)
        self.assertAlmostEqual(e.delta_eV,-rev.delta_eV,12)
        self.assertAlmostEqual(math.log(e.rate_Hz/rev.rate_Hz),-e.delta_eV/(KB*self.cfg.temperature_K),10)
    def test_exchange_reverse_detailed_balance(self):
        s=AtomSimulator(self.cfg);s.initialize(0);s.solve(.5)
        e=next(x for x in s.events() if x.kind=='oxygen_out');s.apply(e)
        rev=next(x for x in s.events() if x.kind=='oxygen_in' and x.source==e.source)
        self.assertAlmostEqual(e.delta_eV,-rev.delta_eV,12)
        self.assertAlmostEqual(math.log(e.rate_Hz/rev.rate_Hz),-e.delta_eV/(KB*self.cfg.temperature_K),10)
    def test_fields(self):
        s=AtomSimulator(self.cfg);s.initialize(.1);s.solve(1.)
        self.assertTrue(np.all(s.phi>=-1e-10));self.assertTrue(np.all(s.phi<=s.voltage+1e-10))
        self.assertLess(s.kcl_error/max(abs(s.current),1e-30),1e-6)
        self.assertLess(abs(s.power_error)/max(s.power,1e-30),1e-6)
        g=s.conductance;s.solve(-1.);self.assertAlmostEqual(s.conductance,g,18);self.assertLess(s.current,0)
    def test_compliance(self):
        s=AtomSimulator(AtomConfig(repeats=(2,2,3),compliance_A=1e-10));s.initialize(1.);s.solve(1.)
        self.assertLess(s.voltage,1.);self.assertLessEqual(abs(s.current),1.000001e-10)
    def test_boundaries_reproducible(self):
        states=[]
        for _ in range(2):
            s=AtomSimulator(self.cfg);s.initialize(.1)
            while s.time_s<1e-7:s.step(1,1e-7)
            self.assertEqual(s.time_s,1e-7);states.append(s.vacancy)
        np.testing.assert_array_equal(*states)
    def test_spanning(self):
        s=AtomSimulator(self.cfg);s.initialize(0);self.assertFalse(s.spanning())
        s.initialize(1.);self.assertTrue(s.spanning())
    def test_cli_config_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);conf=tmp/'config.json';out=tmp/'out'
            conf.write_text(json.dumps(dict(repeats=[2,2,2],seed=123,temperature_K=310.,bond_cutoff_A=2.5)))
            subprocess.run([sys.executable,'-m','atomistic','--config',str(conf),'--pulse','0','1e-12','--output',str(out),'--no-plot'],check=True,capture_output=True)
            data=json.loads((out/'summary.json').read_text())
            self.assertEqual(data['config']['repeats'],[2,2,2]);self.assertEqual(data['config']['seed'],123)
            self.assertEqual(data['config']['temperature_K'],310.)
            self.assertIn('"temperature_K": 310.0',(out/'viewer.html').read_text())
            self.assertIn('"bond_cutoff_A": 2.5',(out/'viewer.html').read_text())
            if shutil.which('node'):
                script=tmp/'viewer.js';script.write_text((out/'viewer.html').read_text().split('<script>')[1].split('</script>')[0])
                subprocess.run(['node','--check',str(script)],check=True,capture_output=True)
    def test_validation(self):
        with self.assertRaises(ValueError):AtomConfig(repeats=(1,2,3))
        with self.assertRaises(ValueError):AtomConfig(temperature_K=float('nan'))
        with self.assertRaises(ValueError):Crystal(AtomConfig(repeats=(2,2,3),hop_cutoff_A=1.))

if __name__=='__main__':unittest.main()
