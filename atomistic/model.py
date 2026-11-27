"""Monoclinic HfO2 oxygen-site KMC with finite top oxygen reservoir.
All crystal distances are Angstrom; energies eV; time s; electrical quantities SI.
The candidate hop graph is geometric, NOT an NEB-validated migration network.
"""
from dataclasses import dataclass, asdict
from itertools import product
import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.sparse.csgraph import connected_components
from vcm3d.model import KB, Event

# Pathak et al., Acta Cryst C76 (2020), combined refinement, ideal Hf occupancy 1.
# Full verified publisher CIF and provenance JSON are included in repository.
SOURCE = 'https://doi.org/10.1107/S2053229620013960'
CELL = (5.1145, 5.1682, 5.2900, 99.211)
ASYMMETRIC = {'Hf': (.2749,.0402,.2084), 'O3': (.0670,.3273,.3428), 'O4': (.4512,.7558,.4840)}

def equivalents(p):
    """P21/c #14, unique axis b, conventional origin choice."""
    x,y,z = p
    return np.mod([(x,y,z),(-x,y+.5,-z+.5),(-x,-y,-z),(x,-y+.5,z+.5)],1)

@dataclass(frozen=True)
class AtomConfig:
    repeats: tuple = (3,3,5)
    seed: int = 17
    temperature_K: float = 600.
    attempt_Hz: float = 1e12
    hop_cutoff_A: float = 3.2
    bond_cutoff_A: float = 2.6
    contact_depth_A: float = 2.0
    migration_barrier_eV: float = .70
    exchange_barrier_eV: float = .75
    site_O3_eV: float = .20
    site_O4_eV: float = .30
    binding_eV: float = .06
    reservoir_stiffness_eV: float = .10
    charge_e: float = 2.
    g_off_S: float = 1e-12
    g_on_S: float = 1e-6
    compliance_A: float = 1e-5
    def __post_init__(self):
        if len(self.repeats)!=3 or any(type(n) is not int or n<2 for n in self.repeats):
            raise ValueError('repeats must contain three integer dimensions >=2')
        for k,v in asdict(self).items():
            if k in ('seed','repeats'): continue
            if not math.isfinite(v): raise ValueError(k+' must be finite')
        for k in ('temperature_K','attempt_Hz','hop_cutoff_A','bond_cutoff_A','contact_depth_A','g_off_S','g_on_S','compliance_A'):
            if getattr(self,k)<=0: raise ValueError(k+' must be positive')
        for k in ('migration_barrier_eV','exchange_barrier_eV','binding_eV','reservoir_stiffness_eV','charge_e'):
            if getattr(self,k)<0: raise ValueError(k+' must be nonnegative')
        if self.g_on_S<self.g_off_S: raise ValueError('g_on must be >= g_off')

class Crystal:
    def __init__(self,cfg=AtomConfig()):
        self.cfg=cfg
        a,b,c,beta=CELL; beta=np.deg2rad(beta)
        self.unit=np.array([[a,0,0],[0,b,0],[c*np.cos(beta),0,c*np.sin(beta)]])
        self.cell=self.unit*np.asarray(cfg.repeats)[:,None]
        positions={k:[] for k in ASYMMETRIC}
        for shift in product(*[range(n) for n in cfg.repeats]):
            for k,p in ASYMMETRIC.items():
                positions[k].extend((equivalents(p)+shift)@self.unit)
        self.hf=np.array(positions['Hf'])
        self.oxygen=np.concatenate([positions['O3'],positions['O4']])
        self.kind=np.r_[np.zeros(len(positions['O3']),int),np.ones(len(positions['O4']),int)]
        self.n=len(self.oxygen)
        self.height=self.cell[2,2]
        # a and b are in-plane vectors. Never wrap through top/bottom contacts.
        delta=self.oxygen[None,:,:]-self.oxygen[:,None,:]
        for axis in (0,1):
            length=self.cell[axis,axis]
            delta[:,:,axis]-=np.rint(delta[:,:,axis]/length)*length
        distance=np.linalg.norm(delta,axis=2)
        ii,jj=np.where(np.triu((distance<=cfg.hop_cutoff_A)&(distance>1e-8),1))
        self.edges=np.array([ii,jj]).T
        self.distance=distance[ii,jj]
        self.displacements=delta[ii,jj]
        self.neighbors=[[] for _ in range(self.n)]
        for i,j in self.edges:
            self.neighbors[i].append(int(j)); self.neighbors[j].append(int(i))
        graph=coo_matrix((np.ones(len(ii)*2),(np.r_[ii,jj],np.r_[jj,ii])),shape=(self.n,self.n))
        if connected_components(graph)[0]!=1: raise ValueError('hop cutoff creates disconnected oxygen graph')
        self.bottom=np.flatnonzero(self.oxygen[:,2]<cfg.contact_depth_A)
        self.top=np.flatnonzero(self.oxygen[:,2]>self.height-cfg.contact_depth_A)
        if not len(self.bottom) or not len(self.top) or np.intersect1d(self.bottom,self.top).size:
            raise ValueError('invalid electrode contact depth')
        self.reservoir_capacity=self.n
        # Physical bulk Hf coordination using periodic images in all directions.
        self.bulk_coordination=[]
        for k in ('O3','O4'):
            o=equivalents(ASYMMETRIC[k])[0]@self.unit
            hf=np.concatenate([(equivalents(ASYMMETRIC['Hf'])+s)@self.unit for s in product((-1,0,1),repeat=3)])
            self.bulk_coordination.append(int((np.linalg.norm(hf-o,axis=1)<cfg.bond_cutoff_A).sum()))

class AtomSimulator:
    def __init__(self,cfg=AtomConfig()):
        self.cfg=cfg; self.crystal=Crystal(cfg); self.n=self.crystal.n
        self.rng=np.random.default_rng(cfg.seed)
        self.vacancy=np.zeros(self.n,bool); self.reservoir=0
        self.time_s=0.; self.event_count=0
        self.phi=np.zeros(self.n); self.voltage=0.; self.conductance=0.; self.current=0.
        self.site_energy=np.where(self.crystal.kind==0,cfg.site_O3_eV,cfg.site_O4_eV)
        self.counts={'hop':0,'oxygen_out':0,'oxygen_in':0}

    def initialize(self,fraction=.03):
        if not 0<=fraction<=1: raise ValueError('fraction outside [0,1]')
        if self.event_count or self.time_s: raise ValueError('already running')
        self.vacancy[:]=False
        n=round(fraction*self.n)
        self.vacancy[self.rng.choice(self.n,n,replace=False)]=True
        self.reservoir=n
        self.check()

    def check(self):
        if int(self.vacancy.sum())!=self.reservoir or not 0<=self.reservoir<=self.crystal.reservoir_capacity:
            raise RuntimeError('oxygen conservation failed')

    def reservoir_energy(self,n):
        # Finite well-mixed reservoir free energy, explicit illustrative convex term.
        return self.cfg.reservoir_stiffness_eV*n*n/(2*self.crystal.reservoir_capacity)

    def solve(self,bias):
        if not math.isfinite(bias): raise ValueError('nonfinite voltage')
        c=self.cfg; cr=self.crystal; i,j=cr.edges.T
        # Effective electronic edge model, NOT a derived HfO2 electronic Hamiltonian.
        g=np.where(self.vacancy[i]&self.vacancy[j],c.g_on_S,c.g_off_S)*3/cr.distance
        gb=np.where(self.vacancy[cr.bottom],c.g_on_S,c.g_off_S)*1.5/cr.oxygen[cr.bottom,2]
        gt=np.where(self.vacancy[cr.top],c.g_on_S,c.g_off_S)*1.5/(cr.height-cr.oxygen[cr.top,2])
        diagonal=np.bincount(np.r_[i,j],weights=np.r_[g,g],minlength=self.n)
        diagonal[cr.bottom]+=gb; diagonal[cr.top]+=gt
        matrix=coo_matrix((np.r_[-g,-g,diagonal],(np.r_[i,j,np.arange(self.n)],np.r_[j,i,np.arange(self.n)])),shape=(self.n,self.n)).tocsr()
        rhs=np.zeros(self.n); rhs[cr.top]=gt
        unit=spsolve(matrix,rhs)
        # Bottom-terminal current avoids cancellation 1-unit near conducting top clusters.
        self.conductance=float(gb@unit[cr.bottom])
        if not self.conductance>0 or not math.isfinite(self.conductance): raise RuntimeError('bad network conductance')
        self.voltage=math.copysign(min(abs(bias),c.compliance_A/self.conductance),bias)
        self.phi=unit*self.voltage; self.current=self.conductance*self.voltage
        self.power=float(np.sum(g*(self.phi[i]-self.phi[j])**2)+np.sum(gb*self.phi[cr.bottom]**2)+np.sum(gt*(self.voltage-self.phi[cr.top])**2))
        self.power_error=self.power-self.current*self.voltage
        self.kcl_error=float(np.max(abs(matrix@self.phi-rhs*self.voltage)))
        return self.conductance

    def rate(self,delta,barrier):
        return self.cfg.attempt_Hz*math.exp(-(barrier+max(0,delta))/(KB*self.cfg.temperature_K))

    def events(self):
        c=self.cfg; cr=self.crystal; v=self.vacancy
        nn=np.array([v[n].sum() for n in cr.neighbors],int)
        ev=[]
        for i in np.flatnonzero(v):
            for j in cr.neighbors[i]:
                if v[j]: continue
                de=self.site_energy[j]-self.site_energy[i]+c.charge_e*(self.phi[j]-self.phi[i])-c.binding_eV*(nn[j]-1-nn[i])
                ev.append(Event('hop',int(i),j,float(de),self.rate(de,c.migration_barrier_eV)))
        for i in cr.top:
            if not v[i] and self.reservoir<cr.reservoir_capacity:
                de=self.site_energy[i]-c.binding_eV*nn[i]+c.charge_e*(self.phi[i]-self.voltage)+self.reservoir_energy(self.reservoir+1)-self.reservoir_energy(self.reservoir)
                ev.append(Event('oxygen_out',int(i),-1,float(de),self.rate(de,c.exchange_barrier_eV)))
            elif v[i] and self.reservoir>0:
                de=-self.site_energy[i]+c.binding_eV*nn[i]-c.charge_e*(self.phi[i]-self.voltage)+self.reservoir_energy(self.reservoir-1)-self.reservoir_energy(self.reservoir)
                ev.append(Event('oxygen_in',int(i),-1,float(de),self.rate(de,c.exchange_barrier_eV)))
        return [e for e in ev if e.rate_Hz>0]

    def apply(self,e):
        i,j=e.source,e.target; v=self.vacancy
        if e.kind=='hop':
            if not v[i] or v[j] or j not in self.crystal.neighbors[i]: raise ValueError('invalid hop')
            v[i]=False; v[j]=True
        elif e.kind=='oxygen_out':
            if v[i] or i not in self.crystal.top or self.reservoir>=self.n: raise ValueError('invalid exchange')
            v[i]=True; self.reservoir+=1
        elif e.kind=='oxygen_in':
            if not v[i] or i not in self.crystal.top or self.reservoir<=0: raise ValueError('invalid exchange')
            v[i]=False; self.reservoir-=1
        else: raise ValueError('unknown event')
        self.event_count+=1; self.counts[e.kind]+=1; self.check()

    def step(self,bias,end):
        if end<self.time_s: raise ValueError('time boundary in past')
        self.solve(bias); ev=self.events(); total=sum(e.rate_Hz for e in ev)
        if not math.isfinite(total): raise RuntimeError('invalid total rate')
        if total==0: self.time_s=end; return None
        dt=-math.log(max(self.rng.random(),np.finfo(float).tiny))/total
        if self.time_s+dt>=end: self.time_s=end; return None
        x=self.rng.random()*total
        chosen=ev[-1]
        for e in ev:
            x-=e.rate_Hz
            if x<0: chosen=e; break
        self.time_s+=dt; self.apply(chosen); return chosen

    def spanning(self):
        cr=self.crystal; seen=set(int(i) for i in cr.bottom if self.vacancy[i]); stack=list(seen); target=set(cr.top)
        while stack:
            i=stack.pop()
            if i in target: return True
            for j in cr.neighbors[i]:
                if self.vacancy[j] and j not in seen: seen.add(j); stack.append(j)
        return False

    def snapshot(self,bias):
        self.solve(bias)
        return dict(time_s=self.time_s,events=self.event_count,applied_V=bias,device_V=self.voltage,current_A=self.current,conductance_S=self.conductance,
                    vacancies=int(self.vacancy.sum()),reservoir_oxygen=self.reservoir,oxide_oxygen=self.n-self.reservoir,total_oxygen=self.n,
                    spanning=self.spanning(),power_W=self.power,kcl_error_A=self.kcl_error,power_error_W=self.power_error,temperature_K=self.cfg.temperature_K)
