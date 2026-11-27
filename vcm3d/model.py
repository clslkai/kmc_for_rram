"""3D closed-oxygen VCM reference model. Energies eV; SI otherwise.

A regular cell has an oxygen site and a coarse-grained interstitial site.
The electrical network is quasi-neutral and does not solve space-charge
Poisson electrostatics. See MODEL.md before interpreting simulation results.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.ndimage import label

KB = 8.617333262145e-5  # eV/K

@dataclass(frozen=True)
class Config:
    shape: tuple[int, int, int] = (7, 7, 12)  # x,y,z; z is electrode normal
    spacing_m: float = 0.5e-9
    ambient_K: float = 500.0
    attempt_Hz: float = 1e12
    vacancy_charge_e: float = 2.0
    pair_energy_eV: float = 0.8
    pair_barrier_eV: float = 0.45
    vacancy_barrier_eV: float = 0.65
    interstitial_barrier_eV: float = 0.70
    vacancy_binding_eV: float = 0.08
    sigma_off_S_m: float = 1e-3
    sigma_on_S_m: float = 1e3
    thermal_W_mK: float = 1.0
    compliance_A: float = 10e-6
    max_temperature_K: float = 1500.0
    seed: int = 7

    def __post_init__(self):
        if len(self.shape) != 3 or any(type(n) is not int or n < 2 for n in self.shape):
            raise ValueError("shape must contain three integer dimensions >= 2")
        positive = ("spacing_m", "ambient_K", "attempt_Hz", "sigma_off_S_m",
                    "sigma_on_S_m", "thermal_W_mK", "compliance_A", "max_temperature_K")
        for name in positive:
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")
        for name in ("pair_energy_eV", "pair_barrier_eV", "vacancy_barrier_eV",
                     "interstitial_barrier_eV", "vacancy_binding_eV", "vacancy_charge_e"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) < 0:
                raise ValueError(f"{name} must be nonnegative and finite")
        if self.max_temperature_K < self.ambient_K:
            raise ValueError("max_temperature_K must be >= ambient_K")
        if self.sigma_on_S_m < self.sigma_off_S_m:
            raise ValueError("sigma_on_S_m must be >= sigma_off_S_m")

@dataclass(frozen=True)
class Event:
    kind: str
    source: int
    target: int
    delta_eV: float
    rate_Hz: float

class Simulator:
    def __init__(self, config: Config = Config()):
        self.cfg = config
        self.rng = np.random.default_rng(config.seed)
        self.n = int(np.prod(config.shape))
        self.vacancy = np.zeros(self.n, dtype=bool)
        self.interstitial = np.zeros(self.n, dtype=bool)
        ids = np.arange(self.n).reshape(config.shape)
        edges = []
        for axis in range(3):
            a = [slice(None)] * 3; b = a.copy()
            a[axis] = slice(None, -1); b[axis] = slice(1, None)
            edges.extend(zip(ids[tuple(a)].ravel(), ids[tuple(b)].ravel()))
        self.edges = np.asarray(edges, dtype=int)
        self.neighbors = [[] for _ in range(self.n)]
        for i, j in self.edges:
            self.neighbors[i].append(int(j)); self.neighbors[j].append(int(i))
        self.bottom = ids[:, :, 0].ravel()
        self.top = ids[:, :, -1].ravel()
        self.time_s = 0.0
        self.event_count = 0
        self.phi = np.zeros(self.n)
        self.temperature = np.full(self.n, config.ambient_K)
        self.conductance_S = 0.0
        self.current_A = 0.0
        self.device_voltage_V = 0.0
        self.power_W = 0.0
        self.energy_residual_W = 0.0
        self.kcl_residual_A = 0.0
        self.thermal_residual_W = 0.0
        self.counts = {k: 0 for k in ("generate", "recombine", "vacancy_hop", "interstitial_hop")}

    def initialize(self, mode="random", fraction=0.04):
        """Equal initial defect counts. Seeded filament is a specified LRS, not forming evidence."""
        if self.event_count or self.time_s:
            raise ValueError("initialize only a fresh simulator")
        if not 0 <= fraction <= 1:
            raise ValueError("fraction must lie in [0,1]")
        self.vacancy[:] = False; self.interstitial[:] = False
        if mode == "pristine":
            return
        if mode == "random":
            count = int(round(self.n * fraction))
            self.vacancy[self.rng.choice(self.n, count, replace=False)] = True
        elif mode == "filament":
            ids = np.arange(self.n).reshape(self.cfg.shape)
            self.vacancy[ids[self.cfg.shape[0] // 2, self.cfg.shape[1] // 2, :]] = True
        else:
            raise ValueError("mode must be pristine, random, or filament")
        count = int(self.vacancy.sum())
        self.interstitial[self.rng.choice(self.n, count, replace=False)] = True
        self.check_invariants()

    def _network(self, conductance, bottom_g, top_g):
        i, j = self.edges.T
        diagonal = np.bincount(np.r_[i, j], weights=np.r_[conductance, conductance], minlength=self.n)
        diagonal[self.bottom] += bottom_g
        diagonal[self.top] += top_g
        rows = np.r_[i, j, np.arange(self.n)]
        cols = np.r_[j, i, np.arange(self.n)]
        data = np.r_[-conductance, -conductance, diagonal]
        return coo_matrix((data, (rows, cols)), shape=(self.n, self.n)).tocsr()

    def solve_fields(self, applied_voltage_V):
        """Solve div(sigma grad phi)=0, ideal compliance, and steady heat balance.

        Electrodes at z faces (half-cell contacts), insulating lateral walls.
        Conductivity is independent of T; only event rates are thermally activated.
        """
        if not math.isfinite(applied_voltage_V):
            raise ValueError("voltage must be finite")
        c = self.cfg
        sigma = np.where(self.vacancy, c.sigma_on_S_m, c.sigma_off_S_m)
        i, j = self.edges.T
        g = 2 * sigma[i] * sigma[j] / (sigma[i] + sigma[j]) * c.spacing_m
        gb = 2 * sigma[self.bottom] * c.spacing_m
        gt = 2 * sigma[self.top] * c.spacing_m
        matrix = self._network(g, gb, gt)
        rhs = np.zeros(self.n); rhs[self.top] = gt
        unit = spsolve(matrix, rhs)
        geq = float(np.sum(gt * (1 - unit[self.top])))
        if not math.isfinite(geq) or geq <= 0:
            raise RuntimeError("nonpositive terminal conductance")
        voltage = math.copysign(min(abs(applied_voltage_V), c.compliance_A / geq), applied_voltage_V)
        self.phi = unit * voltage
        self.device_voltage_V = voltage
        self.conductance_S = geq
        self.current_A = geq * voltage
        p_edge = g * (self.phi[i] - self.phi[j]) ** 2
        p_bottom = gb * self.phi[self.bottom] ** 2
        p_top = gt * (voltage - self.phi[self.top]) ** 2
        self.power_W = float(p_edge.sum() + p_bottom.sum() + p_top.sum())
        self.energy_residual_W = self.power_W - self.current_A * voltage
        self.kcl_residual_A = float(np.max(np.abs(matrix @ self.phi - rhs * voltage)))
        # Half of each contact's dissipated heat is deposited in the electrode.
        p = np.bincount(np.r_[i, j], weights=np.r_[p_edge / 2, p_edge / 2], minlength=self.n)
        p[self.bottom] += p_bottom / 2; p[self.top] += p_top / 2
        kg = c.thermal_W_mK * c.spacing_m
        heat = self._network(np.full(len(i), kg), np.full(len(self.bottom), 2*kg),
                             np.full(len(self.top), 2*kg))
        rise = spsolve(heat, p)
        self.temperature = c.ambient_K + rise
        outflow = 2*kg * (rise[self.bottom].sum() + rise[self.top].sum())
        self.thermal_residual_W = float(outflow - p.sum())
        if not np.all(np.isfinite(self.temperature)) or self.temperature.max() > c.max_temperature_K:
            raise RuntimeError("temperature exceeds model validity limit; reduce bias or compliance")
        return {"conductance_S": geq, "voltage_V": voltage, "current_A": self.current_A,
                "power_W": self.power_W, "max_temperature_K": float(self.temperature.max())}

    def _rate(self, delta, barrier, temp):
        # Metropolis transition-state construction: exact forward/reverse ratio
        # at fixed phi and uniform T; barrier never becomes negative.
        return self.cfg.attempt_Hz * math.exp(-(barrier + max(0.0, delta)) / (KB * temp))

    def events(self):
        """Enumerate every allowed event once; boundaries cannot wrap or lose oxygen."""
        c = self.cfg; v = self.vacancy; o = self.interstitial
        occupied_neighbors = np.array([sum(v[j] for j in ns) for ns in self.neighbors])
        events = []
        def add(kind, i, j, delta, barrier):
            temp = (self.temperature[i] + self.temperature[j]) / 2
            rate = self._rate(float(delta), barrier, temp)
            if rate > 0:
                events.append(Event(kind, i, j, float(delta), rate))
        for i in range(self.n):
            # Distinct oxygen/interstitial sublattices can both occupy a cell.
            # Pair creation and recombination use exactly matching site pairs.
            for j in [i] + self.neighbors[i]:
                delta = c.pair_energy_eV + c.vacancy_charge_e*(self.phi[i]-self.phi[j]) - c.vacancy_binding_eV*occupied_neighbors[i]
                if not v[i] and not o[j]:
                    add("generate", i, j, delta, c.pair_barrier_eV)
                elif v[i] and o[j]:
                    add("recombine", i, j, -delta, c.pair_barrier_eV)
            for j in self.neighbors[i]:
                if v[i] and not v[j]:
                    nb_new = occupied_neighbors[j] - 1  # moving vacancy no longer at i
                    delta = c.vacancy_charge_e*(self.phi[j]-self.phi[i]) - c.vacancy_binding_eV*(nb_new-occupied_neighbors[i])
                    add("vacancy_hop", i, j, delta, c.vacancy_barrier_eV)
                if o[i] and not o[j]:
                    delta = -c.vacancy_charge_e*(self.phi[j]-self.phi[i])
                    add("interstitial_hop", i, j, delta, c.interstitial_barrier_eV)
        return events

    def apply(self, event):
        i, j = event.source, event.target
        if event.kind == "generate":
            if self.vacancy[i] or self.interstitial[j]: raise ValueError("invalid generation")
            self.vacancy[i] = True; self.interstitial[j] = True
        elif event.kind == "recombine":
            if not self.vacancy[i] or not self.interstitial[j]: raise ValueError("invalid recombination")
            self.vacancy[i] = False; self.interstitial[j] = False
        elif event.kind in ("vacancy_hop", "interstitial_hop"):
            arr = self.vacancy if event.kind == "vacancy_hop" else self.interstitial
            if not arr[i] or arr[j] or j not in self.neighbors[i]: raise ValueError("invalid hop")
            arr[i] = False; arr[j] = True
        else:
            raise ValueError("unknown event kind")
        self.event_count += 1; self.counts[event.kind] += 1
        self.check_invariants()

    def check_invariants(self):
        if int(self.vacancy.sum()) != int(self.interstitial.sum()):
            raise RuntimeError("closed-oxygen conservation violated: N_v != N_i")

    def step(self, applied_voltage_V, until_s=math.inf):
        """Exact direct SSA step for a constant-bias interval.

        At a scheduled voltage boundary, advance to the boundary WITHOUT firing
        the sampled future event. A new exponential draw is made at new bias.
        """
        if until_s < self.time_s: raise ValueError("until_s is in the past")
        self.solve_fields(applied_voltage_V)
        catalog = self.events()
        rates = np.array([e.rate_Hz for e in catalog])
        total = float(rates.sum())
        if not math.isfinite(total): raise RuntimeError("nonfinite total rate")
        if total <= 0:
            if math.isfinite(until_s): self.time_s = until_s
            return None
        dt = -math.log(max(self.rng.random(), np.finfo(float).tiny)) / total
        if self.time_s + dt >= until_s:
            self.time_s = until_s
            return None
        index = min(int(np.searchsorted(np.cumsum(rates), self.rng.random()*total, side="right")), len(catalog)-1)
        self.time_s += dt
        chosen = catalog[index]
        self.apply(chosen)
        return chosen

    def spanning(self):
        clusters, _ = label(self.vacancy.reshape(self.cfg.shape))  # six-connected
        a = set(clusters[:, :, 0].ravel()) - {0}
        b = set(clusters[:, :, -1].ravel()) - {0}
        return bool(a & b)

    def snapshot(self, voltage):
        self.solve_fields(voltage)  # diagnostics MUST correspond to post-event state
        return {"time_s": self.time_s, "events": self.event_count, "applied_voltage_V": voltage,
                "device_voltage_V": self.device_voltage_V, "current_A": self.current_A,
                "conductance_S": self.conductance_S, "resistance_ohm": 1/self.conductance_S,
                "frozen_read_0p1V_current_A": self.conductance_S * 0.1,
                "power_W": self.power_W, "max_temperature_K": float(self.temperature.max()),
                "vacancies": int(self.vacancy.sum()), "interstitials": int(self.interstitial.sum()),
                "spanning": self.spanning(), "kcl_residual_A": self.kcl_residual_A,
                "energy_residual_W": self.energy_residual_W, "thermal_residual_W": self.thermal_residual_W}

    def config_dict(self):
        return asdict(self.cfg)
