"""Fixed-charge lattice KMC with traceable reversible channels.

The external uniform field is prescribed directly or by two ideal electrode
potentials, NOT obtained from a resistor graph or a space-charge calculation.
Environment saddle shifts are declared approximations. No charge capture,
space-charge solution, oxygen exchange or electronic conductance is implied.
"""
from collections import Counter
from dataclasses import asdict, dataclass
import math
from numbers import Integral

import numpy as np

from .geometry import Lattice
from .device import DeviceBoundary, DeviceConfig

KB = 8.617333262145e-5  # eV/K


class ModelDomainError(ValueError):
    """Parameters/state exceed the declared activated-event model."""


@dataclass(frozen=True)
class SimulationConfig:
    repeats: tuple = (3, 3, 3)
    charge: int = 2
    temperature_K: float = 600.
    seed: int = 17
    boundary: str = "bulk"
    field_V_nm: tuple = (0., 0., 0.)
    pair_scale: float = 1.
    environment_scale: float = 1.
    device: DeviceConfig | None = None

    def __post_init__(self):
        if self.boundary not in ("bulk", "slab"):
            raise ValueError("boundary must be bulk or slab")
        if isinstance(self.charge, (bool, np.bool_)) or not isinstance(self.charge, Integral):
            raise ValueError("charge must be an integer electron-count charge state")
        if isinstance(self.seed, (bool, np.bool_)) or not isinstance(self.seed, Integral) or self.seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        if not math.isfinite(self.temperature_K) or self.temperature_K <= 0:
            raise ValueError("temperature_K must be finite and positive")
        if len(self.field_V_nm) != 3 or not np.isfinite(self.field_V_nm).all():
            raise ValueError("field_V_nm must be a finite length-three vector")
        if any(not math.isfinite(x) or x < 0 for x in (self.pair_scale, self.environment_scale)):
            raise ValueError("ablation scales must be finite and nonnegative")
        if self.device is not None:
            if not isinstance(self.device, DeviceConfig):
                raise ValueError("device must be a DeviceConfig; use DeviceConfig.from_dict for JSON")
            if self.boundary != "slab":
                raise ValueError("an electrode device requires boundary='slab'")
            if np.any(np.asarray(self.field_V_nm) != 0.):
                raise ValueError("device electrode voltages and a separate uniform field cannot be combined")


@dataclass(frozen=True)
class Event:
    source: int
    target: int
    channel_id: str
    direction: int
    q: int
    rate_Hz: float
    barrier_eV: float
    delta_eV: float
    env_n3: int
    env_n4: int
    displacement_A: tuple
    chemical_delta_eV: float
    field_work_eV: float
    parameter_status: str


class EnergyModel:
    def __init__(self, lattice, catalog, q, pair_scale=1.):
        self.lattice = lattice
        self.site = np.asarray(catalog.site_energies(q), dtype=float)[lattice.basis_id]
        model = catalog.charge_model(q)
        pair_data = model["pair_energies_eV"]
        self.pairs = lattice.interaction_pairs
        self.weights = np.asarray([pair_data["-".join(sorted((lattice.kind[i], lattice.kind[j])))]*pair_scale
                                   for i, j in self.pairs], dtype=float)
        if not np.isfinite(self.site).all() or not np.isfinite(self.weights).all():
            raise ValueError("unresolved site or pair energies")
        self.adjacency = [[] for _ in range(lattice.n)]
        for (i, j), weight in zip(self.pairs, self.weights):
            self.adjacency[i].append((int(j), float(weight)))
            self.adjacency[j].append((int(i), float(weight)))

    def total(self, vacancy):
        v = np.asarray(vacancy, dtype=bool)
        pair_energy = np.sum(self.weights*v[self.pairs[:, 0]]*v[self.pairs[:, 1]]) if len(self.pairs) else 0.
        return float(self.site@v + pair_energy)

    def environment_pair_delta(self, vacancy, i, j):
        # i--j contributes zero both before and after a legal single-vacancy hop.
        removed = sum(weight for k, weight in self.adjacency[i] if k != j and vacancy[k])
        added = sum(weight for k, weight in self.adjacency[j] if k != i and vacancy[k])
        return float(added-removed)

    def delta_hop(self, vacancy, i, j):
        if not vacancy[i] or vacancy[j]:
            raise ValueError("energy requested for an illegal hop")
        return float(self.site[j]-self.site[i]+self.environment_pair_delta(vacancy, i, j))

    def environment_descriptor(self, vacancy, i, j):
        neighbors = (set(self.lattice.interaction_neighbors[i]) |
                     set(self.lattice.interaction_neighbors[j])) - {i, j}
        # Count distinct external vacancies, not the jumping vacancy or twice-shared neighbors.
        kinds = Counter(self.lattice.kind[k] for k in neighbors if vacancy[k])
        return int(kinds["O3"]), int(kinds["O4"])


class Simulator:
    def __init__(self, catalog, config=SimulationConfig()):
        # Direct Python callers must not bypass completeness/geometry checks.
        # CLI loading separately requires explicit opt-in for illustrative data.
        catalog.validate(allow_illustrative=True)
        self.catalog, self.config, self.cfg = catalog, config, config
        if int(config.charge) not in catalog.charges:
            raise ValueError("charge not covered by catalogue")
        charge_model = catalog.charge_model(config.charge)
        # Building from all declared geometric templates preserves independent graphs.
        templates = catalog.data["channels"]
        if not templates:
            raise ValueError("catalogue has no channels")
        migration_cutoff = max(float(t["length_A"]) for t in templates)+1e-8
        self.lattice = Lattice(config.repeats, config.boundary, migration_cutoff,
                               float(charge_model["pair_cutoff_A"]), catalog.structure)
        if self.lattice.structure_id != catalog.structure_id:
            raise ValueError("catalogue structure identity does not match lattice")
        # Channel coverage was checked by the catalogue; never invent a barrier.
        ids = {t["id"] for t in templates}
        self.edges = [edge for edge in self.lattice.migration_edges if edge.channel_id in ids]
        for edge in self.edges:
            catalog.channel(edge.channel_id, config.charge)
        self.n = self.lattice.n
        self.energy = EnergyModel(self.lattice, catalog, config.charge, config.pair_scale)
        self.charge_model = charge_model
        self.device = DeviceBoundary(self.lattice, config.device) if config.device is not None else None
        self.field_V_A = (self.device.field_V_A.copy() if self.device is not None
                          else np.asarray(config.field_V_nm, dtype=float)/10.)
        self.field_V_nm = 10.*self.field_V_A
        self.rng = np.random.default_rng(int(config.seed))
        self.vacancy = np.zeros(self.n, dtype=bool)
        self.labels = np.full(self.n, -1, dtype=int)
        self.displacements_A = np.empty((0, 3), dtype=float)
        self.initial_positions_A = np.empty((0, 3), dtype=float)
        self.time_s = 0.
        self.event_count = 0
        self.counts = Counter()
        self.channel_counts = Counter()
        self._initialized = False
        self._edge_lookup = {(e.source, e.target, e.channel_id) for e in self.edges}

    def initialize(self, n_vacancies=1, vacancy_indices=None):
        if self._initialized:
            raise ValueError("simulation already initialized")
        if vacancy_indices is None:
            if isinstance(n_vacancies, bool) or int(n_vacancies) != n_vacancies or not 0 < n_vacancies < self.n:
                raise ValueError("n_vacancies must be an integer in [1,N-1]")
            indices = np.sort(self.rng.choice(self.n, int(n_vacancies), replace=False))
        else:
            values = np.asarray(vacancy_indices)
            if values.ndim != 1 or not len(values) or not np.issubdtype(values.dtype, np.integer):
                raise ValueError("vacancy_indices must be a nonempty integer vector")
            indices = values.astype(int)
            if len(set(indices)) != len(indices) or (indices < 0).any() or (indices >= self.n).any():
                raise ValueError("vacancy_indices must be distinct valid sites")
        self.vacancy[indices] = True
        self.labels[indices] = np.arange(len(indices))
        self.displacements_A = np.zeros((len(indices), 3))
        self.initial_positions_A = self.lattice.oxygen_A[indices].copy()
        self._initialized = True
        self.check()

    def check(self):
        ids = self.labels[self.vacancy]
        if int(self.vacancy.sum()) != len(self.displacements_A) or not np.array_equal(np.sort(ids), np.arange(len(ids))):
            raise RuntimeError("vacancy count or identity conservation failed")
        if np.any(self.labels[~self.vacancy] != -1) or not np.isfinite(self.displacements_A).all():
            raise RuntimeError("invalid vacancy labels/displacements")

    def events(self):
        if not self._initialized:
            raise ValueError("initialize first")
        result = []
        q, temp = int(self.cfg.charge), self.cfg.temperature_K
        mode = self.charge_model["environment_mode"]
        gamma = self.cfg.environment_scale*float(self.charge_model["environment_saddle_eV_per_neighbor"])
        for edge in self.edges:
            a, b = edge.source, edge.target
            if self.vacancy[a] == self.vacancy[b]:
                continue
            forward = bool(self.vacancy[a])
            i, j = (a, b) if forward else (b, a)
            direction = 1 if forward else -1
            data = self.catalog.channel(edge.channel_id, q)
            n3, n4 = self.energy.environment_descriptor(self.vacancy, i, j)
            if mode == "dilute_only" and n3+n4:
                raise ModelDomainError(f"uncovered local environment for {edge.channel_id}: n3={n3}, n4={n4}")
            if mode not in ("dilute_only", "symmetric_pair"):
                raise ModelDomainError("unknown environment model")
            chemical_delta = self.energy.delta_hop(self.vacancy, i, j)
            pair_delta = self.energy.environment_pair_delta(self.vacancy, i, j)
            displacement = np.asarray(edge.displacement_A)*direction
            field_work = -q*float(self.field_V_A@displacement)
            fraction = float(data["ts_fraction"]) if forward else 1.-float(data["ts_fraction"])
            barrier = (float(data["saddle_eV"])-self.energy.site[i]
                       + .5*pair_delta + gamma*(n3+n4) + fraction*field_work)
            if not math.isfinite(barrier) or barrier < -1e-10:
                raise ModelDomainError(f"negative/nonfinite barrier for {edge.channel_id}; field/environment exceeds parameter domain")
            barrier = max(barrier, 0.)  # floating-point tolerance only, never a physical clipping model
            rate = float(data["attempt_Hz"])*math.exp(-barrier/(KB*temp))
            if not math.isfinite(rate):
                raise ModelDomainError("nonfinite event rate")
            if rate:
                result.append(Event(i, j, edge.channel_id, direction, q, rate, barrier,
                                    chemical_delta+field_work, n3, n4, tuple(displacement),
                                    chemical_delta, field_work, str(data["status"])))
        return result

    def apply(self, event):
        i, j = event.source, event.target
        if not self.vacancy[i] or self.vacancy[j] or event.q != self.cfg.charge:
            raise ValueError("illegal or stale event")
        if (i, j, event.channel_id) not in self._edge_lookup and (j, i, event.channel_id) not in self._edge_lookup:
            raise ValueError("event is not a declared migration channel")
        label = self.labels[i]
        self.vacancy[i], self.vacancy[j] = False, True
        self.labels[i], self.labels[j] = -1, label
        self.displacements_A[label] += event.displacement_A
        self.event_count += 1
        self.counts["hop"] += 1
        self.channel_counts[f"{event.channel_id}|q={event.q}|direction={event.direction}"] += 1
        self.check()

    def step(self, until_s=math.inf):
        if math.isnan(until_s) or until_s < self.time_s:
            raise ValueError("invalid time boundary")
        events = self.events()
        total = math.fsum(e.rate_Hz for e in events)
        if not math.isfinite(total):
            raise ModelDomainError("nonfinite total event rate")
        if not total:
            if not math.isfinite(until_s):
                raise ModelDomainError("no active events and no finite end time")
            self.time_s = until_s
            return None
        dt = -math.log(max(float(self.rng.random()), np.finfo(float).tiny))/total
        if self.time_s+dt >= until_s:
            self.time_s = float(until_s)
            return None
        if self.time_s+dt == self.time_s:
            raise ModelDomainError("waiting time below floating-point time resolution")
        draw = float(self.rng.random())*total
        chosen = events[-1]
        for event in events:
            draw -= event.rate_Hz
            if draw < 0:
                chosen = event
                break
        self.time_s += dt
        self.apply(chosen)
        return chosen

    def _clusters(self):
        remaining = set(int(i) for i in np.flatnonzero(self.vacancy))
        largest = 0
        while remaining:
            first = remaining.pop()
            stack, size = [first], 1
            while stack:
                node = stack.pop()
                neighbors = set(self.lattice.interaction_neighbors[node]) & remaining
                remaining -= neighbors
                size += len(neighbors)
                stack.extend(neighbors)
            largest = max(largest, size)
        pairs = self.lattice.interaction_pairs
        pair_count = int(np.sum(self.vacancy[pairs[:, 0]] & self.vacancy[pairs[:, 1]])) if len(pairs) else 0
        return pair_count, largest

    def snapshot(self):
        if self.device is not None and not self._initialized:
            raise ValueError("initialize first before requesting a device-state snapshot")
        pair_count, largest = self._clusters()
        d = self.displacements_A
        tensor = d.T@d/len(d) if len(d) else np.zeros((3, 3))
        result = dict(time_s=self.time_s, events=self.event_count, vacancies=int(self.vacancy.sum()),
                      energy_eV=self.energy.total(self.vacancy),
                      mean_displacement_A=d.mean(axis=0).tolist() if len(d) else [0., 0., 0.],
                      msd_A2=float(np.trace(tensor)), msd_tensor_A2=tensor.tolist(),
                      pair_count=pair_count, largest_cluster=largest)
        if self.device is not None:
            device = self.device
            depths = device.depth_A[self.vacancy]
            external_energy = device.electrostatic_energy_eV(self.vacancy, self.cfg.charge)
            bottom = int(np.count_nonzero(self.vacancy & device.bottom_mask))
            top = int(np.count_nonzero(self.vacancy & device.top_mask))
            result.update(
                applied_voltage_V=device.config.top.potential_V-device.config.bottom.potential_V,
                bottom_potential_V=device.config.bottom.potential_V,
                top_potential_V=device.config.top.potential_V,
                bottom_contact_vacancies=bottom, top_contact_vacancies=top,
                interior_vacancies=int(self.vacancy.sum())-bottom-top,
                mean_vacancy_depth_A=float(depths.mean()),
                min_vacancy_depth_A=float(depths.min()), max_vacancy_depth_A=float(depths.max()),
                electrostatic_energy_relative_bottom_eV=external_energy,
                total_energy_relative_bottom_eV=result["energy_eV"]+external_energy,
            )
        return result

    def run(self, duration_s, max_events=10000, sample_every=10):
        if not self._initialized:
            raise ValueError("initialize first")
        if not math.isfinite(duration_s) or duration_s <= 0:
            raise ValueError("duration_s must be finite and positive")
        if any(isinstance(x, bool) or int(x) != x or x < 1 for x in (max_events, sample_every)):
            raise ValueError("event limit and sampling interval must be positive integers")
        end, limit = self.time_s+duration_s, self.event_count+int(max_events)
        records, logs, states, displacement_history = [], [], [], []
        def capture():
            records.append(self.snapshot())
            states.append(self.vacancy.copy())
            displacement_history.append(self.displacements_A.copy())
        capture()
        while self.time_s < end and self.event_count < limit:
            event = self.step(end)
            if event is not None:
                logs.append(dict(time_s=self.time_s, **asdict(event)))
                if self.event_count % sample_every == 0:
                    capture()
        if records[-1]["time_s"] != self.time_s or records[-1]["events"] != self.event_count:
            capture()
        status = "completed_duration" if self.time_s == end else "event_limit_before_end"
        return dict(records=records, events=logs, states=states,
                    displacements=displacement_history, status=status)
