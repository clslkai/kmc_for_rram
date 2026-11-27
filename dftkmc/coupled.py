"""Variable-vacancy KMC with finite oxygen stores and effective transport.

At each state a positive resistor network is solved. Its node voltages are used
as a phenomenological ionic electrostatic proxy and FROZEN while enumerating
and selecting one event. The network is refreshed afterwards. Shared saddles
give local detailed balance at one frozen potential; this is NOT a claim of
global detailed balance or a state-function electrostatic energy at finite bias.
The original fixed-N bulk/slab Simulator remains a separate reference engine.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from dataclasses import asdict, dataclass, field
import json
import math
from numbers import Integral
from pathlib import Path

import numpy as np

from .device import DeviceConfig, ElectrodeConfig, _exact_fields, _finite_real, _json_constant, _json_pairs
from .exchange import ExchangeConfig, FiniteReservoir, exchange_pair
from .model import KB, ModelDomainError, SimulationConfig, Simulator
from .transport import ElectronicNetwork, TransportConfig


FIELD_CLOSURE = "quasistatic_kirchhoff_frozen_during_event"


def _integer(value, name, minimum=0):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


@dataclass(frozen=True)
class CoupledConfig:
    repeats: tuple = (2, 2, 4)
    charge: int = 2
    temperature_K: float = 800.
    seed: int = 17
    bottom_potential_V: float = 0.
    top_potential_V: float = .5
    bottom_contact_depth_A: float = 2.
    top_contact_depth_A: float = 2.
    pair_scale: float = 1.
    environment_scale: float = 1.
    exchange: ExchangeConfig = field(default_factory=ExchangeConfig)
    transport: TransportConfig = field(default_factory=TransportConfig)
    formation_min_conductance_S: float = 1e-9
    formation_min_gain: float = 1000.
    kinetic_field_closure: str = FIELD_CLOSURE
    schema_version: int = 1

    def __post_init__(self):
        if len(self.repeats) != 3:
            raise ValueError("repeats requires three integers")
        object.__setattr__(self, "repeats", tuple(_integer(x, "repeats", 2) for x in self.repeats))
        if isinstance(self.charge, (bool, np.bool_)) or not isinstance(self.charge, Integral):
            raise ValueError("charge must be an integer fixed charge branch")
        object.__setattr__(self, "charge", int(self.charge))
        object.__setattr__(self, "seed", _integer(self.seed, "seed"))
        for name in ("temperature_K", "bottom_contact_depth_A", "top_contact_depth_A",
                     "formation_min_conductance_S", "formation_min_gain"):
            object.__setattr__(self, name, _finite_real(getattr(self, name), name, positive=True))
        if self.formation_min_gain <= 1:
            raise ValueError("formation_min_gain must exceed one")
        for name in ("bottom_potential_V", "top_potential_V", "pair_scale", "environment_scale"):
            object.__setattr__(self, name, _finite_real(getattr(self, name), name))
        if not math.isfinite(self.top_potential_V-self.bottom_potential_V):
            raise ValueError("the applied voltage difference must be finite")
        if min(self.pair_scale, self.environment_scale) < 0:
            raise ValueError("ablation scales must be nonnegative")
        if not isinstance(self.exchange, ExchangeConfig) or not isinstance(self.transport, TransportConfig):
            raise ValueError("exchange and transport must be their validated configuration objects")
        if self.kinetic_field_closure != FIELD_CLOSURE:
            raise ValueError(f"only {FIELD_CLOSURE} is implemented")
        if _integer(self.schema_version, "schema_version", 1) != 1:
            raise ValueError("coupled schema_version must be 1")

    @classmethod
    def from_dict(cls, raw):
        _exact_fields(raw, cls.__dataclass_fields__, "coupled configuration")
        values = dict(raw)
        values["exchange"] = ExchangeConfig.from_dict(values["exchange"])
        values["transport"] = TransportConfig.from_dict(values["transport"])
        return cls(**values)

    @classmethod
    def load(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8"),
                                        object_pairs_hook=_json_pairs, parse_constant=_json_constant))

    def to_dict(self):
        result = asdict(self)
        result["repeats"] = list(self.repeats)
        return result

    def write(self, path):
        destination = Path(path)
        if destination.exists():
            raise FileExistsError(f"Configuration already exists: {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(self.to_dict(), indent=2, allow_nan=False)+"\n", encoding="utf-8")


@dataclass(frozen=True)
class CoupledEvent:
    kind: str
    source: int
    target: int
    reservoir: str | None
    rate_Hz: float
    barrier_eV: float
    chemical_delta_eV: float
    oxide_delta_eV: float
    reservoir_delta_eV: float
    electrical_work_eV: float
    delta_eV: float
    channel_id: str
    displacement_A: tuple
    state_epoch: int
    q: int
    direction: int
    source_potential_V: float
    target_potential_V: float
    electrode_potential_V: float | None
    attempt_Hz: float
    parameter_status: str


class CoupledDeviceSimulator:
    def __init__(self, catalog, config=CoupledConfig()):
        if not isinstance(config, CoupledConfig):
            raise ValueError("config must be CoupledConfig")
        self.catalog, self.cfg, self.config = catalog, config, config
        # Reuse immutable geometry and the declared ionic/pair graph, not the
        # legacy state engine. DeviceBoundary is used only for contact geometry;
        # its linear field and blocking-oxygen metadata do not govern this run.
        geometry = DeviceConfig(
            bottom=ElectrodeConfig(config.bottom_potential_V, config.bottom_contact_depth_A),
            top=ElectrodeConfig(config.top_potential_V, config.top_contact_depth_A),
        )
        base = Simulator(catalog, SimulationConfig(
            repeats=config.repeats, charge=config.charge, temperature_K=config.temperature_K,
            seed=config.seed, boundary="slab", device=geometry,
            pair_scale=config.pair_scale, environment_scale=config.environment_scale))
        self.lattice, self.energy, self.edges = base.lattice, base.energy, base.edges
        self.geometry_boundary = base.device
        self.charge_model = base.charge_model
        self.n = self.lattice.n
        if not 0 <= config.exchange.reference_basis < self.lattice.basis_size:
            raise ValueError("exchange.reference_basis is outside the catalogue basis")
        reference_energy = float(catalog.site_energies(config.charge)[config.exchange.reference_basis])
        self.formation_site_eV = (self.energy.site-reference_energy+
                                  config.exchange.vacancy_formation_offset_eV)
        self.network = ElectronicNetwork(self.lattice, self.geometry_boundary, config.transport)
        self.reservoirs = {name: FiniteReservoir(getattr(config.exchange, name)) for name in ("bottom", "top")}
        self.reservoir_oxygen_ids = {name: [] for name in self.reservoirs}
        self.electrons_to_electrode = dict.fromkeys(self.reservoirs, 0)
        self.rng = np.random.default_rng(config.seed)
        self.vacancy = np.zeros(self.n, dtype=bool)
        self.oxygen_ids = np.arange(self.n, dtype=np.int64)
        self.vacancy_ids = np.full(self.n, -1, dtype=np.int64)
        self._next_vacancy_id = 0
        self._retired_vacancy_ids = set()
        self._initial_oxygen_ids = set()
        self._initialized = False
        self.state_epoch = 0
        self.event_count, self.time_s = 0, 0.
        self.counts = Counter()
        self.requested_voltage_V = config.top_potential_V-config.bottom_potential_V
        self.fields = None
        self.first_spanning_time_s = None
        self.first_formation_time_s = None
        self.formation_observations = []
        self.peak_conductance_S = 0.
        self._last_criterion = False
        self._criterion_start_s = None
        self._criterion_durations = []
        self._pulse_index = 0
        self._edge_lookup = {(e.source, e.target, e.channel_id) for e in self.edges}

    def initialize(self, n_vacancies=0, vacancy_indices=None):
        if self._initialized:
            raise ValueError("simulation already initialized")
        if vacancy_indices is None:
            count = _integer(n_vacancies, "n_vacancies")
            if count > self.n:
                raise ValueError("n_vacancies exceeds the number of oxygen sites")
            indices = np.sort(self.rng.choice(self.n, count, replace=False))
        else:
            raw = np.asarray(vacancy_indices)
            if raw.ndim != 1 or (raw.size and not np.issubdtype(raw.dtype, np.integer)):
                raise ValueError("vacancy_indices must be a one-dimensional integer array")
            indices = raw.astype(int)
            if len(set(indices)) != len(indices) or (indices < 0).any() or (indices >= self.n).any():
                raise ValueError("vacancy_indices must be distinct valid sites")
        self.vacancy[indices] = True
        self.oxygen_ids[indices] = -1
        self.vacancy_ids[indices] = np.arange(len(indices))
        self._next_vacancy_id = len(indices)
        self.initial_vacancies = len(indices)
        self.initial_vacancy_indices = indices.tolist()
        next_oxygen = self.n
        for name, reservoir in self.reservoirs.items():
            self.reservoir_oxygen_ids[name] = list(range(next_oxygen, next_oxygen+reservoir.count))
            next_oxygen += reservoir.count
        all_ids = self.oxygen_ids[self.oxygen_ids >= 0].tolist()
        for pool in self.reservoir_oxygen_ids.values():
            all_ids.extend(pool)
        self._initial_oxygen_ids = set(all_ids)
        self.initial_total_oxygen = len(all_ids)
        self._initialized = True
        self._refresh_fields()
        self.initial_conductance_S = self.fields.equivalent_conductance_S
        self.initial_spanning = self.fields.vacancy_graph_spanning
        if self.initial_spanning:
            self.first_spanning_time_s = self.time_s
        self._observe()
        self.check()

    def _refresh_fields(self):
        self.fields = self.network.solve(
            self.vacancy, bottom_V=self.cfg.bottom_potential_V,
            top_V=self.cfg.bottom_potential_V+self.requested_voltage_V)

    def set_bias(self, applied_voltage_V):
        if not self._initialized:
            raise ValueError("initialize first")
        requested = _finite_real(applied_voltage_V, "applied_voltage_V")
        fields = self.network.solve(self.vacancy, bottom_V=self.cfg.bottom_potential_V,
                                    top_V=self.cfg.bottom_potential_V+requested)
        self.requested_voltage_V, self.fields = requested, fields
        self.state_epoch += 1
        self._observe()

    def check(self):
        if not self._initialized:
            raise ValueError("initialize first")
        if not np.array_equal(self.vacancy, self.oxygen_ids < 0):
            raise RuntimeError("oxygen identity and vacancy occupancy disagree")
        if not np.array_equal(self.vacancy, self.vacancy_ids >= 0):
            raise RuntimeError("vacancy identity and occupancy disagree")
        active = self.vacancy_ids[self.vacancy].tolist()
        if len(active) != len(set(active)) or self._retired_vacancy_ids.intersection(active):
            raise RuntimeError("vacancy birth/death identity violation")
        oxygen_ids = self.oxygen_ids[~self.vacancy].tolist()
        for name, reservoir in self.reservoirs.items():
            pool = self.reservoir_oxygen_ids[name]
            if len(pool) != reservoir.count or not 0 <= reservoir.count <= reservoir.config.capacity:
                raise RuntimeError("reservoir atom identities disagree with finite oxygen inventory")
            if reservoir.count != reservoir.initial_count+reservoir.out_count-reservoir.in_count:
                raise RuntimeError("per-electrode oxygen transfer ledger failed")
            if self.electrons_to_electrode[name] != self.cfg.charge*(reservoir.out_count-reservoir.in_count):
                raise RuntimeError("per-electrode reaction electron ledger failed")
            oxygen_ids.extend(pool)
        if len(oxygen_ids) != len(set(oxygen_ids)) or set(oxygen_ids) != self._initial_oxygen_ids:
            raise RuntimeError("closed oxide-plus-reservoir oxygen atom conservation failed")
        if self.cfg.charge*(int(self.vacancy.sum())-self.initial_vacancies) != sum(self.electrons_to_electrode.values()):
            raise RuntimeError("reaction charge compensation ledger failed")

    def chemical_free_energy_eV(self):
        pairs = self.energy.pairs
        pair_energy = float(np.sum(self.energy.weights*self.vacancy[pairs[:, 0]]*self.vacancy[pairs[:, 1]])) if len(pairs) else 0.
        return float(self.formation_site_eV@self.vacancy + pair_energy +
                     math.fsum(pool.free_energy_eV(self.cfg.temperature_K) for pool in self.reservoirs.values()))

    def _frozen_potentials(self, potential_V=None, electrode_potentials_V=None):
        if (potential_V is None) != (electrode_potentials_V is None):
            raise ValueError("a frozen-field override requires both site and electrode potentials")
        if potential_V is None:
            return (self.fields.relative_potential_V,
                    {"bottom": 0., "top": self.fields.applied_voltage_V},
                    self.fields.potential_V,
                    {"bottom": self.fields.bottom_potential_V, "top": self.fields.top_potential_V})
        phi = np.asarray(potential_V, dtype=float)
        if phi.shape != (self.n,) or not np.isfinite(phi).all():
            raise ValueError("frozen site potential must be one finite value per site")
        _exact_fields(electrode_potentials_V, {"bottom", "top"}, "frozen electrode potentials")
        electrodes = {name: _finite_real(value, name) for name, value in electrode_potentials_V.items()}
        reference = electrodes["bottom"]
        return phi-reference, {name: value-reference for name, value in electrodes.items()}, phi, electrodes

    def frozen_field_state_function_eV(self, potential_V, electrode_potentials_V):
        """Oracle at ONE held-fixed potential; not a finite-bias feedback energy."""
        phi, electrodes, _, _ = self._frozen_potentials(potential_V, electrode_potentials_V)
        return self.chemical_free_energy_eV() + self.cfg.charge*(float(phi@self.vacancy)-
               math.fsum(pool.count*electrodes[name] for name, pool in self.reservoirs.items()))

    def events(self, potential_V=None, electrode_potentials_V=None):
        if not self._initialized:
            raise ValueError("initialize first")
        phi, electrode_phi, absolute_phi, absolute_electrodes = self._frozen_potentials(potential_V, electrode_potentials_V)
        q, temp = self.cfg.charge, self.cfg.temperature_K
        mode = self.charge_model["environment_mode"]
        gamma = self.cfg.environment_scale*float(self.charge_model["environment_saddle_eV_per_neighbor"])
        events = []
        for edge in self.edges:
            a, b = edge.source, edge.target
            if self.vacancy[a] == self.vacancy[b]:
                continue
            forward = bool(self.vacancy[a])
            i, j = (a, b) if forward else (b, a)
            direction = 1 if forward else -1
            datum = self.catalog.channel(edge.channel_id, q)
            n3, n4 = self.energy.environment_descriptor(self.vacancy, i, j)
            if mode == "dilute_only" and n3+n4:
                raise ModelDomainError("uncovered neighboring-vacancy environment in coupled hopping")
            if mode not in ("dilute_only", "symmetric_pair"):
                raise ModelDomainError("unknown environment closure")
            chemical = self.energy.delta_hop(self.vacancy, i, j)
            pair_delta = self.energy.environment_pair_delta(self.vacancy, i, j)
            electrical = q*float(phi[j]-phi[i])
            fraction = float(datum["ts_fraction"]) if forward else 1.-float(datum["ts_fraction"])
            barrier = (float(datum["saddle_eV"])-self.energy.site[i] + .5*pair_delta+
                       gamma*(n3+n4)+fraction*electrical)
            if not math.isfinite(barrier) or barrier < -1e-10:
                raise ModelDomainError("negative/nonfinite coupled hop barrier; voltage/environment exceeds catalogue domain")
            barrier = max(barrier, 0.)  # existing floating-point tolerance only
            attempt = float(datum["attempt_Hz"])
            rate = attempt*math.exp(-barrier/(KB*temp))
            if rate:
                events.append(CoupledEvent("hop", i, j, None, rate, barrier, chemical, chemical, 0.,
                    electrical, chemical+electrical, edge.channel_id,
                    tuple(np.asarray(edge.displacement_A)*direction), self.state_epoch, q, direction,
                    float(absolute_phi[i]), float(absolute_phi[j]), None, attempt, str(datum["status"])))
        for name, reservoir in self.reservoirs.items():
            if not reservoir.config.enabled:
                continue
            mask = getattr(self.geometry_boundary, name+"_mask")
            for site in np.flatnonzero(mask):
                site = int(site)
                out = not self.vacancy[site]
                direction = 1 if out else -1
                if not reservoir.available(direction):
                    continue
                oxide_cost = float(self.formation_site_eV[site] +
                    math.fsum(weight for neighbor, weight in self.energy.adjacency[site] if self.vacancy[neighbor]))
                pair = exchange_pair(oxide_cost, reservoir.config, reservoir.count if out else reservoir.count-1,
                                     q, float(phi[site]), electrode_phi[name], temp)
                prefix = "out" if out else "in"
                rate = pair[prefix+"_rate_Hz"]
                if rate:
                    events.append(CoupledEvent("oxygen_"+prefix, site, -1, name, rate,
                        pair[prefix+"_barrier_eV"], direction*pair["chemical_delta_eV"],
                        direction*pair["oxide_delta_eV"], direction*pair["reservoir_delta_eV"],
                        direction*pair["electrical_work_eV"], direction*pair["delta_eV"],
                        f"interface:{name}:basis{int(self.lattice.basis_id[site])}", (0., 0., 0.),
                        self.state_epoch, q, direction, float(absolute_phi[site]),
                        float(absolute_phi[site]), float(absolute_electrodes[name]), reservoir.config.attempt_Hz,
                        self.cfg.exchange.calibration_status))
        if any(not math.isfinite(e.rate_Hz) or e.rate_Hz <= 0 for e in events):
            raise ModelDomainError("invalid coupled event rate")
        return events

    def apply(self, event):
        if not self._initialized:
            raise ValueError("initialize first")
        if not isinstance(event, CoupledEvent) or event.state_epoch != self.state_epoch or event.q != self.cfg.charge:
            raise ValueError("illegal or stale coupled event")
        i, j = event.source, event.target
        if not 0 <= i < self.n:
            raise ValueError("invalid event source")
        self.check()
        candidate = self.vacancy.copy()
        if event.kind == "hop":
            if not 0 <= j < self.n or not self.vacancy[i] or self.vacancy[j]:
                raise ValueError("illegal vacancy hop")
            if ((i, j, event.channel_id) not in self._edge_lookup and
                    (j, i, event.channel_id) not in self._edge_lookup):
                raise ValueError("undeclared ionic hop")
            candidate[i], candidate[j] = False, True
        elif event.kind in ("oxygen_out", "oxygen_in"):
            name = event.reservoir
            if name not in self.reservoirs or not getattr(self.geometry_boundary, name+"_mask")[i]:
                raise ValueError("exchange must belong to its electrode contact region")
            out = event.kind == "oxygen_out"
            direction = 1 if out else -1
            reservoir = self.reservoirs[name]
            if out == bool(self.vacancy[i]) or not reservoir.available(direction):
                raise ValueError("exchange violates vacancy state or finite reservoir inventory")
            candidate[i] = out
        else:
            raise ValueError("unknown coupled event kind")
        # Solve before committing atom identities or ledgers. If the candidate
        # network is outside its numerical domain, all saved state still belongs
        # to the last accepted event, including its matching electronic fields.
        next_fields = self.network.solve(candidate, bottom_V=self.cfg.bottom_potential_V,
                                        top_V=self.cfg.bottom_potential_V+self.requested_voltage_V)
        before_vacancies = int(self.vacancy.sum())
        before_pools = {name: pool.count for name, pool in self.reservoirs.items()}
        if event.kind == "hop":
            oxygen_id, vacancy_id = int(self.oxygen_ids[j]), int(self.vacancy_ids[i])
            self.oxygen_ids[i], self.oxygen_ids[j] = oxygen_id, -1
            self.vacancy_ids[i], self.vacancy_ids[j] = -1, vacancy_id
            self.vacancy[i], self.vacancy[j] = False, True
            oxygen_from, oxygen_to = f"oxide:{j}", f"oxide:{i}"
            identity_action = "moved"
        else:
            if out:
                oxygen_id = int(self.oxygen_ids[i])
                vacancy_id = self._next_vacancy_id
                self._next_vacancy_id += 1
                self.reservoir_oxygen_ids[name].append(oxygen_id)
                self.oxygen_ids[i], self.vacancy_ids[i], self.vacancy[i] = -1, vacancy_id, True
                oxygen_from, oxygen_to, identity_action = f"oxide:{i}", f"reservoir:{name}", "born"
            else:
                # Atom labels are bookkeeping; reservoir microscopic trajectories
                # are unresolved. A deterministic last-in label convention avoids
                # spending physical KMC random numbers on an unmodelled label.
                oxygen_id = self.reservoir_oxygen_ids[name].pop()
                vacancy_id = int(self.vacancy_ids[i])
                self._retired_vacancy_ids.add(vacancy_id)
                self.oxygen_ids[i], self.vacancy_ids[i], self.vacancy[i] = oxygen_id, -1, False
                oxygen_from, oxygen_to, identity_action = f"reservoir:{name}", f"oxide:{i}", "removed"
            reservoir.apply(direction)
            self.electrons_to_electrode[name] += self.cfg.charge*direction
        self.event_count += 1
        self.state_epoch += 1
        self.counts[event.kind] += 1
        self.check()
        self.fields = next_fields
        self._observe()
        return dict(oxygen_id=oxygen_id, oxygen_from=oxygen_from, oxygen_to=oxygen_to,
                    vacancy_id=vacancy_id, vacancy_identity_action=identity_action,
                    vacancies_before=before_vacancies, vacancies_after=int(self.vacancy.sum()),
                    reservoir_oxygen_before=before_pools,
                    reservoir_oxygen_after={name: pool.count for name, pool in self.reservoirs.items()})

    def _gain(self):
        if self.initial_conductance_S <= 0:
            return None
        return self.fields.equivalent_conductance_S/self.initial_conductance_S

    def _criterion(self):
        gain = self._gain()
        return bool(self.fields.vacancy_graph_spanning and gain is not None and
                    gain >= self.cfg.formation_min_gain and
                    self.fields.equivalent_conductance_S >= self.cfg.formation_min_conductance_S)

    def _observe(self):
        self.peak_conductance_S = max(self.peak_conductance_S, self.fields.equivalent_conductance_S)
        if self.fields.vacancy_graph_spanning and self.first_spanning_time_s is None:
            self.first_spanning_time_s = self.time_s
        criterion = self._criterion()
        if criterion != self._last_criterion:
            self.formation_observations.append(dict(time_s=self.time_s, event_count=self.event_count,
                criterion_met=criterion, conductance_S=self.fields.equivalent_conductance_S,
                vacancy_graph_spanning=self.fields.vacancy_graph_spanning))
            if criterion:
                self._criterion_start_s = self.time_s
                if self.first_formation_time_s is None:
                    self.first_formation_time_s = self.time_s
            else:
                self._criterion_durations.append(self.time_s-self._criterion_start_s)
                self._criterion_start_s = None
        self._last_criterion = criterion

    def snapshot(self):
        if not self._initialized:
            raise ValueError("initialize first")
        counts = {name: pool.count for name, pool in self.reservoirs.items()}
        nv = int(self.vacancy.sum())
        total = self.n-nv+sum(counts.values())
        bottom = int(np.count_nonzero(self.vacancy & self.geometry_boundary.bottom_mask))
        top = int(np.count_nonzero(self.vacancy & self.geometry_boundary.top_mask))
        durations = self._criterion_durations + ([self.time_s-self._criterion_start_s] if self._criterion_start_s is not None else [])
        return dict(time_s=self.time_s, events=self.event_count, pulse_index=self._pulse_index,
            vacancies=nv, oxide_oxygen=self.n-nv, bottom_reservoir_oxygen=counts["bottom"],
            top_reservoir_oxygen=counts["top"], total_oxygen=total,
            mass_residual=total-self.initial_total_oxygen, oxygen_balance_residual=total-self.initial_total_oxygen,
            reaction_charge_residual_e=self.cfg.charge*(nv-self.initial_vacancies)-sum(self.electrons_to_electrode.values()),
            reaction_electrons_to_bottom=self.electrons_to_electrode["bottom"],
            reaction_electrons_to_top=self.electrons_to_electrode["top"],
            bottom_contact_vacancies=bottom, top_contact_vacancies=top, interior_vacancies=nv-bottom-top,
            chemical_free_energy_eV=self.chemical_free_energy_eV(),
            conductance_gain=self._gain(), model_formation_criterion=self._criterion(),
            first_spanning_time_s=self.first_spanning_time_s, first_formation_time_s=self.first_formation_time_s,
            longest_criterion_duration_s=max(durations, default=0.), peak_conductance_S=self.peak_conductance_S,
            **self.fields.summary())

    def step(self, until_s=math.inf):
        if math.isnan(until_s) or until_s < self.time_s:
            raise ValueError("invalid time boundary")
        choices = self.events()
        total = math.fsum(e.rate_Hz for e in choices)
        if not math.isfinite(total):
            raise ModelDomainError("nonfinite total event rate")
        if total == 0:
            if not math.isfinite(until_s):
                raise ModelDomainError("no active events and no finite stopping time")
            self.time_s = float(until_s)
            return None
        rng_before = deepcopy(self.rng.bit_generator.state)
        dt = -math.log(max(float(self.rng.random()), np.finfo(float).tiny))/total
        if self.time_s+dt >= until_s:
            self.time_s = float(until_s)
            return None
        if self.time_s+dt == self.time_s:
            raise ModelDomainError("waiting time is below floating-point resolution")
        draw = float(self.rng.random())*total
        chosen = choices[-1]
        for event in choices:
            draw -= event.rate_Hz
            if draw < 0:
                chosen = event
                break
        time_before = self.time_s
        before = self.fields
        self.time_s += dt
        try:
            ledger = self.apply(chosen)
        except (ModelDomainError, RuntimeError, ValueError):
            self.time_s = time_before
            self.rng.bit_generator.state = rng_before
            raise
        return dict(time_s=self.time_s, time_before_s=time_before, dt_s=dt,
                    event_index=self.event_count, pulse_index=self._pulse_index,
                    requested_voltage_V=self.requested_voltage_V,
                    applied_voltage_before_V=before.applied_voltage_V,
                    applied_voltage_after_V=self.fields.applied_voltage_V,
                    conductance_before_S=before.equivalent_conductance_S,
                    conductance_after_S=self.fields.equivalent_conductance_S,
                    current_after_A=self.fields.terminal_current_A,
                    potential_change_max_V=float(np.max(np.abs(self.fields.relative_potential_V-before.relative_potential_V))),
                    vacancy_graph_spanning_after=self.fields.vacancy_graph_spanning,
                    model_formation_criterion_after=self._criterion(),
                    mass_residual=self.n-int(self.vacancy.sum())+
                        sum(pool.count for pool in self.reservoirs.values())-self.initial_total_oxygen,
                    reaction_charge_residual_e=self.cfg.charge*(int(self.vacancy.sum())-self.initial_vacancies)-
                        sum(self.electrons_to_electrode.values()),
                    **asdict(chosen), **ledger)

    def run(self, duration_s, max_events=10000, sample_every=10):
        return self.run_waveform([(self.requested_voltage_V, duration_s)], max_events=max_events, sample_every=sample_every)

    def run_waveform(self, pulses, max_events=10000, sample_every=10):
        if not self._initialized:
            raise ValueError("initialize first")
        max_events = _integer(max_events, "max_events", 1)
        sample_every = _integer(sample_every, "sample_every", 1)
        if not pulses:
            raise ValueError("at least one voltage segment is required")
        segments = []
        for bias, duration in pulses:
            segments.append((_finite_real(bias, "pulse voltage"), _finite_real(duration, "pulse duration", positive=True)))
        start, limit = self.time_s, self.event_count+max_events
        target = start+math.fsum(duration for _, duration in segments)
        if not math.isfinite(target) or target <= start:
            raise ValueError("invalid total waveform duration")
        records, logs, states, oxygen_states, vacancy_id_states, potentials = [], [], [], [], [], []
        edge_currents, edge_fields = [], []

        def capture():
            records.append(self.snapshot())
            states.append(self.vacancy.copy())
            oxygen_states.append(self.oxygen_ids.copy())
            vacancy_id_states.append(self.vacancy_ids.copy())
            potentials.append(self.fields.potential_V.copy())
            edge_currents.append(self.fields.edge_current_A.copy())
            edge_fields.append(self.fields.edge_field_V_nm.copy())

        status, error = "completed_waveform", None
        for index, (bias, duration) in enumerate(segments):
            end = start+math.fsum(value[1] for value in segments[:index+1])
            try:
                self.set_bias(bias)
                self._pulse_index = index
                capture()  # both sides of voltage discontinuities are explicit samples
                while self.time_s < end and self.event_count < limit:
                    old_span = self.fields.vacancy_graph_spanning
                    old_criterion = self._criterion()
                    row = self.step(end)
                    if row is not None:
                        logs.append(row)
                        if (self.event_count % sample_every == 0 or
                                old_span != self.fields.vacancy_graph_spanning or old_criterion != self._criterion()):
                            capture()
            except (ModelDomainError, RuntimeError, ValueError) as exc:
                status, error = "model_domain_error", str(exc)
            if not records or records[-1]["time_s"] != self.time_s or records[-1]["events"] != self.event_count:
                capture()
            if error is not None:
                break
            if self.time_s < end:
                status = "event_limit_before_waveform_end"
                break
        return dict(records=records, events=logs, states=states, oxygen_ids=oxygen_states,
                    vacancy_ids=vacancy_id_states, potentials_V=potentials,
                    edge_current_A=edge_currents, edge_field_V_nm=edge_fields,
                    status=status, error=error, requested_end_time_s=target,
                    target_reached=self.time_s == target and error is None,
                    pulses=[dict(requested_voltage_V=bias, duration_s=duration) for bias, duration in segments])
