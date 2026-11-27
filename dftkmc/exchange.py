"""Finite oxygen reservoirs and reversible exchange at a *frozen* potential.

Direction +1 means lattice oxygen enters this reservoir and a vacancy is
created.  Direction -1 returns reservoir oxygen to the lattice.  This module
counts oxygen transfers; the coupled engine owns actual oxygen identities and
the reaction-electron ledger (+q electrons to the participating electrode for
creation of a q-positive vacancy).  That ledger is stoichiometric bookkeeping,
not a self-consistent solution of electrode charge or electronic occupation.

The oxide creation energy supplied to ``exchange_pair`` must include a separately
declared formation offset plus ``site_i - site_reference`` and its pair-energy
increment.  ``ExchangeConfig.reference_basis`` fixes this reference.  Migration
site energies alone cannot establish a vacancy formation energy.  A common shift
of every site and saddle energy must therefore leave exchange unchanged.

For C storage positions with m oxygen atoms the declared reservoir free energy is
F(m)=mu0*m+kappa*m*(m-1)/(2*C)-k_B*T*ln(binomial(C,m)); the last term is optional
configurational mixing entropy, NOT a vibrational free-energy correction.  Do not
also multiply rates by m or C-m: that would count the same mixing entropy twice.

At identical frozen site/electrode potentials, the canonical creation pair has
Delta = Delta_oxide + F(m+1)-F(m) + q*(phi_site-V_electrode).  The shared saddle is
(F_initial+F_final)/2+B0, so the two barriers are B0 +/- Delta/2 with one shared
attempt frequency.  Negative barriers are unsupported and raise ModelDomainError;
there is no clipping.  This guarantees the paired rate ratio at that frozen field.
If a driven Kirchhoff potential is recomputed after changing the vacancy state,
these expressions do NOT establish global detailed balance or a global electric
free energy for the coupled process.  Exact equilibrium tests use zero bias or a
prescribed potential independent of the reacting state.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
from numbers import Integral, Real

import numpy as np

from .model import KB, ModelDomainError


SYNTHETIC_PROVENANCE = (
    "SYNTHETIC NON-LITERATURE INTERFACE MODEL. Formation offsets, finite-reservoir "
    "storage parameters and exchange barriers are illustrative assumptions; "
    "no DFT/NEB interface calibration or electrode-specific redox chemistry was supplied."
)
_RESERVOIR_FIELDS = {
    "capacity", "initial_oxygen", "mu0_eV", "kappa_eV", "barrier_eV", "attempt_Hz",
    "enabled", "mixing_entropy",
}
_EXCHANGE_FIELDS = {
    "vacancy_formation_offset_eV", "reference_basis", "bottom", "top",
    "calibration_status", "provenance",
}


def _finite(value, name, *, positive=False, nonnegative=False):
    if not isinstance(value, Real) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite real number, not a boolean or string")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _integer(value, name, *, minimum=None, maximum=None):
    if not isinstance(value, Integral) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer, not a boolean or fractional value")
    result = int(value)
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    if maximum is not None and result > maximum:
        raise ValueError(f"{name} must be at most {maximum}")
    return result


def _fields(raw, expected, name):
    if not isinstance(raw, dict):
        raise ValueError(f"{name} must be an object")
    missing, unknown = set(expected) - set(raw), set(raw) - set(expected)
    if missing:
        raise ValueError(f"{name} is missing fields: {', '.join(sorted(missing))}")
    if unknown:
        raise ValueError(f"{name} has unknown fields: {', '.join(sorted(map(str, unknown)))}")


def _thermal_energy(temperature_K):
    temperature = _finite(temperature_K, "temperature_K", positive=True)
    kt = KB * temperature
    if not math.isfinite(kt) or kt <= 0:
        raise ModelDomainError("k_B*T is outside representable positive energy scales")
    return kt


def _direction(direction):
    direction = _integer(direction, "direction")
    if direction not in (-1, 1):
        raise ValueError("direction must be +1 (oxide oxygen out) or -1 (oxygen in)")
    return direction


@dataclass(frozen=True)
class ReservoirConfig:
    capacity: int = 64
    initial_oxygen: int = 0
    mu0_eV: float = 0.0
    kappa_eV: float = 0.2
    barrier_eV: float = 0.85
    attempt_Hz: float = 1e12
    enabled: bool = True
    mixing_entropy: bool = True

    def __post_init__(self):
        capacity = _integer(self.capacity, "capacity", minimum=1)
        initial = _integer(self.initial_oxygen, "initial_oxygen", minimum=0, maximum=capacity)
        object.__setattr__(self, "capacity", capacity)
        object.__setattr__(self, "initial_oxygen", initial)
        for name in ("mu0_eV", "kappa_eV", "barrier_eV", "attempt_Hz"):
            object.__setattr__(self, name, _finite(getattr(self, name), name,
                                                 positive=name == "attempt_Hz",
                                                 nonnegative=name == "barrier_eV"))
        for name in ("enabled", "mixing_entropy"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a boolean")

    @classmethod
    def from_dict(cls, raw):
        _fields(raw, _RESERVOIR_FIELDS, "reservoir")
        return cls(**raw)

    def to_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class ExchangeConfig:
    vacancy_formation_offset_eV: float = 0.4
    reference_basis: int = 0
    bottom: ReservoirConfig = field(default_factory=lambda: ReservoirConfig(enabled=False))
    top: ReservoirConfig = field(default_factory=ReservoirConfig)
    calibration_status: str = "illustrative"
    provenance: str = SYNTHETIC_PROVENANCE

    def __post_init__(self):
        object.__setattr__(self, "vacancy_formation_offset_eV",
                           _finite(self.vacancy_formation_offset_eV, "vacancy_formation_offset_eV"))
        object.__setattr__(self, "reference_basis", _integer(self.reference_basis, "reference_basis", minimum=0, maximum=7))
        if not isinstance(self.bottom, ReservoirConfig) or not isinstance(self.top, ReservoirConfig):
            raise ValueError("bottom and top must be ReservoirConfig instances; use ExchangeConfig.from_dict for JSON")
        if self.calibration_status != "illustrative":
            raise ValueError("this exchange schema accepts only illustrative synthetic parameters; DFT calibration cannot be invented")
        if not isinstance(self.provenance, str) or not self.provenance.strip() or "SYNTHETIC" not in self.provenance.upper():
            raise ValueError("exchange provenance must be a nonempty explicit SYNTHETIC statement")

    @classmethod
    def from_dict(cls, raw):
        _fields(raw, _EXCHANGE_FIELDS, "exchange")
        return cls(vacancy_formation_offset_eV=raw["vacancy_formation_offset_eV"],
                   reference_basis=raw["reference_basis"],
                   bottom=ReservoirConfig.from_dict(raw["bottom"]),
                   top=ReservoirConfig.from_dict(raw["top"]),
                   calibration_status=raw["calibration_status"], provenance=raw["provenance"])

    def to_dict(self):
        return asdict(self)


def _creation_delta(config, temperature_K, count_before):
    kt = _thermal_energy(temperature_K)
    m = _integer(count_before, "count_before", minimum=0, maximum=config.capacity - 1)
    try:
        result = config.mu0_eV + config.kappa_eV * (m / config.capacity)
        if config.mixing_entropy:
            # Exact finite-C ratio, including m=0 and m=C-1; no log(0) at a wall.
            result -= kt * (math.log(config.capacity - m) - math.log(m + 1))
    except (ValueError, OverflowError) as exc:
        raise ModelDomainError("reservoir free-energy increment is outside the numerical domain") from exc
    if not math.isfinite(result):
        raise ModelDomainError("reservoir free-energy increment is non-finite")
    return result


class FiniteReservoir:
    """Finite, well-mixed oxygen count with an auditable transfer ledger.

    Thermodynamic queries may inspect an otherwise disabled reservoir.  Actual
    availability and ``apply`` always respect ``enabled`` and both capacity walls.
    Oxygen atom identities are intentionally left to the coupled simulator.
    """

    def __init__(self, config: ReservoirConfig):
        if not isinstance(config, ReservoirConfig):
            raise ValueError("config must be ReservoirConfig")
        self.config = config
        self._initial_count = config.initial_oxygen
        self._count = config.initial_oxygen
        self._out_count = 0
        self._in_count = 0

    @property
    def count(self):
        return self._count

    @property
    def initial_count(self):
        return self._initial_count

    @property
    def out_count(self):
        return self._out_count

    @property
    def in_count(self):
        return self._in_count

    def available(self, direction):
        direction = _direction(direction)
        return self.config.enabled and (self.count < self.config.capacity if direction == 1 else self.count > 0)

    def apply(self, direction):
        direction = _direction(direction)
        if not self.available(direction):
            reason = "disabled" if not self.config.enabled else ("full" if direction == 1 else "empty")
            raise ValueError(f"reservoir is {reason}; oxygen transfer is unavailable")
        self._count += direction
        if direction == 1:
            self._out_count += 1
        else:
            self._in_count += 1
        if self.count != self.initial_count + self.out_count - self.in_count:
            raise RuntimeError("reservoir oxygen ledger is inconsistent")

    def free_energy_eV(self, T, count=None):
        kt = _thermal_energy(T)
        m = self.count if count is None else _integer(count, "count", minimum=0, maximum=self.config.capacity)
        c = self.config.capacity
        try:
            result = self.config.mu0_eV * m + 0.5 * self.config.kappa_eV * (m * (m - 1) / c)
            if self.config.mixing_entropy and 0 < m < c:
                log_multiplicity = math.lgamma(c + 1) - math.lgamma(m + 1) - math.lgamma(c - m + 1)
                result -= kt * log_multiplicity
        except (ValueError, OverflowError) as exc:
            raise ModelDomainError("reservoir free energy is outside the numerical domain") from exc
        if not math.isfinite(result):
            raise ModelDomainError("reservoir free energy is non-finite")
        return result

    def creation_delta_eV(self, T, count_before=None):
        return _creation_delta(self.config, T, self.count if count_before is None else count_before)


def exchange_pair(oxide_creation_energy_eV, reservoir_config, count_before, charge,
                  site_potential_V, electrode_potential_V, temperature_K):
    """Return both directions of one canonical oxygen-out reaction pair.

    ``count_before`` is the reservoir occupancy in the oxygen-filled-oxide state.
    For an oxygen-in event from current occupancy m, pass count_before=m-1 and
    select in_rate_Hz/in_barrier_eV; its signed state/work increments are the
    negatives of the returned canonical creation deltas.  No rate multiplicity
    must be applied outside this function when mixing_entropy=True.
    """
    if not isinstance(reservoir_config, ReservoirConfig):
        raise ValueError("reservoir_config must be ReservoirConfig")
    if not reservoir_config.enabled:
        raise ValueError("reservoir exchange is disabled")
    oxide_delta = _finite(oxide_creation_energy_eV, "oxide_creation_energy_eV")
    q = _integer(charge, "charge")
    phi = _finite(site_potential_V, "site_potential_V")
    voltage = _finite(electrode_potential_V, "electrode_potential_V")
    kt = _thermal_energy(temperature_K)
    reservoir_delta = _creation_delta(reservoir_config, temperature_K, count_before)
    try:
        chemical_delta = oxide_delta + reservoir_delta
        electrical_work = q * (phi - voltage)
        delta = chemical_delta + electrical_work
        out_barrier = reservoir_config.barrier_eV + 0.5 * delta
        in_barrier = reservoir_config.barrier_eV - 0.5 * delta
    except (ValueError, OverflowError) as exc:
        raise ModelDomainError("exchange energies are outside the numerical domain") from exc
    values = (chemical_delta, electrical_work, delta, out_barrier, in_barrier)
    if not all(math.isfinite(x) for x in values):
        raise ModelDomainError("exchange work, free-energy delta and activation barriers must be finite")
    if out_barrier < 0 or in_barrier < 0:
        raise ModelDomainError(
            f"negative exchange barrier: out={out_barrier:.12g} eV, in={in_barrier:.12g} eV; "
            "formation, storage or frozen-field drive exceeds the shared-saddle domain (no clipping)"
        )
    out_rate = reservoir_config.attempt_Hz * math.exp(-out_barrier / kt)
    in_rate = reservoir_config.attempt_Hz * math.exp(-in_barrier / kt)
    return {
        "chemical_delta_eV": chemical_delta,
        "oxide_delta_eV": oxide_delta,
        "reservoir_delta_eV": reservoir_delta,
        "electrical_work_eV": electrical_work,
        "delta_eV": delta,
        "out_barrier_eV": out_barrier,
        "in_barrier_eV": in_barrier,
        "out_rate_Hz": out_rate,
        "in_rate_Hz": in_rate,
    }
