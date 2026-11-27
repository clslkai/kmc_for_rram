"""Planar blocking electrodes around a fixed-lattice oxide slab.

The A/B cell faces are periodic.  The two electrodes bound the finite C
coordinate, so their unit normal is (A x B)/|A x B| and their separation is
C dot normal, not |C|.  Electrode planes are ideal potential boundary conditions;
no metal atoms, electronic transport, oxygen exchange or Poisson solver are
introduced here.  Contact masks label oxide sites near the two planes.  Those
sites retain the linear oxide potential, rather than becoming metal nodes.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import math
from numbers import Integral, Real
from pathlib import Path
import tempfile

import numpy as np


_DEPTH_ATOL_A = 1e-9


def _finite_real(value, name, *, positive=False):
    if not isinstance(value, Real) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite real number, not a boolean or string")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(value) or (positive and value <= 0):
        qualifier = "finite and positive" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return value


def _exact_fields(raw, fields, name):
    if not isinstance(raw, dict):
        raise ValueError(f"{name} must be a JSON object")
    missing, unknown = set(fields) - set(raw), set(raw) - set(fields)
    if missing:
        raise ValueError(f"{name} is missing fields: {', '.join(sorted(missing))}")
    if unknown:
        raise ValueError(f"{name} has unknown fields: {', '.join(sorted(map(str, unknown)))}")


def _json_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _json_constant(value):
    raise ValueError(f"non-finite JSON constant is forbidden: {value}")


def _array(value, name, shape=None):
    try:
        raw = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite numeric array") from exc
    if shape is not None and raw.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    for entry in raw.flat:
        _finite_real(entry, name)
    result = np.asarray(raw, dtype=float)
    return result


def _readonly(value, dtype=None):
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class ElectrodeConfig:
    potential_V: float = 0.0
    contact_depth_A: float = 2.0
    oxygen_boundary: str = "blocking"

    def __post_init__(self):
        object.__setattr__(self, "potential_V", _finite_real(self.potential_V, "potential_V"))
        object.__setattr__(self, "contact_depth_A", _finite_real(self.contact_depth_A, "contact_depth_A", positive=True))
        if self.oxygen_boundary != "blocking":
            raise ValueError("only oxygen_boundary='blocking' is implemented; interface oxygen exchange requires a separate calibrated model")


@dataclass(frozen=True)
class DeviceConfig:
    bottom: ElectrodeConfig = field(default_factory=ElectrodeConfig)
    top: ElectrodeConfig = field(default_factory=lambda: ElectrodeConfig(potential_V=1.0))
    potential_model: str = "linear_dielectric"
    schema_version: int = 1

    def __post_init__(self):
        if not isinstance(self.bottom, ElectrodeConfig) or not isinstance(self.top, ElectrodeConfig):
            raise ValueError("bottom and top must be ElectrodeConfig instances; use DeviceConfig.from_dict for JSON")
        if self.potential_model != "linear_dielectric":
            raise ValueError("only potential_model='linear_dielectric' is implemented; no self-consistent Poisson model is available")
        if (not isinstance(self.schema_version, Integral) or isinstance(self.schema_version, (bool, np.bool_))
                or self.schema_version != 1):
            raise ValueError("device schema_version must be integer 1")
        object.__setattr__(self, "schema_version", int(self.schema_version))
        if not math.isfinite(self.top.potential_V - self.bottom.potential_V):
            raise ValueError("top-minus-bottom electrode potential difference must be finite")

    @classmethod
    def from_dict(cls, raw):
        """Parse the complete version-1 schema; no missing electrode bias is guessed."""
        _exact_fields(raw, {"bottom", "top", "potential_model", "schema_version"}, "device")
        electrodes = {}
        for name in ("bottom", "top"):
            _exact_fields(raw[name], {"potential_V", "contact_depth_A", "oxygen_boundary"}, name)
            electrodes[name] = ElectrodeConfig(**raw[name])
        return cls(bottom=electrodes["bottom"], top=electrodes["top"],
                   potential_model=raw["potential_model"], schema_version=raw["schema_version"])

    @classmethod
    def load(cls, path):
        with Path(path).open(encoding="utf-8") as handle:
            raw = json.load(handle, object_pairs_hook=_json_pairs, parse_constant=_json_constant)
        return cls.from_dict(raw)

    def to_dict(self):
        return asdict(self)

    def write(self, path):
        """Atomically write the complete validated configuration as strict JSON."""
        payload = json.dumps(self.to_dict(), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=destination.parent,
                                             prefix=destination.name + ".", suffix=".tmp", delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(payload)
            temporary.replace(destination)
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()


class DeviceBoundary:
    """Geometry and imposed linear electrostatic potential of a blocking slab.

    Coordinates are Cartesian angstrom; unit and supercell vectors are rows.
    ``relative_potential_V`` uses the bottom electrode as zero, whereas
    ``site_potential_V`` and ``potential_V`` include its configured absolute
    potential.  The static field depends only on top-minus-bottom voltage.
    No lattice sites, migration edges or interaction pairs are modified.
    """

    def __init__(self, lattice, config: DeviceConfig):
        if not isinstance(config, DeviceConfig):
            raise ValueError("config must be DeviceConfig")
        if getattr(lattice, "boundary", None) != "slab" or tuple(getattr(lattice, "pbc", ())) != (True, True, False):
            raise ValueError("DeviceBoundary requires a slab with A/B periodic and C finite (pbc=True,True,False)")
        self.lattice, self.config = lattice, config
        self._cell_A = _readonly(_array(lattice.cell_A, "lattice.cell_A", (3, 3)))
        a, b, c = self._cell_A
        cross = np.cross(a, b)
        area = float(np.linalg.norm(cross))
        if not math.isfinite(area) or area <= 0:
            raise ValueError("electrode-plane area must be finite and positive")
        self.area_A2 = area
        self.normal = _readonly(cross / area)
        self.thickness_A = float(c @ self.normal)
        if not math.isfinite(self.thickness_A) or self.thickness_A <= 0:
            raise ValueError("C dot electrode normal must be finite and positive")
        self._depth_tolerance_A = _DEPTH_ATOL_A + 1e-12 * self.thickness_A
        if config.bottom.contact_depth_A + config.top.contact_depth_A >= self.thickness_A:
            raise ValueError("bottom and top contact depths must sum to less than the oxide thickness; overlapping contact layers are unsupported")
        positions = _array(lattice.oxygen_A, "lattice.oxygen_A")
        if positions.ndim != 2 or positions.shape != (lattice.n, 3) or lattice.n < 1:
            raise ValueError("lattice.oxygen_A must contain one Cartesian position per nonempty oxide site")
        self.depth_A = _readonly(self._depths(positions))
        bottom = self.depth_A <= config.bottom.contact_depth_A + self._depth_tolerance_A
        top = self.thickness_A - self.depth_A <= config.top.contact_depth_A + self._depth_tolerance_A
        if not np.any(bottom) or not np.any(top):
            raise ValueError("each electrode contact region must contain oxide sites; empty layers are not automatically enlarged")
        if np.any(bottom & top):
            raise ValueError("bottom and top contact site sets overlap")
        self.bottom_mask, self.top_mask = _readonly(bottom), _readonly(top)
        regions = np.full(lattice.n, "interior", dtype="<U14")
        regions[bottom], regions[top] = "bottom_contact", "top_contact"
        self.region = _readonly(regions)
        voltage = config.top.potential_V - config.bottom.potential_V
        self.relative_potential_V = _readonly(voltage * (self.depth_A / self.thickness_A))
        self.site_potential_V = _readonly(config.bottom.potential_V + self.relative_potential_V)
        self.field_V_A = _readonly(-(voltage / self.thickness_A) * self.normal)
        self.field_V_nm = _readonly(10.0 * self.field_V_A)
        if not (np.isfinite(self.relative_potential_V).all() and np.isfinite(self.site_potential_V).all()
                and np.isfinite(self.field_V_A).all() and np.isfinite(self.field_V_nm).all()):
            raise ValueError("derived device potentials and field must be finite")

    def _depths(self, positions_A):
        depth = np.asarray(positions_A, dtype=float) @ self.normal
        if np.any(depth < -self._depth_tolerance_A) or np.any(depth > self.thickness_A + self._depth_tolerance_A):
            raise ValueError("positions lie outside the finite slab between its electrode planes")
        # Only remove roundoff at a face, never project a real out-of-slab site.
        return np.clip(depth, 0.0, self.thickness_A)

    def potential_V(self, positions_A):
        """Absolute potential at a point (3,) or points (N,3) inside the slab.

        Lateral periodic images are valid inputs: translating by A or B leaves
        potential unchanged.  Positions outside the finite normal interval fail.
        """
        positions = _array(positions_A, "positions_A")
        scalar = positions.shape == (3,)
        if not scalar and (positions.ndim != 2 or positions.shape[1] != 3):
            raise ValueError("positions_A must have shape (3,) or (N,3)")
        depth = self._depths(positions)
        voltage = self.config.top.potential_V - self.config.bottom.potential_V
        result = self.config.bottom.potential_V + voltage * (depth / self.thickness_A)
        return float(result) if scalar else np.asarray(result, dtype=float)

    def electrostatic_energy_eV(self, vacancy, q):
        """External energy q*sum(phi_i-phi_bottom), in eV, for occupied vacancies.

        q is an integer charge in units of positive elementary charge.  This
        helper neither adds a 1/2 factor nor includes self-energy/space charge.
        Using the relative potential makes its reported value gauge invariant.
        """
        if not isinstance(q, Integral) or isinstance(q, (bool, np.bool_)):
            raise ValueError("q must be an integral charge in units of elementary charge")
        occupancy = np.asarray(vacancy)
        if occupancy.shape != (self.lattice.n,) or occupancy.dtype.kind != "b":
            raise ValueError("vacancy must be a boolean occupancy vector with one entry per oxide site")
        if q == 0 or not np.any(occupancy):
            return 0.0
        try:
            value = float(int(q) * math.fsum(float(x) for x in self.relative_potential_V[occupancy]))
        except (ValueError, OverflowError) as exc:
            raise ValueError("external electrostatic energy is non-finite") from exc
        if not math.isfinite(value):
            raise ValueError("external electrostatic energy is non-finite")
        return value

    def summary(self):
        """Serializable boundary metadata; contact tags do not imply conduction."""
        bottom_count, top_count = int(self.bottom_mask.sum()), int(self.top_mask.sum())
        return {
            "model": "electrode_oxide_electrode_boundary",
            "config": self.config.to_dict(),
            "boundary": "slab", "pbc": [True, True, False],
            "boundary_convention": "Periodic translations along cell A and B; finite C faces at projected depths 0 and C dot normal.",
            "normal": self.normal.tolist(), "thickness_A": self.thickness_A,
            "thickness_nm": self.thickness_A / 10.0, "area_A2": self.area_A2,
            "applied_voltage_V": self.config.top.potential_V - self.config.bottom.potential_V,
            "voltage_convention": "V_applied=V_top-V_bottom; E=-grad(phi). Positive top bias drives positive vacancies toward the bottom.",
            "field_V_A": self.field_V_A.tolist(), "field_V_nm": self.field_V_nm.tolist(),
            "potential_reference": "Absolute potentials include V_bottom; relative site potentials and external energy use V_bottom=0 as reference.",
            "planes": {
                "bottom": {"point_A": [0.0, 0.0, 0.0], "depth_A": 0.0, "potential_V": self.config.bottom.potential_V},
                "top": {"point_A": self._cell_A[2].tolist(), "depth_A": self.thickness_A, "potential_V": self.config.top.potential_V},
                "span_A": self._cell_A[0].tolist(), "span_B": self._cell_A[1].tolist(),
            },
            "contacts": {
                "role": "Oxide-adjacent site tags only, not equipotential metal nodes or calibrated interface states.",
                "bottom": {"contact_depth_A": self.config.bottom.contact_depth_A,
                           "n_sites": bottom_count, "site_indices": np.flatnonzero(self.bottom_mask).tolist(),
                           "oxygen_boundary": "blocking"},
                "top": {"contact_depth_A": self.config.top.contact_depth_A,
                        "n_sites": top_count, "site_indices": np.flatnonzero(self.top_mask).tolist(),
                        "oxygen_boundary": "blocking"},
            },
            "counts": {"oxide_sites": int(self.lattice.n), "bottom_contact": bottom_count,
                       "interior": int(self.lattice.n) - bottom_count - top_count, "top_contact": top_count},
            "explicit_metal_atoms": False, "oxygen_exchange_enabled": False,
            "electronic_current_model": None, "poisson_solver": None, "space_charge_feedback": False,
            "interface_chemical_calibration": None,
            "assumptions": [
                "Ideal planar electrodes impose a static linear dielectric potential; microscopic electrode material and work-function effects are unspecified.",
                "Both interfaces block oxygen exchange; the boundary adds no vacancy creation, removal or charge-state changes.",
                "The same classical catalogue is used up to the truncated oxide surfaces; contact tags do not supply metal-specific interface energies or barriers.",
                "The prescribed field has no defect-charge feedback, image-charge correction or self-consistent Poisson solution.",
                "No electronic transport, current, heating or switching-filament prediction is made by this boundary module.",
                "For fixed vacancy count and charge, a common shift of both electrode potentials leaves rates and relative external energy unchanged.",
            ],
        }
