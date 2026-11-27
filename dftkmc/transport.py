"""Illustrative isothermal electronic resistor network on an oxide slab.

This graph is constructed from distances, independently of the ionic migration
and vacancy-interaction graphs.  It is a classical effective conductance model,
not a DFT-calibrated transport law, a Poisson solver, quantum transport, or a
temperature solver.  Contact *oxide* sites retain their solved potential and
connect to ideal metal nodes through finite resistors.

Only translations along cell A and B are periodic.  Cell vectors are rows and
all distances are in angstrom.  An interior conductance is ``g0*r_ref/r``; g0 is
the vacancy value only when BOTH sites contain vacancies.  A contact uses its
site's vacancy state and the effective length ``plane_distance+contact_offset``.
The strictly positive background admits leakage without a vacancy-only path.

Local potentials may be used by a separately documented driven KMC closure.
Because conductances depend on the defect configuration, this module does not
claim that q*sum(n_i*phi_i) is a thermodynamic state energy or that recomputing
the field after an ionic event establishes global detailed balance.
"""
from __future__ import annotations

from collections import deque
from dataclasses import asdict, dataclass, fields
import math
from numbers import Real
import warnings

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from .device import DeviceBoundary


def _finite_real(value, name, *, positive=False):
    if not isinstance(value, Real) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite real number, not a boolean or string")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return result


def _readonly(value, dtype=None):
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _finite_array(value, name, shape):
    raw = np.asarray(value)
    if raw.shape != shape or raw.dtype.kind not in "iuf" or not np.isfinite(raw).all():
        raise ValueError(f"{name} must be a finite numeric array of shape {shape}")
    return _readonly(raw, dtype=float)


@dataclass(frozen=True)
class TransportConfig:
    """Complete parameters of the illustrative distance-resistor law.

    ``read_voltage_V`` is a positive, independent frozen-configuration probe.
    Its linear current is reported without applying the write-current limiter;
    no ionic evolution or physical read pulse is performed.
    """

    electronic_cutoff_A: float = 4.0
    background_conductance_S: float = 1e-12
    vacancy_conductance_S: float = 1e-6
    reference_distance_A: float = 1.0
    contact_offset_A: float = 1.0
    current_compliance_A: float | None = None
    read_voltage_V: float = 0.1
    calibration_status: str = "illustrative"
    provenance: str = "Synthetic effective resistor parameters; not DFT- or measurement-calibrated."

    def __post_init__(self):
        for name in ("electronic_cutoff_A", "background_conductance_S", "vacancy_conductance_S",
                     "reference_distance_A", "contact_offset_A", "read_voltage_V"):
            object.__setattr__(self, name, _finite_real(getattr(self, name), name, positive=True))
        if self.vacancy_conductance_S < self.background_conductance_S:
            raise ValueError("vacancy_conductance_S must be at least background_conductance_S")
        if self.current_compliance_A is not None:
            object.__setattr__(self, "current_compliance_A",
                               _finite_real(self.current_compliance_A, "current_compliance_A", positive=True))
        if self.calibration_status != "illustrative":
            raise ValueError("this effective resistor law currently supports calibration_status='illustrative' only")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("transport provenance must be a nonempty string")

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, raw):
        if not isinstance(raw, dict):
            raise ValueError("transport configuration must be a JSON object")
        required = {item.name for item in fields(cls)}
        missing, unknown = required - set(raw), set(raw) - required
        if missing:
            raise ValueError("transport configuration is missing fields: " + ", ".join(sorted(missing)))
        if unknown:
            raise ValueError("transport configuration has unknown fields: " + ", ".join(sorted(map(str, unknown))))
        return cls(**raw)


@dataclass(frozen=True)
class TransportSolution:
    """One frozen-defect DC solution, with explicit current orientations.

    ``terminal_current_A`` is positive from the TOP electrode into the oxide;
    ``bottom_terminal_current_A`` is positive from BOTTOM into the oxide.  Thus
    their sum is zero in steady state.  Interior currents are source->target
    according to ``network.pairs``.  Every contact current is electrode->site.
    ``edge_field_V_nm`` is the signed potential gradient along source->target,
    not a three-dimensional local electric-field reconstruction.
    """

    potential_V: np.ndarray
    relative_potential_V: np.ndarray
    requested_voltage_V: float
    applied_voltage_V: float
    bottom_potential_V: float
    top_potential_V: float
    equivalent_conductance_S: float
    terminal_current_A: float
    bottom_terminal_current_A: float
    read_voltage_V: float
    read_current_A: float
    current_compliance_A: float | None
    compliance_active: bool
    edge_conductance_S: np.ndarray
    edge_current_A: np.ndarray
    edge_power_W: np.ndarray
    edge_field_V_nm: np.ndarray
    bottom_contact_conductance_S: np.ndarray
    top_contact_conductance_S: np.ndarray
    bottom_contact_current_A: np.ndarray
    top_contact_current_A: np.ndarray
    bottom_contact_power_W: np.ndarray
    top_contact_power_W: np.ndarray
    node_kcl_residual_A: np.ndarray
    max_abs_kcl_residual_A: float
    relative_kcl_residual: float
    terminal_current_imbalance_A: float
    relative_terminal_current_imbalance: float
    total_power_W: float
    power_balance_residual_W: float
    relative_power_balance_residual: float
    terminal_network_connected: bool
    vacancy_graph_spanning: bool
    conducting_path_site_ids: tuple[int, ...]
    calibration_status: str
    provenance: str

    @property
    def top_terminal_current_A(self):
        return self.terminal_current_A

    @property
    def vacancy_spanning(self):
        """Alias; this is a graph property, not a calibrated electronic ON test."""
        return self.vacancy_graph_spanning

    def summary(self):
        result = {}
        for item in fields(self):
            value = getattr(self, item.name)
            if not isinstance(value, np.ndarray):
                result[item.name] = list(value) if isinstance(value, tuple) else value
        result.update(
            model="illustrative_isothermal_distance_resistor_network",
            current_convention="Top and bottom terminal/contact currents are positive electrode->oxide; internal edges use source->target.",
            edge_field_convention="Signed E along source->target: 10*(phi_source-phi_target)/distance_A, in V/nm.",
            path_definition="A bottom-to-top path of vacancy sites on the independent electronic distance graph; not an experimental or calibrated filament criterion.",
            read_current_definition="G of this frozen defect state times positive read_voltage_V; no ionic evolution, heating, or read-current compliance is applied.",
            electronic_on_threshold=None,
            potential_model="Kirchhoff resistor network; not space-charge Poisson",
            isothermal=True,
            temperature_feedback=False,
            quantum_transport=False,
        )
        return result

    def to_dict(self, include_arrays=True):
        result = self.summary()
        if include_arrays:
            result.update({item.name: getattr(self, item.name).tolist()
                           for item in fields(self) if isinstance(getattr(self, item.name), np.ndarray)})
        return result


class ElectronicNetwork:
    """Independent slab graph and two-terminal Dirichlet Kirchhoff solver.

    The positive graph is fixed by geometry.  A component attached to exactly
    one electrode is equipotential with it and carries no current.  Components
    attached to neither electrode are rejected: their DC potential is not
    defined and must not be silently replaced by a guessed field for KMC.
    """

    def __init__(self, lattice, device: DeviceBoundary, config: TransportConfig | None = None):
        if not isinstance(device, DeviceBoundary) or device.lattice is not lattice:
            raise ValueError("device must be a DeviceBoundary bound to this lattice")
        config = TransportConfig() if config is None else config
        if not isinstance(config, TransportConfig):
            raise ValueError("config must be TransportConfig")
        self.lattice, self.device, self.config = lattice, device, config
        self.n = int(lattice.n)
        self.cell_A = _finite_array(lattice.cell_A, "lattice.cell_A", (3, 3))
        self.positions_A = _finite_array(lattice.oxygen_A, "lattice.oxygen_A", (self.n, 3))
        self.bottom_sites = _readonly(np.flatnonzero(device.bottom_mask), dtype=int)
        self.top_sites = _readonly(np.flatnonzero(device.top_mask), dtype=int)
        self.bottom_contact_distances_A = _readonly(device.depth_A[self.bottom_sites])
        self.top_contact_distances_A = _readonly(device.thickness_A - device.depth_A[self.top_sites])
        lateral = self.cell_A[:2]
        self.unique_image_bound_A = float(np.linalg.svd(lateral, compute_uv=False).min()) / 2
        tolerance = 1e-10 * max(1.0, config.electronic_cutoff_A)
        if config.electronic_cutoff_A + tolerance >= self.unique_image_bound_A:
            raise ValueError("lateral supercell too small for unique electronic periodic images; increase A/B repeats or reduce electronic_cutoff_A")
        self.pairs, self.displacements_A, self.distances_A = self._build_graph(lateral, tolerance)
        # Only independent electronic pairs enter this adjacency.  Ionic
        # migration and interaction attributes are deliberately never read.
        self._neighbors = [[] for _ in range(self.n)]
        for i, j in self.pairs:
            self._neighbors[i].append(int(j))
            self._neighbors[j].append(int(i))
        self._classify_components()
        with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
            self._edge_geometry_factor = config.reference_distance_A / self.distances_A
            self._bottom_geometry_factor = config.reference_distance_A / (self.bottom_contact_distances_A + config.contact_offset_A)
            self._top_geometry_factor = config.reference_distance_A / (self.top_contact_distances_A + config.contact_offset_A)
        # The entire allowable conductance interval must remain representable.
        for factor in (self._edge_geometry_factor, self._bottom_geometry_factor, self._top_geometry_factor):
            for base in (config.background_conductance_S, config.vacancy_conductance_S):
                with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                    values = factor * base
                if not np.isfinite(values).all() or np.any(values <= 0):
                    raise ValueError("derived electronic conductances must be finite and strictly positive")

    def _build_graph(self, lateral, tolerance):
        # Project to the lateral span instead of wrapping Cartesian x/y or C.
        # The cutoff < sigma_min(A,B)/2 bound ensures that any accepted image
        # lies inside the coefficient rounding cell, including skew A/B cells.
        projection = lateral.T @ np.linalg.inv(lateral @ lateral.T)
        pairs, displacements, distances = [], [], []
        for i in range(self.n - 1):
            raw = self.positions_A[i + 1:] - self.positions_A[i]
            image = np.rint(raw @ projection)
            delta = raw - image @ lateral
            length = np.linalg.norm(delta, axis=1)
            if np.any(length <= 1e-12):
                raise ValueError("duplicate oxide positions under lateral periodicity are not supported")
            for offset in np.flatnonzero(length <= self.config.electronic_cutoff_A + tolerance):
                pairs.append((i, int(i + 1 + offset)))
                displacements.append(delta[offset])
                distances.append(length[offset])
        return (_readonly(np.asarray(pairs, dtype=int).reshape(-1, 2)),
                _readonly(np.asarray(displacements, dtype=float).reshape(-1, 3)),
                _readonly(distances, dtype=float))

    def _classify_components(self):
        i, j = self.pairs.T
        adjacency = coo_matrix((np.ones(2 * len(i)), (np.r_[i, j], np.r_[j, i])),
                               shape=(self.n, self.n)).tocsr()
        count, labels = connected_components(adjacency, directed=False)
        bottom = np.zeros(count, dtype=bool)
        top = np.zeros(count, dtype=bool)
        bottom[labels[self.bottom_sites]], top[labels[self.top_sites]] = True, True
        floating = ~bottom & ~top
        if np.any(floating):
            nodes = np.flatnonzero(floating[labels]).tolist()
            raise ValueError(f"electronic graph has floating oxide component(s) with no electrode contact; potentials undefined at site IDs {nodes}")
        self.n_components = int(count)
        self.terminal_network_connected = bool(np.any(bottom & top))
        self._bridging_nodes = np.flatnonzero((bottom & top)[labels])
        self._bottom_only_nodes = np.flatnonzero((bottom & ~top)[labels])
        self._top_only_nodes = np.flatnonzero((top & ~bottom)[labels])
        self._labels = labels

    def _vacancy_path(self, vacancy):
        targets = set(int(i) for i in self.top_sites if vacancy[i])
        starts = [int(i) for i in self.bottom_sites if vacancy[i]]
        queue, predecessor = deque(starts), {i: None for i in starts}
        while queue:
            i = queue.popleft()
            if i in targets:
                path = [i]
                while predecessor[path[-1]] is not None:
                    path.append(predecessor[path[-1]])
                return tuple(reversed(path))
            for j in self._neighbors[i]:
                if vacancy[j] and j not in predecessor:
                    predecessor[j] = i
                    queue.append(j)
        return ()

    def _unit_solution(self, edge_g, bottom_g, top_g):
        """Solve bottom=0, top=1 independent of applied bias and compliance."""
        i, j = self.pairs.T
        unit = np.zeros(self.n)
        unit[self._top_only_nodes] = 1.0
        if not len(self._bridging_nodes):
            return unit
        # Scaling all conductances leaves potentials unchanged and avoids a
        # needless small dimensional scale in the sparse linear system.
        scale = float(max(np.max(bottom_g), np.max(top_g), np.max(edge_g, initial=0.0)))
        ge, gb, gt = edge_g / scale, bottom_g / scale, top_g / scale
        diagonal = np.zeros(self.n)
        np.add.at(diagonal, i, ge)
        np.add.at(diagonal, j, ge)
        np.add.at(diagonal, self.bottom_sites, gb)
        np.add.at(diagonal, self.top_sites, gt)
        row = np.r_[np.arange(self.n), i, j]
        col = np.r_[np.arange(self.n), j, i]
        matrix = coo_matrix((np.r_[diagonal, -ge, -ge], (row, col)), shape=(self.n, self.n)).tocsr()
        rhs = np.zeros(self.n)
        np.add.at(rhs, self.top_sites, gt)
        active = self._bridging_nodes
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", MatrixRankWarning)
                values = spsolve(matrix[active][:, active].tocsc(), rhs[active])
        except (MatrixRankWarning, RuntimeError) as exc:
            raise ValueError("electronic Kirchhoff system could not be solved") from exc
        if not np.isfinite(values).all():
            raise ValueError("electronic Kirchhoff solution is non-finite")
        if np.any(values < -1e-9) or np.any(values > 1 + 1e-9):
            raise ValueError("electronic Kirchhoff solution violates the Dirichlet maximum principle")
        unit[active] = values
        return unit

    def solve(self, vacancy, bottom_V=None, top_V=None):
        """Solve at the requested bias, reducing its magnitude for compliance.

        The bottom potential remains fixed.  ``applied_voltage_V`` and
        ``top_potential_V`` are the ACTUAL values after ideal current limiting;
        a coupled ionic/interface model must use these returned values.
        """
        vacancy = np.asarray(vacancy)
        if vacancy.shape != (self.n,) or vacancy.dtype.kind != "b":
            raise ValueError("vacancy must be a boolean occupancy vector with one entry per oxide site")
        bottom = _finite_real(self.device.config.bottom.potential_V if bottom_V is None else bottom_V, "bottom_V")
        top = _finite_real(self.device.config.top.potential_V if top_V is None else top_V, "top_V")
        requested = top - bottom
        if not math.isfinite(requested):
            raise ValueError("top-minus-bottom voltage must be finite")
        i, j = self.pairs.T
        low, high = self.config.background_conductance_S, self.config.vacancy_conductance_S
        edge_g = self._edge_geometry_factor * np.where(vacancy[i] & vacancy[j], high, low)
        bottom_g = self._bottom_geometry_factor * np.where(vacancy[self.bottom_sites], high, low)
        top_g = self._top_geometry_factor * np.where(vacancy[self.top_sites], high, low)
        unit = self._unit_solution(edge_g, bottom_g, top_g)
        unit_edge_drop = unit[i] - unit[j]
        unit_bottom_drop, unit_top_drop = -unit[self.bottom_sites], 1.0 - unit[self.top_sites]
        unit_edge_current = edge_g * unit_edge_drop
        unit_bottom_current, unit_top_current = bottom_g * unit_bottom_drop, top_g * unit_top_drop
        # At unit bias the bottom-side sum avoids subtraction of numbers near
        # 1 when a high-conductance component is attached strongly to the top.
        conductance = -float(math.fsum(float(x) for x in unit_bottom_current))
        if not self.terminal_network_connected:
            conductance = 0.0
        if not math.isfinite(conductance) or (self.terminal_network_connected and conductance <= 0):
            raise ValueError("equivalent conductance must be finite and positive for a connected network")
        # Audit the unit-voltage solve even at requested zero bias, where all
        # physical currents vanish and would otherwise conceal a bad G.
        unit_kcl = np.zeros(self.n)
        np.add.at(unit_kcl, i, unit_edge_current)
        np.add.at(unit_kcl, j, -unit_edge_current)
        np.add.at(unit_kcl, self.bottom_sites, -unit_bottom_current)
        np.add.at(unit_kcl, self.top_sites, -unit_top_current)
        unit_top_terminal = float(math.fsum(float(x) for x in unit_top_current))
        unit_power = float(math.fsum(float(x) for values in
                           (unit_edge_current * unit_edge_drop, unit_bottom_current * unit_bottom_drop,
                            unit_top_current * unit_top_drop) for x in values))
        unit_scale = max(conductance, abs(unit_top_terminal),
                         float(np.max(np.abs(unit_edge_current), initial=0.0)))
        if unit_scale and max(float(np.max(np.abs(unit_kcl), initial=0.0)),
                              abs(unit_top_terminal - conductance), abs(unit_power - unit_top_terminal)) / unit_scale > 1e-7:
            raise ValueError("unit-bias electronic solve failed current/power conservation tolerance; reduce conductance contrast or improve numerical conditioning")
        voltage = requested
        limit = self.config.current_compliance_A
        compliance = False
        if limit is not None and conductance > 0:
            maximum = limit / conductance
            if abs(voltage) > maximum:
                voltage = math.copysign(maximum, voltage)
                compliance = True
        relative = voltage * unit
        potential = bottom + relative
        actual_top = bottom + voltage
        # All differences use the relative solution to retain gauge invariance
        # numerically, even when a nonzero common voltage offset is supplied.
        edge_drop = voltage * unit_edge_drop
        bottom_drop = voltage * unit_bottom_drop
        top_drop = voltage * unit_top_drop
        edge_current = voltage * unit_edge_current
        bottom_current, top_current = voltage * unit_bottom_current, voltage * unit_top_current
        edge_power = edge_current * edge_drop
        bottom_power, top_power = bottom_current * bottom_drop, top_current * top_drop
        kcl = np.zeros(self.n)
        np.add.at(kcl, i, edge_current)
        np.add.at(kcl, j, -edge_current)
        np.add.at(kcl, self.bottom_sites, -bottom_current)
        np.add.at(kcl, self.top_sites, -top_current)
        current = float(math.fsum(float(x) for x in top_current))
        bottom_terminal = float(math.fsum(float(x) for x in bottom_current))
        imbalance = current + bottom_terminal
        total_power = float(math.fsum(float(x) for values in (edge_power, bottom_power, top_power) for x in values))
        power_residual = total_power - current * voltage
        max_kcl = float(np.max(np.abs(kcl), initial=0.0))
        current_scale = max(abs(current), abs(bottom_terminal),
                            float(np.max(np.abs(edge_current), initial=0.0)),
                            float(np.max(np.abs(bottom_current), initial=0.0)),
                            float(np.max(np.abs(top_current), initial=0.0)))
        power_scale = max(abs(total_power), abs(current * voltage))
        relative_kcl = max_kcl / current_scale if current_scale else 0.0
        relative_imbalance = abs(imbalance) / current_scale if current_scale else 0.0
        relative_power = abs(power_residual) / power_scale if power_scale else 0.0
        read_current = conductance * self.config.read_voltage_V
        edge_field = 10 * edge_drop / self.distances_A
        if not (np.isfinite(potential).all() and math.isfinite(actual_top)
                and math.isfinite(current) and math.isfinite(bottom_terminal)
                and math.isfinite(total_power) and math.isfinite(read_current)
                and np.isfinite(edge_power).all() and np.isfinite(bottom_power).all()
                and np.isfinite(top_power).all() and np.isfinite(kcl).all()
                and np.isfinite(edge_field).all()):
            raise ValueError("derived transport potentials, currents and powers must be finite")
        if max(relative_kcl, relative_imbalance, relative_power) > 1e-7:
            raise ValueError("electronic solve failed current/power conservation tolerance; reduce conductance contrast or improve numerical conditioning")
        path = self._vacancy_path(vacancy)
        return TransportSolution(
            potential_V=_readonly(potential), relative_potential_V=_readonly(relative),
            requested_voltage_V=requested, applied_voltage_V=voltage,
            bottom_potential_V=bottom, top_potential_V=actual_top,
            equivalent_conductance_S=conductance, terminal_current_A=current,
            bottom_terminal_current_A=bottom_terminal,
            read_voltage_V=self.config.read_voltage_V, read_current_A=read_current,
            current_compliance_A=limit, compliance_active=compliance,
            edge_conductance_S=_readonly(edge_g), edge_current_A=_readonly(edge_current),
            edge_power_W=_readonly(edge_power), edge_field_V_nm=_readonly(edge_field),
            bottom_contact_conductance_S=_readonly(bottom_g), top_contact_conductance_S=_readonly(top_g),
            bottom_contact_current_A=_readonly(bottom_current), top_contact_current_A=_readonly(top_current),
            bottom_contact_power_W=_readonly(bottom_power), top_contact_power_W=_readonly(top_power),
            node_kcl_residual_A=_readonly(kcl), max_abs_kcl_residual_A=max_kcl,
            relative_kcl_residual=relative_kcl, terminal_current_imbalance_A=imbalance,
            relative_terminal_current_imbalance=relative_imbalance, total_power_W=total_power,
            power_balance_residual_W=power_residual, relative_power_balance_residual=relative_power,
            terminal_network_connected=self.terminal_network_connected,
            vacancy_graph_spanning=bool(path), conducting_path_site_ids=path,
            calibration_status=self.config.calibration_status, provenance=self.config.provenance,
        )

    def summary(self):
        """Static graph/law metadata; no vacancy-state conduction claim."""
        return {
            "model": "illustrative_isothermal_distance_resistor_network",
            "config": self.config.to_dict(),
            "n_sites": self.n, "n_edges": len(self.pairs), "n_components": self.n_components,
            "n_bottom_contacts": len(self.bottom_sites), "n_top_contacts": len(self.top_sites),
            "bottom_sites": self.bottom_sites.tolist(), "top_sites": self.top_sites.tolist(),
            "terminal_network_connected": self.terminal_network_connected,
            "pbc": [True, True, False], "unique_image_bound_A": self.unique_image_bound_A,
            "independent_of_migration_and_interaction_graphs": True,
            "interior_conductance_law": "g0*reference_distance_A/r_A; g0=vacancy_conductance_S only for a vacancy-vacancy edge, background_conductance_S otherwise.",
            "contact_conductance_law": "g0*reference_distance_A/(normal_plane_distance_A+contact_offset_A); g0 selected by the contact site's vacancy state.",
            "contact_nodes": "Finite contact resistors connect oxide sites to ideal electrode nodes; oxide contacts are not hard-clamped to the metal potential.",
            "floating_component_policy": "Reject components without either electrode; one-terminal components are equipotential with that electrode and carry zero current.",
            "compliance_model": "Ideal instantaneous voltage reduction at fixed bottom potential; not a series circuit or transient current controller.",
            "isothermal": True, "temperature_feedback": False, "poisson_solver": None,
            "space_charge_feedback": False, "quantum_transport": False,
            "electronic_on_threshold": None,
        }
