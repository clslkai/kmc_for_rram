"""Crystallographic labels and independent ionic/interaction graphs.

Distances are Angstrom. Cell vectors are ROWS. Channel equivalence here removes
translations and reversal only; it does not claim a space-group irreducible set.
"""
from dataclasses import dataclass
from itertools import product
import hashlib
import json
import math

import numpy as np


def default_reference():
    from atomistic.model import CELL, ASYMMETRIC, equivalents
    a, b, c, beta = CELL
    beta = np.deg2rad(beta)
    return {
        "cell_A": [[a, 0., 0.], [0., b, 0.],
                   [float(c*np.cos(beta)), 0., float(c*np.sin(beta))]],
        "hf_fractional": equivalents(ASYMMETRIC["Hf"]).tolist(),
        "oxygen_fractional": np.concatenate([
            equivalents(ASYMMETRIC["O3"]), equivalents(ASYMMETRIC["O4"])
        ]).tolist(),
        "oxygen_kinds": ["O3"]*4 + ["O4"]*4,
        "provenance": "Pathak et al. 2020 combined experimental refinement; "
                      "idealized Hf occupancy; not DFT-relaxed. DOI:10.1107/S2053229620013960",
    }


def normalized_reference(reference):
    required = ("cell_A", "hf_fractional", "oxygen_fractional", "oxygen_kinds")
    if not isinstance(reference, dict) or any(k not in reference for k in required):
        raise ValueError("reference requires cell_A, hf_fractional, oxygen_fractional, oxygen_kinds")
    unit = np.asarray(reference["cell_A"], dtype=float)
    hf = np.asarray(reference["hf_fractional"], dtype=float)
    oxy = np.asarray(reference["oxygen_fractional"], dtype=float)
    kinds = list(reference["oxygen_kinds"])
    if unit.shape != (3, 3) or not np.isfinite(unit).all() or np.linalg.det(unit) <= 1e-8:
        raise ValueError("cell_A must be a finite right-handed nonsingular 3x3 matrix")
    for name, array in (("hf_fractional", hf), ("oxygen_fractional", oxy)):
        if array.ndim != 2 or array.shape[1] != 3 or not len(array) or not np.isfinite(array).all():
            raise ValueError(name + " must be a finite nonempty Nx3 array")
    if len(kinds) != len(oxy) or any(k not in ("O3", "O4") for k in kinds):
        raise ValueError("one O3/O4 label is required per oxygen basis site")
    hf, oxy = np.mod(hf, 1.), np.mod(oxy, 1.)
    if len({tuple(np.round(row, 10)) for row in oxy}) != len(oxy):
        raise ValueError("duplicate oxygen basis site")
    return dict(cell_A=unit.tolist(), hf_fractional=hf.tolist(),
                oxygen_fractional=oxy.tolist(), oxygen_kinds=kinds,
                provenance=str(reference.get("provenance", "unspecified reference")))


def reference_hash(reference):
    ref = normalized_reference(reference)
    identity = {k: ref[k] for k in ("cell_A", "hf_fractional", "oxygen_fractional", "oxygen_kinds")}
    return hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def channel_id(a, b, shift):
    return f"b{int(a)}-b{int(b)}@{int(shift[0])},{int(shift[1])},{int(shift[2])}"


@dataclass(frozen=True)
class HopEdge:
    source: int
    target: int
    channel_id: str
    displacement_A: tuple
    forward: bool = True


class Lattice:
    def __init__(self, repeats=(3, 3, 3), boundary="bulk", migration_cutoff_A=3.2,
                 interaction_cutoff_A=3.6, reference=None):
        if len(repeats) != 3 or any(isinstance(n, bool) or int(n) != n or n < 2 for n in repeats):
            raise ValueError("repeats must be three integers >=2")
        if boundary not in ("bulk", "slab"):
            raise ValueError("boundary must be bulk or slab")
        if any(not math.isfinite(x) or x <= 0 for x in (migration_cutoff_A, interaction_cutoff_A)):
            raise ValueError("graph cutoffs must be finite and positive")
        self.repeats = tuple(int(n) for n in repeats)
        self.boundary = boundary
        self.pbc = (True, True, boundary == "bulk")
        self.reference = normalized_reference(default_reference() if reference is None else reference)
        self.structure_id = reference_hash(self.reference)
        self.unit_A = np.asarray(self.reference["cell_A"])
        self.cell_A = self.unit_A * np.asarray(self.repeats)[:, None]
        self.migration_cutoff_A = float(migration_cutoff_A)
        self.interaction_cutoff_A = float(interaction_cutoff_A)
        # A sufficient (conservative) uniqueness bound, also valid for skew cells:
        # any displacement below s_min(cell)/2 has a unique rounded fractional image.
        bound = np.linalg.svd(self.cell_A, compute_uv=False).min()/2
        if max(migration_cutoff_A, interaction_cutoff_A) >= bound-1e-10:
            raise ValueError("supercell too small for unique periodic images; increase repeats or reduce cutoffs")
        oxy = np.asarray(self.reference["oxygen_fractional"])
        hf = np.asarray(self.reference["hf_fractional"])
        cells = np.asarray(list(product(*[range(n) for n in self.repeats])), dtype=int)
        self.basis_size = len(oxy)
        self.basis_id = np.tile(np.arange(len(oxy)), len(cells))
        self.cell_index = np.repeat(cells, len(oxy), axis=0)
        self.oxygen_A = (self.cell_index + oxy[self.basis_id]) @ self.unit_A
        self.hf_A = np.concatenate([(cell+hf)@self.unit_A for cell in cells])
        self.kind = np.asarray(self.reference["oxygen_kinds"])[self.basis_id]
        self.n = len(self.oxygen_A)
        self._templates = self._make_templates()
        self.migration_edges = self._instantiate(cells)
        self.interaction_pairs = self._pairs(interaction_cutoff_A)
        self.interaction_neighbors = [[] for _ in range(self.n)]
        for i, j in self.interaction_pairs:
            self.interaction_neighbors[i].append(int(j))
            self.interaction_neighbors[j].append(int(i))

    def _make_templates(self):
        basis = np.asarray(self.reference["oxygen_fractional"])
        kinds = self.reference["oxygen_kinds"]
        radius = int(np.ceil(self.migration_cutoff_A/np.linalg.svd(self.unit_A, compute_uv=False).min()))+1
        result = []
        for a, b in product(range(self.basis_size), repeat=2):
            for shift in product(range(-radius, radius+1), repeat=3):
                key, reverse = (a, b, *shift), (b, a, *(-x for x in shift))
                if key > reverse:
                    continue
                displacement = (basis[b]+shift-basis[a])@self.unit_A
                length = float(np.linalg.norm(displacement))
                if not 1e-8 < length <= self.migration_cutoff_A+1e-10:
                    continue
                result.append(dict(id=channel_id(a, b, shift), from_basis=a, to_basis=b,
                                   cell_shift=list(shift), displacement_A=displacement.tolist(),
                                   length_A=length, from_kind=kinds[a], to_kind=kinds[b]))
        return sorted(result, key=lambda record: record["id"])

    def _instantiate(self, cells):
        edges = []
        for cell in cells:
            cell_flat = np.ravel_multi_index(tuple(cell), self.repeats)
            for template in self._templates:
                target_cell = cell+np.asarray(template["cell_shift"])
                valid = True
                for axis in range(3):
                    if self.pbc[axis]:
                        target_cell[axis] %= self.repeats[axis]
                    elif not 0 <= target_cell[axis] < self.repeats[axis]:
                        valid = False
                if not valid:
                    continue
                i = int(cell_flat*self.basis_size+template["from_basis"])
                j = int(np.ravel_multi_index(tuple(target_cell), self.repeats)*self.basis_size+template["to_basis"])
                if i == j:
                    raise ValueError("periodic self hop: enlarge the supercell")
                edges.append(HopEdge(i, j, template["id"], tuple(template["displacement_A"])))
        return edges

    def _pairs(self, cutoff):
        inv = np.linalg.inv(self.cell_A)
        pairs = []
        # O(N) temporary storage, rather than the legacy N x N x 3 allocation.
        for i in range(self.n-1):
            fractional = (self.oxygen_A[i+1:]-self.oxygen_A[i])@inv
            for axis, periodic in enumerate(self.pbc):
                if periodic:
                    fractional[:, axis] -= np.rint(fractional[:, axis])
            distances = np.linalg.norm(fractional@self.cell_A, axis=1)
            pairs.extend((i, int(j+i+1)) for j in np.flatnonzero((distances <= cutoff+1e-10)&(distances > 1e-8)))
        return np.asarray(pairs, dtype=int).reshape(-1, 2)

    def channel_templates(self):
        return json.loads(json.dumps(self._templates))

    def to_reference_json(self):
        return json.loads(json.dumps(self.reference))
