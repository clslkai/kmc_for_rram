"""Validated, provenance-bearing input for classical fixed-charge lattice KMC.

Version 1 is scoped to a conventional hafnia basis of four Hf and eight O.

``saddle_eV`` is an absolute transition-state energy on the *same reference*
as the eight vacancy-site energies for its charge branch.  It is neither a
forward activation barrier nor a formation energy to be added to that barrier.
The v1 schema has ONE ``attempt_Hz`` per charge/channel, shared by forward and
reverse jumps.  It therefore describes classical energies with no vibrational
entropy.  Unequal directional HTST prefactors require a future free-energy
schema; they must not be inserted as silently ignored additional fields.

The factory generates either explicitly synthetic examples or incomplete
worklists.  It never generates DFT results.  ``audit`` checks completeness and
internal consistency, not the scientific authenticity of supplied provenance.
"""
from __future__ import annotations

from copy import deepcopy
import json
import math
from numbers import Integral, Real
from pathlib import Path
import re
import tempfile
from typing import Any

import numpy as np

SCHEMA_VERSION = "dftkmc-1"
PAIR_KEYS = ("O3-O3", "O3-O4", "O4-O4")
STATUSES = ("illustrative", "pending", "dft")
RATE_CONVENTION = "classical_shared_saddle_equal_prefactor"
RATE_FIELDS = {"saddle_eV", "attempt_Hz", "ts_fraction", "status", "provenance"}
_CHARGE_RE = re.compile(r"(?:0|-[1-9][0-9]*|[1-9][0-9]*)\Z")


def _finite(value: Any) -> bool:
    if not isinstance(value, Real) or isinstance(value, (bool, np.bool_)):
        return False
    try:
        return math.isfinite(value)
    except (ValueError, OverflowError):
        return False


def _charge_key(q: Any) -> str:
    if not isinstance(q, Integral) or isinstance(q, (bool, np.bool_)):
        raise ValueError("charge must be an integer, not a fractional value or boolean")
    return str(int(q))


def _nonempty(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _json_pairs(pairs: list[tuple[str, Any]]) -> dict:
    result: dict = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number is not allowed: {value}")


def _reference_lattice(reference: dict, migration_cutoff_A: float, pair_cutoff_A: float):
    """Use the geometry owner's canonical hash/IDs, never a second convention."""
    from .geometry import Lattice
    singular_min = float(np.linalg.svd(np.asarray(reference["cell_A"], dtype=float), compute_uv=False).min())
    repeat = max(2, int(math.floor(2 * max(migration_cutoff_A, pair_cutoff_A) / singular_min)) + 1)
    return Lattice(
        repeats=(repeat, repeat, repeat), boundary="bulk", reference=reference,
        migration_cutoff_A=migration_cutoff_A,
        interaction_cutoff_A=pair_cutoff_A,
    )


class Catalog:
    """A catalogue which can also hold incomplete, non-runnable DFT worklists.

    Construction deliberately does not imply validation.  ``load`` validates by
    default, while ``load(validate=False)`` enables inspection of pending data.
    Call ``validate`` before constructing a rate model.  The public ``data`` is
    editable, so validation is recalculated rather than cached.
    """

    def __init__(self, data: dict):
        if not isinstance(data, dict):
            raise ValueError("catalogue root must be a dictionary")
        self.data = deepcopy(data)

    @classmethod
    def load(cls, path, allow_illustrative: bool = False, validate: bool = True) -> "Catalog":
        with Path(path).open(encoding="utf-8") as handle:
            data = json.load(handle, object_pairs_hook=_json_pairs, parse_constant=_json_constant)
        catalog = cls(data)
        if validate:
            catalog.validate(allow_illustrative=allow_illustrative)
        return catalog

    @property
    def charges(self) -> tuple[int, ...]:
        energies = self.data.get("site_energies_eV")
        if not isinstance(energies, dict) or not energies:
            raise ValueError("site_energies_eV must contain at least one charge branch")
        if any(not isinstance(key, str) or not _CHARGE_RE.fullmatch(key) for key in energies):
            raise ValueError("charge keys must be canonical integer strings, e.g. '0', '2', '-1'")
        return tuple(sorted(int(key) for key in energies))

    @property
    def structure(self) -> dict:
        value = self.data.get("structure")
        if not isinstance(value, dict):
            raise ValueError("catalogue structure is missing or invalid")
        return value

    @property
    def structure_id(self) -> str:
        metadata = self.data.get("metadata")
        value = metadata.get("structure_id") if isinstance(metadata, dict) else None
        if not _nonempty(value):
            raise ValueError("metadata.structure_id is missing")
        return value

    def site_energies(self, q: int) -> np.ndarray:
        key = _charge_key(q)
        mapping = self.data.get("site_energies_eV", {})
        if not isinstance(mapping, dict) or key not in mapping:
            raise ValueError(f"site energies do not cover charge {key}")
        values = mapping[key]
        if not isinstance(values, list) or len(values) != 8 or not all(_finite(v) for v in values):
            raise ValueError(f"site energies for charge {key} are incomplete or non-finite")
        return np.asarray(values, dtype=float).copy()

    def charge_model(self, q: int) -> dict:
        key = _charge_key(q)
        mapping = self.data.get("charge_models", {})
        if not isinstance(mapping, dict) or not isinstance(mapping.get(key), dict):
            raise ValueError(f"charge model does not cover charge {key}")
        return mapping[key]

    def template(self, channel_id: str) -> dict:
        channels = self.data.get("channels")
        if not isinstance(channels, list):
            raise ValueError("channels must be a list")
        matches = [channel for channel in channels if isinstance(channel, dict) and channel.get("id") == channel_id]
        if len(matches) != 1:
            reason = "not covered" if not matches else "ambiguous (duplicate ID)"
            raise ValueError(f"channel {channel_id!r} is {reason}; no fallback barrier is permitted")
        return matches[0]

    def channel(self, channel_id: str, q: int) -> dict:
        key = _charge_key(q)
        mapping = self.template(channel_id).get("data")
        if not isinstance(mapping, dict) or not isinstance(mapping.get(key), dict):
            raise ValueError(f"channel {channel_id!r} does not cover charge {key}; no fallback barrier is permitted")
        return mapping[key]

    def audit(self) -> dict:
        """Report all detectable missing/inconsistent input without filling it.

        ``ready`` means numerically complete and internally consistent.  An
        illustrative catalogue still needs explicit permission at ``validate``.
        A true report does not certify that calculations were performed or that
        a geometry is physically appropriate for a particular experiment.
        """
        errors: list[str] = []
        warnings: list[str] = []
        missing: list[str] = []

        def error(message: str) -> None:
            errors.append(message)

        def number(value: Any, path: str, *, positive=False, nonnegative=False) -> bool:
            if value is None:
                missing.append(path)
                error(f"{path} is null or missing (pending data)")
                return False
            if not _finite(value):
                error(f"{path} must be a finite number, not a boolean")
                return False
            if positive and value <= 0:
                error(f"{path} must be positive")
                return False
            if nonnegative and value < 0:
                error(f"{path} must be nonnegative")
                return False
            return True

        def provenance(value: Any, path: str) -> None:
            if not _nonempty(value) or value.strip().lower() in {"unspecified", "unspecified reference", "unknown", "n/a"}:
                error(f"{path} must explicitly identify the source or pending/synthetic nature")

        def charge_keys(value: Any, path: str) -> set[str]:
            if not isinstance(value, dict):
                error(f"{path} must be a charge-keyed dictionary")
                return set()
            valid = set()
            for key in value:
                if not isinstance(key, str) or not _CHARGE_RE.fullmatch(key):
                    error(f"{path} contains noncanonical/integral charge key {key!r}")
                else:
                    valid.add(key)
            return valid

        if self.data.get("schema_version") != SCHEMA_VERSION:
            error(f"schema_version must be {SCHEMA_VERSION!r}")
        metadata = self.data.get("metadata")
        if not isinstance(metadata, dict):
            error("metadata must be a dictionary")
            metadata = {}
        status = metadata.get("calibration_status")
        if status not in STATUSES:
            error(f"metadata.calibration_status must be one of {STATUSES}")
        if not _nonempty(metadata.get("material")):
            error("metadata.material must be nonempty")
        for key, expected in (("energy_unit", "eV"), ("length_unit", "angstrom")):
            if metadata.get(key) != expected:
                error(f"metadata.{key} must be {expected!r}")
        provenance(metadata.get("provenance"), "metadata.provenance")
        exclusions = metadata.get("channel_exclusions", {})
        excluded_ids: set[str] = set()
        if not isinstance(exclusions, dict):
            error("metadata.channel_exclusions must map canonical IDs to explicit truncation reasons")
        else:
            for excluded_id, reason in exclusions.items():
                if not _nonempty(excluded_id):
                    error("metadata.channel_exclusions keys must be nonempty canonical IDs")
                    continue
                excluded_ids.add(excluded_id)
                if not _nonempty(reason):
                    error(f"channel exclusion {excluded_id!r} needs a nonempty scientific/model-truncation reason")
        if excluded_ids:
            warnings.append("Explicitly truncated channel graph: dynamics are conditional on the declared exclusions and cannot establish full diffusion in the material.")
        if metadata.get("rate_convention", RATE_CONVENTION) != RATE_CONVENTION:
            error("only classical shared-saddle energies with equal forward/reverse prefactors are supported")
        if metadata.get("vibrational_entropy_included", False) is not False:
            error("vibrational entropy is unsupported by dftkmc-1; endpoint energies must be classical")
        if status == "dft" and metadata.get("synthetic_data") is True:
            error("synthetic_data cannot be promoted to DFT by relabelling the calibration status")
        if status == "illustrative":
            warnings.append("Synthetic, non-literature demonstration parameters: not DFT-calibrated predictions; explicit allow_illustrative is required.")
        elif status == "pending":
            error("calibration_status is pending; this catalogue cannot run")
        warnings.append("Classical energies only; one attempt_Hz is shared by both directions of each charge/channel. No vibrational entropy is included.")

        reference = self.data.get("structure")
        valid_reference = isinstance(reference, dict)
        arrays: dict[str, np.ndarray] = {}
        if not valid_reference:
            error("structure must be a reference dictionary")
            reference = {}
        provenance(reference.get("provenance"), "structure.provenance")
        for key, shape in (("cell_A", (3, 3)), ("hf_fractional", (4, 3)), ("oxygen_fractional", (8, 3))):
            value = reference.get(key)
            try:
                array = np.asarray(value, dtype=object)
                ok = array.shape == shape and all(_finite(v) for v in array.flat)
            except (TypeError, ValueError, OverflowError):
                ok = False
            if not ok:
                error(f"structure.{key} must be a finite array of shape {shape}")
                valid_reference = False
            else:
                arrays[key] = np.asarray(array, dtype=float)
        if "cell_A" in arrays and abs(float(np.linalg.det(arrays["cell_A"]))) <= 1e-12:
            error("structure.cell_A must be nonsingular")
            valid_reference = False
        kinds = reference.get("oxygen_kinds")
        if not isinstance(kinds, list) or len(kinds) != 8 or any(v not in ("O3", "O4") for v in kinds):
            error("structure.oxygen_kinds must contain eight explicit O3/O4 labels")
            valid_reference = False

        sites = self.data.get("site_energies_eV")
        charges = charge_keys(sites, "site_energies_eV")
        if not charges:
            error("site_energies_eV must cover at least one integral charge")
        site_values: dict[str, list] = {}
        if isinstance(sites, dict):
            for key in sorted(charges, key=int):
                values = sites[key]
                if not isinstance(values, list) or len(values) != 8:
                    error(f"site_energies_eV[{key}] must contain exactly eight site energies")
                    continue
                site_values[key] = values
                for basis, value in enumerate(values):
                    number(value, f"site_energies_eV[{key}][{basis}]")

        models = self.data.get("charge_models")
        model_charges = charge_keys(models, "charge_models")
        if charges != model_charges:
            error(f"charge_models charge coverage differs from site energies: missing={sorted(charges-model_charges)}, extra={sorted(model_charges-charges)}")
        pair_cutoffs = []
        if isinstance(models, dict):
            for key in sorted(model_charges, key=int):
                path = f"charge_models[{key}]"
                model = models[key]
                if not isinstance(model, dict):
                    error(f"{path} must be a dictionary")
                    continue
                if number(model.get("pair_cutoff_A"), f"{path}.pair_cutoff_A", positive=True):
                    pair_cutoffs.append(float(model["pair_cutoff_A"]))
                provenance(model.get("provenance"), f"{path}.provenance")
                mode = model.get("environment_mode")
                if mode not in ("dilute_only", "symmetric_pair"):
                    error(f"{path}.environment_mode must be dilute_only or symmetric_pair")
                pair = model.get("pair_energies_eV")
                if not isinstance(pair, dict):
                    error(f"{path}.pair_energies_eV must be a dictionary with all three pair types")
                else:
                    if set(pair) != set(PAIR_KEYS):
                        error(f"{path}.pair_energies_eV must contain exactly {PAIR_KEYS}")
                    for pair_key in PAIR_KEYS:
                        number(pair.get(pair_key), f"{path}.pair_energies_eV[{pair_key}]")
                gamma_ok = number(model.get("environment_saddle_eV_per_neighbor"), f"{path}.environment_saddle_eV_per_neighbor")
                if mode == "dilute_only" and gamma_ok and model["environment_saddle_eV_per_neighbor"] != 0:
                    error(f"{path}: dilute_only requires zero environment-saddle correction")
                if mode == "symmetric_pair":
                    warnings.append(f"charge {key}: symmetric_pair is a declared model approximation, not environment-resolved NEB data.")

        channels = self.data.get("channels")
        if not isinstance(channels, list) or not channels:
            error("channels must be a nonempty list; missing channels have no fallback")
            channels = []
        seen_ids = set()
        seen_paths = set()
        observed_templates: list[tuple[dict, tuple, str]] = []
        lengths = []
        status_counts = {name: 0 for name in STATUSES}
        present_by_charge = {key: 0 for key in sorted(charges, key=int)}
        for index, channel in enumerate(channels):
            path = f"channels[{index}]"
            if not isinstance(channel, dict):
                error(f"{path} must be a dictionary")
                continue
            channel_id = channel.get("id")
            if not _nonempty(channel_id):
                error(f"{path}.id must be nonempty")
            elif channel_id in seen_ids:
                error(f"duplicate channel ID {channel_id!r}")
            else:
                seen_ids.add(channel_id)
            endpoints = [channel.get("from_basis"), channel.get("to_basis")]
            endpoint_ok = all(isinstance(v, Integral) and not isinstance(v, (bool, np.bool_)) and 0 <= v < 8 for v in endpoints)
            if not endpoint_ok:
                error(f"{path} endpoint basis IDs must be integers in [0, 7]")
            shift = channel.get("cell_shift")
            shift_ok = isinstance(shift, (list, tuple)) and len(shift) == 3 and all(isinstance(v, Integral) and not isinstance(v, (bool, np.bool_)) for v in shift)
            if not shift_ok:
                error(f"{path}.cell_shift must contain three integers")
            displacement = channel.get("displacement_A")
            displacement_ok = isinstance(displacement, (list, tuple)) and len(displacement) == 3 and all(_finite(v) for v in displacement)
            if not displacement_ok:
                error(f"{path}.displacement_A must contain three finite values")
            length_ok = number(channel.get("length_A"), f"{path}.length_A", positive=True)
            if length_ok:
                lengths.append(float(channel["length_A"]))
            if endpoint_ok and shift_ok:
                a, b = (int(v) for v in endpoints)
                direct = (a, b, *[int(v) for v in shift])
                reverse = (b, a, *[-int(v) for v in shift])
                canonical = min(direct, reverse)
                if direct == reverse:
                    error(f"{path} is a zero-displacement self-channel")
                if canonical in seen_paths:
                    error(f"{path} duplicates an existing channel or its reverse")
                seen_paths.add(canonical)
                observed_templates.append((channel, direct, path))
                if valid_reference:
                    expected_d = (arrays["oxygen_fractional"][b] + np.asarray(shift, dtype=float) - arrays["oxygen_fractional"][a]) @ arrays["cell_A"]
                    if displacement_ok and not np.allclose(displacement, expected_d, rtol=1e-9, atol=1e-8):
                        error(f"{path}.displacement_A is inconsistent with the reference and cell_shift")
                    if length_ok and not np.isclose(channel["length_A"], np.linalg.norm(expected_d), rtol=1e-9, atol=1e-8):
                        error(f"{path}.length_A is inconsistent with the reference displacement")
                    if channel.get("from_kind") != kinds[a] or channel.get("to_kind") != kinds[b]:
                        error(f"{path} O3/O4 endpoint labels disagree with the explicit structure")
            data = channel.get("data")
            data_charges = charge_keys(data, f"{path}.data")
            if data_charges != charges:
                error(f"{path}.data charge coverage differs from site energies: missing={sorted(charges-data_charges)}, extra={sorted(data_charges-charges)}")
            if not isinstance(data, dict):
                continue
            for key in sorted(data_charges, key=int):
                rate = data[key]
                rpath = f"{path}.data[{key}]"
                if not isinstance(rate, dict):
                    error(f"{rpath} must be a dictionary")
                    continue
                unexpected = set(rate) - RATE_FIELDS
                if unexpected:
                    error(f"{rpath} has unsupported fields {sorted(map(str, unexpected))}; use one shared attempt_Hz, never directional prefactors")
                provenance(rate.get("provenance"), f"{rpath}.provenance")
                rate_status = rate.get("status")
                if key in present_by_charge:
                    present_by_charge[key] += 1
                if rate_status in STATUSES:
                    status_counts[rate_status] += 1
                if rate_status not in STATUSES:
                    error(f"{rpath}.status must be one of {STATUSES}")
                elif rate_status != status:
                    error(f"{rpath}.status disagrees with metadata.calibration_status")
                if rate_status == "pending":
                    error(f"{rpath} is pending; no rate can be constructed")
                saddle_ok = number(rate.get("saddle_eV"), f"{rpath}.saddle_eV")
                number(rate.get("attempt_Hz"), f"{rpath}.attempt_Hz", positive=True)
                if number(rate.get("ts_fraction"), f"{rpath}.ts_fraction") and not 0 <= rate["ts_fraction"] <= 1:
                    error(f"{rpath}.ts_fraction must lie in [0, 1]")
                if saddle_ok and endpoint_ok and key in site_values:
                    endpoint_energies = [site_values[key][int(v)] for v in endpoints]
                    if all(_finite(v) for v in endpoint_energies) and rate["saddle_eV"] < max(endpoint_energies) - 1e-12:
                        error(f"{rpath}.saddle_eV is below an endpoint energy; saddle is an absolute shared energy, not a barrier")

        # Reconstruct from the declared reference using the single authoritative
        # geometry implementation, so neither IDs nor hashes can drift silently.
        cutoff = metadata.get("migration_cutoff_A")
        if excluded_ids and cutoff is None:
            error("metadata.migration_cutoff_A must be explicit when channel_exclusions are declared")
        if cutoff is None and lengths:
            cutoff = max(lengths) + 1e-8
            warnings.append("migration_cutoff_A is not declared; geometric coverage is checked through the longest supplied channel only. Engine channels outside this range must fail closed.")
        cutoff_ok = _finite(cutoff) and cutoff > 0
        if not cutoff_ok:
            error("metadata.migration_cutoff_A must be positive and finite, or derivable from valid channel lengths")
        if valid_reference and cutoff_ok:
            try:
                lattice = _reference_lattice(reference, float(cutoff), max(pair_cutoffs, default=3.6))
                if metadata.get("structure_id") != lattice.structure_id:
                    error("metadata.structure_id does not match the declared cell, basis and oxygen kinds")
                expected = {item["id"]: item for item in lattice.channel_templates()}
                expected_paths = {
                    (item["from_basis"], item["to_basis"], *item["cell_shift"]): item["id"]
                    for item in expected.values()
                }
                for channel, direct, path in observed_templates:
                    channel_id = channel.get("id")
                    if direct not in expected_paths:
                        error(f"{path} is not a canonical geometry channel within the declared cutoff")
                    elif channel_id != expected_paths[direct]:
                        error(f"{path}.id disagrees with its canonical geometry channel")
                overlap = sorted(excluded_ids & seen_ids)
                if overlap:
                    error(f"channel exclusions overlap actual channels: {', '.join(overlap)}")
                unknown_exclusions = sorted(excluded_ids - set(expected))
                if unknown_exclusions:
                    error(f"unknown excluded geometry channels: {', '.join(unknown_exclusions)}")
                missing_ids = sorted(set(expected) - seen_ids - excluded_ids)
                extra_ids = sorted(seen_ids - set(expected))
                if missing_ids:
                    error(f"uncovered geometry channels: {', '.join(missing_ids)}")
                if extra_ids:
                    error(f"unknown geometry channels: {', '.join(extra_ids)}")
            except (ValueError, TypeError, KeyError, OverflowError) as exc:
                error(f"cannot verify reference/channel geometry: {exc}")

        return {
            "ready": not errors,
            "errors": errors,
            "warnings": warnings,
            "missing": missing,
            "schema_version": self.data.get("schema_version"),
            "calibration_status": status,
            "charges": sorted((int(key) for key in charges)),
            "n_channels": len(channels),
            "n_excluded_channels": len(excluded_ids),
            "declared_dft_calibrated": not errors and status == "dft",
            "calibration_verification": "Declaration only: validation checks input consistency and cannot certify that DFT/NEB calculations were performed or scientifically validated.",
            "channel_status_counts": status_counts,
            "present_channel_entries_by_charge": present_by_charge,
            "n_missing_values": len(missing),
            "requires_allow_illustrative": status == "illustrative",
            "rate_convention": RATE_CONVENTION,
        }

    def validate(self, allow_illustrative: bool = False) -> "Catalog":
        report = self.audit()
        errors = list(report["errors"])
        if report["calibration_status"] == "illustrative" and not allow_illustrative:
            errors.insert(0, "illustrative catalogue requires explicit allow_illustrative=True; these are synthetic, non-literature parameters")
        if errors:
            shown = errors[:15]
            suffix = f"\n... {len(errors)-len(shown)} further issues; call audit() for the complete report" if len(errors) > len(shown) else ""
            raise ValueError("Catalogue is not runnable:\n- " + "\n- ".join(shown) + suffix)
        return self

    def write(self, path) -> None:
        """Write strict JSON atomically, including intentionally pending nulls."""
        destination = Path(path)
        payload = json.dumps(self.data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
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


def make_catalog(lattice, status: str = "illustrative", charges=(0, 2)) -> Catalog:
    """Create a clearly synthetic example or a null-valued DFT worklist.

    ``status='dft'`` is deliberately unavailable: real energies, charge choices,
    saddle positions and provenance must come from actual supplied calculations.
    For pending data the midpoint ``ts_fraction`` is only a marked placeholder.
    ``symmetric_pair`` in the example is an explicitly declared approximation.
    Optional ``metadata.channel_exclusions`` maps intentionally omitted canonical
    paths to explicit truncation reasons.  The default catalogue excludes none;
    excluded and supplied paths must be disjoint and cover the declared cutoff.
    """
    if status not in ("illustrative", "pending"):
        raise ValueError("make_catalog only creates illustrative or pending data; it cannot manufacture DFT results")
    keys = [_charge_key(q) for q in charges]
    if not keys or len(set(keys)) != len(keys):
        raise ValueError("charges must be nonempty and unique")
    keys.sort(key=int)
    reference = deepcopy(lattice.to_reference_json())
    if (not isinstance(reference, dict) or not _nonempty(reference.get("provenance"))
            or reference["provenance"].strip().lower() in {"unspecified", "unspecified reference", "unknown", "n/a"}):
        raise ValueError("the lattice must supply explicit structure provenance; it cannot be invented by make_catalog")
    structure_id = lattice.structure_id
    if not _nonempty(structure_id):
        raise ValueError("the lattice must supply its structure_id")
    templates = deepcopy(lattice.channel_templates())
    if not isinstance(templates, list) or not templates:
        raise ValueError("the lattice must provide nonempty canonical channel templates")
    pair_cutoff = getattr(lattice, "interaction_cutoff_A", None)
    if not _finite(pair_cutoff) or pair_cutoff <= 0:
        raise ValueError("the lattice must declare positive interaction_cutoff_A; it cannot be inferred from pair energies")
    cutoff = getattr(lattice, "migration_cutoff_A", None)
    if cutoff is None:
        cutoff = max(item["length_A"] for item in templates) + 1e-8
    if not _finite(cutoff) or cutoff <= 0:
        raise ValueError("the lattice must provide a finite positive migration cutoff")
    kinds = reference.get("oxygen_kinds")
    if not isinstance(kinds, list) or len(kinds) != 8 or any(k not in ("O3", "O4") for k in kinds):
        raise ValueError("the lattice must supply eight explicit O3/O4 basis labels")
    illustrative = status == "illustrative"
    origin = (
        "SYNTHETIC NON-LITERATURE DEMONSTRATION. Arbitrary parameters generated solely to exercise software and detailed balance; no DFT calculations or material predictions."
        if illustrative else
        "PENDING DFT/NEB WORKLIST. No computed energies or attempt frequencies have been supplied. Charges are requested branches, not established stable charge states."
    )
    data = {
        "schema_version": SCHEMA_VERSION,
        "metadata": {
            "material": "HfO2", "structure_id": structure_id,
            "calibration_status": status, "energy_unit": "eV", "length_unit": "angstrom",
            "provenance": origin, "synthetic_data": illustrative,
            "migration_cutoff_A": float(cutoff), "channel_exclusions": {},
            "rate_convention": RATE_CONVENTION, "vibrational_entropy_included": False,
            "energy_reference": "Common arbitrary classical-energy reference within each fixed charge branch; saddle_eV is absolute on that reference. Different q branches are not thermodynamically compared.",
            "prefactor_convention": "One equal forward/reverse attempt_Hz per charge/channel; unequal HTST prefactors are unsupported.",
            "charge_selection": "Requested independent fixed-charge branches; no charge-state population or electron capture model.",
            "channel_reduction": "Translation and reverse equivalence only; no full space-group reduction is claimed.",
        },
        "structure": reference,
        "site_energies_eV": {}, "charge_models": {}, "channels": templates,
    }
    for key in keys:
        data["site_energies_eV"][key] = [0.0 if kind == "O3" else 0.10 for kind in kinds] if illustrative else [None] * 8
        data["charge_models"][key] = {
            "pair_cutoff_A": float(pair_cutoff),
            "pair_energies_eV": dict(zip(PAIR_KEYS, (-0.05, -0.04, -0.03))) if illustrative else dict.fromkeys(PAIR_KEYS),
            "environment_mode": "symmetric_pair" if illustrative else "dilute_only",
            "environment_saddle_eV_per_neighbor": 0.015 if illustrative else 0.0,
            "provenance": origin + (" Symmetric pair and saddle-neighbor terms are declared illustrative approximations, not environment-resolved NEB." if illustrative else " Dilute-only eligibility; environment saddle correction is disabled. Pair energies remain unknown, not assumed zero."),
        }
    for template in templates:
        if "data" in template:
            raise ValueError("lattice channel templates must not already contain calibration data")
        template["data"] = {}
        for key in keys:
            energies = data["site_energies_eV"][key]
            saddle = None
            if illustrative:
                saddle = max(energies[template["from_basis"]], energies[template["to_basis"]]) + 0.80 + 0.015 * abs(int(key)) + 0.02 * float(template["length_A"])
            template["data"][key] = {
                "saddle_eV": saddle, "attempt_Hz": 1e13 if illustrative else None,
                "ts_fraction": 0.5, "status": status,
                "provenance": origin + (" Arbitrary saddle=max(endpoint energies)+0.80+0.015*abs(q)+0.02*length_A eV; midpoint saddle and equal 1e13 Hz prefactors are illustrative." if illustrative else " ts_fraction=0.5 is an uncalibrated midpoint placeholder; replace with a documented transition-state projection before DFT calibration."),
            }
    return Catalog(data)
