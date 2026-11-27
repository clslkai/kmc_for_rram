"""Synthetic raw-data fixtures and a strict, auditable import path.

This module does not run, parse, or certify DFT/NEB calculations.  Its deliberately
synthetic ``dftkmc-mock-1`` JSON bundle mimics records one could receive from a
calculation workflow.  The complete geometry is covered; no missing barrier is
inferred, interpolated, or replaced by a fallback value.

Bundle schema (all numbers in the units explicitly declared in metadata)::

    {
      "schema_version": "dftkmc-mock-1",
      "metadata": {
        "material": "HfO2", "charges": [0, 2], "structure_id": "<sha256>",
        "calibration_status": "illustrative", "synthetic_data": true,
        "energy_unit": "eV", "length_unit": "angstrom", "frequency_unit": "Hz",
        "migration_cutoff_A": 3.2, "provenance": "SYNTHETIC ...",
        "energy_reference": "<common reference within each fixed-q branch>"
      },
      "structure": {"cell_A": [[...]], "hf_fractional": [[...]],
                    "oxygen_fractional": [[...]], "oxygen_kinds": [...],
                    "provenance": "<actual source of the geometric reference>"},
      "site_energies": {"0": [
        {"basis_id": 0, "kind": "O3", "energy_eV": 0.0,
         "status": "illustrative", "synthetic_data": true,
         "provenance": "SYNTHETIC ..."}, "<seven further records>"
      ], "2": ["<eight records>"]},
      "charge_models": {"0": {
        "pair_cutoff_A": 3.6,
        "pair_energies_eV": {"O3-O3": -0.04, "O3-O4": -0.05, "O4-O4": -0.06},
        "environment_mode": "symmetric_pair",
        "environment_saddle_eV_per_neighbor": 0.02,
        "status": "illustrative", "synthetic_data": true,
        "provenance": "SYNTHETIC ..."
      }, "2": {"<same fields>"}},
      "neb_paths": [
        {"<all fields from Lattice.channel_templates()>": "...", "data": {
          "0": {"barrier_forward_eV": 0.8, "barrier_reverse_eV": 0.7,
                "attempt_forward_Hz": 1e12, "attempt_reverse_Hz": 1e12,
                "ts_fraction": 0.5, "status": "illustrative",
                "synthetic_data": true, "provenance": "SYNTHETIC ..."},
          "2": {"<same fields>"}}}
      ]
    }

Optional metadata: bundle_id, generator, rate_convention,
vibrational_entropy_included and charge_selection.  Unsupported fields are
rejected, so a DFT status or a directional prefactor cannot be silently ignored.
The importer computes S = epsilon_initial + B_forward and independently checks
S = epsilon_final + B_reverse (absolute tolerance 1e-9 eV).  It requires exactly
equal forward/reverse prefactors, consistent with the classical-energy model.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from numbers import Integral, Real
from pathlib import Path
import tempfile

import numpy as np

from .catalog import Catalog, PAIR_KEYS, RATE_CONVENTION, SCHEMA_VERSION as CATALOG_SCHEMA
from .geometry import Lattice

SCHEMA_VERSION = "dftkmc-mock-1"
ENERGY_CLOSURE_TOL_EV = 1e-9
SYNTHETIC_PROVENANCE = (
    "SYNTHETIC NON-LITERATURE VIRTUAL DATA. Arbitrary numerical fixtures for "
    "testing import, KMC rates and statistics; no DFT/NEB calculation was performed. "
    "Different q baselines are arbitrary and do not assert a material or literature trend."
)
META_REQUIRED = {
    "material", "charges", "structure_id", "calibration_status", "synthetic_data",
    "energy_unit", "length_unit", "frequency_unit", "migration_cutoff_A",
    "provenance", "energy_reference",
}
META_OPTIONAL = {
    "bundle_id", "generator", "rate_convention", "vibrational_entropy_included", "charge_selection",
}
TEMPLATE_FIELDS = {
    "id", "from_basis", "to_basis", "cell_shift", "displacement_A", "length_A", "from_kind", "to_kind",
}
SITE_FIELDS = {"basis_id", "kind", "energy_eV", "status", "synthetic_data", "provenance"}
MODEL_FIELDS = {
    "pair_cutoff_A", "pair_energies_eV", "environment_mode",
    "environment_saddle_eV_per_neighbor", "status", "synthetic_data", "provenance",
}
BARRIER_FIELDS = {
    "barrier_forward_eV", "barrier_reverse_eV", "attempt_forward_Hz",
    "attempt_reverse_Hz", "ts_fraction", "status", "synthetic_data", "provenance",
}


def _fields(value, required, path, optional=()):
    if not isinstance(value, dict):
        raise ValueError(f"{path} must be an object")
    missing, unknown = set(required) - set(value), set(value) - set(required) - set(optional)
    if missing:
        raise ValueError(f"{path} is missing required fields: {', '.join(sorted(missing))}")
    if unknown:
        raise ValueError(f"{path} has unsupported fields: {', '.join(sorted(map(str, unknown)))}")


def _text(value, path, synthetic=False):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{path} must be a nonempty provenance/description string")
    if value.strip().lower() in {"unknown", "unspecified", "unspecified reference", "n/a"}:
        raise ValueError(f"{path} must identify its source, not an unspecified placeholder")
    if synthetic and "SYNTHETIC" not in value.upper():
        raise ValueError(f"{path} must explicitly contain the SYNTHETIC provenance label")
    return value


def _number(value, path, positive=False, nonnegative=False):
    if value is None:
        raise ValueError(f"{path} is null: missing scientific values cannot be filled during import")
    if not isinstance(value, Real) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{path} must be a finite number, not a boolean or string")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{path} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{path} must be finite")
    if positive and result <= 0:
        raise ValueError(f"{path} must be positive")
    if nonnegative and result < 0:
        raise ValueError(f"{path} must be nonnegative")
    return result


def _index(value, path):
    if not isinstance(value, Integral) or isinstance(value, (bool, np.bool_)) or not 0 <= value < 8:
        raise ValueError(f"{path} must be an integer basis index in [0, 7]")
    return int(value)


def _synthetic(record, path, metadata=False):
    status_key = "calibration_status" if metadata else "status"
    if record.get(status_key) != "illustrative" or record.get("synthetic_data") is not True:
        raise ValueError(f"{path} must declare {status_key}='illustrative' and synthetic_data=true; DFT status is forbidden")
    _text(record.get("provenance"), f"{path}.provenance", synthetic=True)


def _charge_mapping(mapping, keys, path):
    if not isinstance(mapping, dict) or set(mapping) != set(keys):
        raise ValueError(f"{path} must cover exactly the declared canonical charge keys {list(keys)}")


def source_digest(bundle: dict) -> str:
    """SHA256 of the complete canonical JSON content, not an authenticity claim.

    Object keys are sorted, arrays retain their order, and compact UTF-8 JSON is
    hashed without an added newline.  NaN/Infinity and non-JSON objects fail.
    """
    if not isinstance(bundle, dict):
        raise ValueError("mock bundle must be an object")
    try:
        payload = json.dumps(bundle, sort_keys=True, separators=(",", ":"),
                             ensure_ascii=False, allow_nan=False).encode("utf-8")
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"mock bundle must be strict finite JSON: {exc}") from exc
    return hashlib.sha256(payload).hexdigest()


def make_mock_bundle(lattice=None) -> dict:
    """Return deterministic synthetic raw data for every supplied geometry path.

    The default reference has 34 canonical paths and q=0/q=2 branches.  All
    energies and rates are arbitrary.  The small squared z-direction correction
    is invariant under path reversal; it does not represent fitted anisotropy.
    """
    lattice = Lattice(repeats=(2, 2, 2)) if lattice is None else lattice
    reference = deepcopy(lattice.to_reference_json())
    _text(reference.get("provenance"), "structure.provenance")
    kinds = reference.get("oxygen_kinds")
    if not isinstance(kinds, list) or len(kinds) != 8 or any(kind not in ("O3", "O4") for kind in kinds):
        raise ValueError("mock schema requires eight explicitly labelled O3/O4 oxygen basis sites")
    migration_cutoff = _number(lattice.migration_cutoff_A, "lattice.migration_cutoff_A", positive=True)
    pair_cutoff = _number(lattice.interaction_cutoff_A, "lattice.interaction_cutoff_A", positive=True)
    templates = deepcopy(lattice.channel_templates())
    if not isinstance(templates, list) or not templates:
        raise ValueError("lattice must provide nonempty canonical channel templates")
    bundle = {
        "schema_version": SCHEMA_VERSION,
        "metadata": {
            "bundle_id": "hafnia-virtual-fixed-charge-v1", "generator": "dftkmc.mockdata.make_mock_bundle:v1",
            "material": "HfO2", "charges": [0, 2], "structure_id": lattice.structure_id,
            "calibration_status": "illustrative", "synthetic_data": True,
            "energy_unit": "eV", "length_unit": "angstrom", "frequency_unit": "Hz",
            "migration_cutoff_A": migration_cutoff, "provenance": SYNTHETIC_PROVENANCE,
            "energy_reference": "One common arbitrary classical-energy zero per fixed-q branch. Absolute energies from different q branches must not be compared.",
            "rate_convention": RATE_CONVENTION, "vibrational_entropy_included": False,
            "charge_selection": "q=0 and q=2 are requested independent synthetic branches, not predicted charge-state populations.",
        },
        "structure": reference, "site_energies": {}, "charge_models": {}, "neb_paths": templates,
    }
    for q, o4_energy in ((0, 0.10), (2, 0.06)):
        key = str(q)
        bundle["site_energies"][key] = [
            {"basis_id": basis, "kind": kind, "energy_eV": 0.0 if kind == "O3" else o4_energy,
             "status": "illustrative", "synthetic_data": True,
             "provenance": SYNTHETIC_PROVENANCE + f" Site recipe for q={q}: O3=0 eV and O4={o4_energy} eV."}
            for basis, kind in enumerate(kinds)
        ]
        bundle["charge_models"][key] = {
            "pair_cutoff_A": pair_cutoff,
            "pair_energies_eV": dict(zip(PAIR_KEYS, (-0.04, -0.05, -0.06))),
            "environment_mode": "symmetric_pair", "environment_saddle_eV_per_neighbor": 0.02,
            "status": "illustrative", "synthetic_data": True,
            "provenance": SYNTHETIC_PROVENANCE + " Pair interactions (-0.04,-0.05,-0.06 eV) and +0.02 eV per-neighbor saddle term are declared model approximations; not environment-resolved NEB.",
        }
    for template in templates:
        _fields(template, TEMPLATE_FIELDS, "lattice channel template")
        length = _number(template["length_A"], "template.length_A", positive=True)
        direction_z = float(template["displacement_A"][2]) / length
        correction = 0.025 * (length - 2.5) + 0.020 * direction_z ** 2
        template["data"] = {}
        for q, baseline in ((0, 0.78), (2, 0.68)):
            key = str(q)
            sites = bundle["site_energies"][key]
            initial = sites[template["from_basis"]]["energy_eV"]
            final = sites[template["to_basis"]]["energy_eV"]
            saddle = max(initial, final) + baseline + correction
            template["data"][key] = {
                "barrier_forward_eV": saddle - initial,
                "barrier_reverse_eV": saddle - final,
                "attempt_forward_Hz": 1e12, "attempt_reverse_Hz": 1e12,
                "ts_fraction": 0.5, "status": "illustrative", "synthetic_data": True,
                "provenance": SYNTHETIC_PROVENANCE + f" For q={q}, S=max(eps_i,eps_j)+{baseline}+[0.025*(length_A-2.5)+0.020*(dz/length_A)^2] eV. Bf=S-eps_i; Br=S-eps_j. Shared 1e12 Hz prefactor and midpoint TS are synthetic assumptions.",
            }
    # Validate generated content through the same import gate as external fixtures.
    import_mock_bundle(bundle)
    return bundle


def import_mock_bundle(bundle: dict, *, expected_sha256: str | None = None) -> Catalog:
    """Convert consistent synthetic raw barriers into a validated Catalog.

    Inputs are never mutated.  ``expected_sha256`` optionally checks an externally
    retained content digest before conversion.  The result is always illustrative
    and preserves the input digest and source provenance.  This checksum does not
    authenticate who produced the data or establish scientific calibration.
    """
    _fields(bundle, {"schema_version", "metadata", "structure", "site_energies", "charge_models", "neb_paths"}, "bundle")
    if bundle["schema_version"] != SCHEMA_VERSION:
        raise ValueError(f"schema_version must be {SCHEMA_VERSION!r}")
    digest = source_digest(bundle)
    if expected_sha256 is not None:
        if not isinstance(expected_sha256, str) or expected_sha256 != digest:
            raise ValueError("source SHA256 mismatch; the bundle differs from the retained source digest")
    metadata = bundle["metadata"]
    _fields(metadata, META_REQUIRED, "metadata", META_OPTIONAL)
    _synthetic(metadata, "metadata", metadata=True)
    if metadata["material"] != "HfO2":
        raise ValueError("metadata.material must be 'HfO2' for the fixed conventional hafnia mock schema")
    _text(metadata["structure_id"], "metadata.structure_id")
    _text(metadata["energy_reference"], "metadata.energy_reference")
    for field in ("bundle_id", "generator", "charge_selection"):
        if field in metadata:
            _text(metadata[field], f"metadata.{field}")
    for field, unit in (("energy_unit", "eV"), ("length_unit", "angstrom"), ("frequency_unit", "Hz")):
        if metadata[field] != unit:
            raise ValueError(f"metadata.{field} must be {unit!r}; no implicit unit conversion is performed")
    _number(metadata["migration_cutoff_A"], "metadata.migration_cutoff_A", positive=True)
    if metadata.get("rate_convention", RATE_CONVENTION) != RATE_CONVENTION:
        raise ValueError("only shared-saddle classical energies with equal directional prefactors are supported")
    if metadata.get("vibrational_entropy_included", False) is not False:
        raise ValueError("vibrational entropy is unsupported by the mock classical-energy importer")
    charges = metadata["charges"]
    if (not isinstance(charges, list) or not charges
            or any(not isinstance(q, Integral) or isinstance(q, (bool, np.bool_)) for q in charges)
            or len(set(charges)) != len(charges)):
        raise ValueError("metadata.charges must be a nonempty list of distinct integral charges")
    keys = tuple(str(int(q)) for q in sorted(charges))
    reference = bundle["structure"]
    _fields(reference, {"cell_A", "hf_fractional", "oxygen_fractional", "oxygen_kinds", "provenance"}, "structure")
    _text(reference["provenance"], "structure.provenance")
    kinds = reference["oxygen_kinds"]
    if not isinstance(kinds, list) or len(kinds) != 8 or any(kind not in ("O3", "O4") for kind in kinds):
        raise ValueError("structure.oxygen_kinds must contain eight explicit O3/O4 labels")

    site_energies = {}
    _charge_mapping(bundle["site_energies"], keys, "site_energies")
    for key in keys:
        records = bundle["site_energies"][key]
        if not isinstance(records, list) or len(records) != 8:
            raise ValueError(f"site_energies[{key}] must contain exactly eight labelled records")
        values, seen = [None] * 8, set()
        for record in records:
            path = f"site_energies[{key}]"
            _fields(record, SITE_FIELDS, path)
            _synthetic(record, path)
            basis = _index(record["basis_id"], path + ".basis_id")
            if basis in seen:
                raise ValueError(f"duplicate site-energy basis ID {basis} for charge {key}")
            seen.add(basis)
            if record["kind"] != kinds[basis]:
                raise ValueError(f"site_energies[{key}][basis={basis}].kind disagrees with structure")
            values[basis] = _number(record["energy_eV"], f"site_energies[{key}][basis={basis}].energy_eV")
        site_energies[key] = values

    charge_models = {}
    _charge_mapping(bundle["charge_models"], keys, "charge_models")
    for key in keys:
        model = bundle["charge_models"][key]
        _fields(model, MODEL_FIELDS, f"charge_models[{key}]")
        _synthetic(model, f"charge_models[{key}]")
        charge_models[key] = {field: deepcopy(value) for field, value in model.items()
                              if field not in ("status", "synthetic_data")}

    channels, closure_errors = [], []
    paths = bundle["neb_paths"]
    if not isinstance(paths, list) or not paths:
        raise ValueError("neb_paths must be a nonempty complete path list")
    for index, path in enumerate(paths):
        location = f"neb_paths[{index}]"
        _fields(path, TEMPLATE_FIELDS | {"data"}, location)
        a = _index(path["from_basis"], location + ".from_basis")
        b = _index(path["to_basis"], location + ".to_basis")
        _charge_mapping(path["data"], keys, location + ".data")
        channel = {field: deepcopy(path[field]) for field in TEMPLATE_FIELDS}
        channel["data"] = {}
        for key in keys:
            data = path["data"][key]
            entry = location + f".data[{key}]"
            _fields(data, BARRIER_FIELDS, entry)
            _synthetic(data, entry)
            forward = _number(data["barrier_forward_eV"], entry + ".barrier_forward_eV", nonnegative=True)
            reverse = _number(data["barrier_reverse_eV"], entry + ".barrier_reverse_eV", nonnegative=True)
            prefactor_forward = _number(data["attempt_forward_Hz"], entry + ".attempt_forward_Hz", positive=True)
            prefactor_reverse = _number(data["attempt_reverse_Hz"], entry + ".attempt_reverse_Hz", positive=True)
            if prefactor_forward != prefactor_reverse:
                raise ValueError(f"{entry}: unequal directional prefactors are unsupported; one shared attempt_Hz is required")
            fraction = _number(data["ts_fraction"], entry + ".ts_fraction")
            if not 0 <= fraction <= 1:
                raise ValueError(f"{entry}.ts_fraction must lie in [0, 1]")
            initial, final = site_energies[key][a], site_energies[key][b]
            saddle_forward, saddle_reverse = initial + forward, final + reverse
            if not math.isfinite(saddle_forward) or not math.isfinite(saddle_reverse):
                raise ValueError(f"{entry}: derived saddle energy is non-finite")
            closure_error = abs(saddle_forward - saddle_reverse)
            closure_errors.append(closure_error)
            if closure_error > ENERGY_CLOSURE_TOL_EV:
                raise ValueError(
                    f"{entry}: inconsistent forward/reverse barriers: eps_i+Bf={saddle_forward:.12g} eV, "
                    f"eps_j+Br={saddle_reverse:.12g} eV (tolerance {ENERGY_CLOSURE_TOL_EV:g} eV); "
                    "a common saddle and Bf-Br=eps_j-eps_i are required"
                )
            channel["data"][key] = {
                "saddle_eV": saddle_forward, "attempt_Hz": prefactor_forward,
                "ts_fraction": fraction, "status": "illustrative", "provenance": data["provenance"],
            }
        channels.append(channel)

    catalog_metadata = {
        "material": metadata["material"], "structure_id": metadata["structure_id"],
        "calibration_status": "illustrative", "synthetic_data": True,
        "energy_unit": "eV", "length_unit": "angstrom",
        "migration_cutoff_A": float(metadata["migration_cutoff_A"]), "channel_exclusions": {},
        "provenance": metadata["provenance"], "energy_reference": metadata["energy_reference"],
        "rate_convention": RATE_CONVENTION, "vibrational_entropy_included": False,
        "prefactor_convention": "One equal forward/reverse attempt_Hz per charge/channel; no vibrational entropy.",
        "mock_source": {
            "schema_version": SCHEMA_VERSION, "sha256": digest, "synthetic_data": True,
            "provenance": metadata["provenance"], "source_metadata": deepcopy(metadata),
            "digest_scope": "Complete canonical UTF-8 JSON bundle, sorted object keys, compact separators, no added newline.",
            "interpretation": "Content checksum only; not source authentication or evidence of completed DFT/NEB.",
        },
        "site_energy_provenance": deepcopy(bundle["site_energies"]),
        "import_consistency": {
            "conversion": "saddle_eV=epsilon_initial+barrier_forward_eV; independent reverse closure checked",
            "energy_closure_tolerance_eV": ENERGY_CLOSURE_TOL_EV,
            "max_energy_closure_error_eV": max(closure_errors, default=0.0),
            "checked_channel_charge_records": len(closure_errors),
            "equal_directional_prefactors": True,
        },
    }
    if "charge_selection" in metadata:
        catalog_metadata["charge_selection"] = metadata["charge_selection"]
    catalog = Catalog({
        "schema_version": CATALOG_SCHEMA, "metadata": catalog_metadata,
        "structure": deepcopy(reference), "site_energies_eV": site_energies,
        "charge_models": charge_models, "channels": channels,
    })
    # Includes structure hash, full geometric coverage, canonical IDs/displacements,
    # pair model, provenance and saddle/endpoint bounds.  No inferred exclusions.
    catalog.validate(allow_illustrative=True)
    return catalog


def _unique_json_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _invalid_constant(value):
    raise ValueError(f"non-finite JSON constant {value!r} is forbidden")


def load_mock_bundle(path, *, expected_sha256=None) -> dict:
    """Read strict JSON and validate the complete synthetic import contract."""
    with Path(path).open(encoding="utf-8") as handle:
        bundle = json.load(handle, object_pairs_hook=_unique_json_pairs,
                           parse_constant=_invalid_constant)
    import_mock_bundle(bundle, expected_sha256=expected_sha256)
    return bundle


def write_mock_bundle(bundle: dict, path) -> None:
    """Validate and atomically write a synthetic source bundle as strict JSON."""
    import_mock_bundle(bundle)
    payload = json.dumps(bundle, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
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


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Generate/import explicitly synthetic raw data; no DFT calculation or simulation is run.")
    commands = parser.add_subparsers(dest="command", required=True)
    generate = commands.add_parser("generate", help="write the deterministic synthetic raw bundle")
    generate.add_argument("--output", required=True)
    convert = commands.add_parser("import", help="validate raw barrier pairs and write an illustrative catalogue")
    convert.add_argument("--input", required=True)
    convert.add_argument("--output", required=True)
    convert.add_argument("--expected-sha256")
    args = parser.parse_args(argv)
    try:
        if args.command == "generate":
            bundle = make_mock_bundle()
            write_mock_bundle(bundle, args.output)
        else:
            bundle = load_mock_bundle(args.input, expected_sha256=args.expected_sha256)
            import_mock_bundle(bundle, expected_sha256=args.expected_sha256).write(args.output)
        print(json.dumps({"status": "illustrative", "synthetic_data": True,
                          "output": str(args.output), "source_sha256": source_digest(bundle),
                          "paths": len(bundle["neb_paths"]), "charges": bundle["metadata"]["charges"]}, indent=2))
    except (ValueError, OSError, TypeError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
