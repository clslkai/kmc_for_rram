"""Command-line entry points for a traceable DFT-parameterization workflow."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Sequence

from .catalog import Catalog, make_catalog
from .device import DeviceConfig, ElectrodeConfig
from .coupled import CoupledConfig
from .coupled_workflow import run_coupled
from .geometry import Lattice
from .model import SimulationConfig
from .workflow import ABLATIONS, export_neb, run_simulation, run_study


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be a nonnegative integer")
    return parsed


def _reference(path: str | None) -> dict | None:
    if path is None:
        return None
    with Path(path).open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("A structure reference must be a JSON object.")
    return value


def _geometry_options(parser: argparse.ArgumentParser, *, optional_cutoffs: bool = False) -> None:
    parser.add_argument("--repeats", nargs=3, type=_positive_int, default=(3, 3, 3), metavar=("NX", "NY", "NZ"),
                        help="Conventional-cell repeats; must be large enough for unique periodic hop images.")
    parser.add_argument("--migration-cutoff", type=float, default=None if optional_cutoffs else 3.2,
                        help="Candidate ionic-hop cutoff in angstrom (default 3.2, or catalogue value).")
    parser.add_argument("--interaction-cutoff", type=float, default=None if optional_cutoffs else 3.6,
                        help="Independent vacancy interaction cutoff in angstrom (default 3.6, or catalogue value).")


def _run_options(parser: argparse.ArgumentParser, *, study: bool = False) -> None:
    parser.add_argument("--catalog", required=True, help="Validated JSON catalogue; pending parameters cannot run.")
    parser.add_argument("--allow-illustrative", action="store_true", help="Explicitly permit synthetic illustrative parameters.")
    parser.add_argument("--repeats", nargs=3, type=_positive_int, default=(3, 3, 3), metavar=("NX", "NY", "NZ"))
    parser.add_argument("--boundary", choices=("bulk", "slab"), default=None,
                        help="Default bulk; --device selects slab and rejects an explicit bulk boundary.")
    parser.add_argument("--device", help="Electrode boundary JSON; prescribes voltage on ideal bottom/top planes.")
    parser.add_argument("--temperature", "--T", "-T", dest="temperature_K", type=float, default=600.0, metavar="K")
    parser.add_argument("--seed", type=_nonnegative_int, default=17, help="Run seed, or base seed for independent study streams.")
    parser.add_argument("--nvac", "--n-vacancies", dest="n_vacancies", type=_positive_int, default=1)
    parser.add_argument("--field", nargs=3, type=float, default=None, metavar=("FX", "FY", "FZ"),
                        help="Prescribed uniform field in V/nm (default zero); mutually exclusive with --device.")
    parser.add_argument("--duration", "--duration-s", dest="duration_s", type=float, required=True,
                        help="Requested physical duration in seconds; event caps are reported as censoring.")
    parser.add_argument("--max-events", type=_positive_int, default=100000)
    parser.add_argument("--sample-every", type=_positive_int, default=100)
    parser.add_argument("--pair-scale", type=float, default=1.0)
    parser.add_argument("--environment-scale", type=float, default=1.0)
    parser.add_argument("--output", required=True, help="New or empty output directory.")
    parser.add_argument("--plot", action="store_true", help="Write an offline overview.png for each trajectory.")
    if study:
        parser.add_argument("--charges", nargs="+", type=int, help="Fixed charge branches to compare (default: every catalogue branch).")
        parser.add_argument("--ablations", nargs="+", choices=tuple(ABLATIONS), default=list(ABLATIONS),
                            help="Each group receives independent trajectory seeds. Closure is never silently changed.")
        parser.add_argument("--replicas", type=_positive_int, default=8)
        parser.add_argument("--confidence", type=float, default=0.95)
        parser.add_argument("--bootstrap-samples", type=_positive_int, default=2000)
    else:
        parser.add_argument("--charge", type=int, default=2)
        parser.add_argument("--vacancy-indices", nargs="+", type=_nonnegative_int,
                            help="Optional explicit initial vacancy sites; overrides random placement and --nvac.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="dftkmc", description="Fixed-charge HfO2 lattice KMC and traceable DFT/NEB preparation.")
    subcommands = parser.add_subparsers(dest="command", required=True)
    coupled_template = subcommands.add_parser("coupled-template", help="Write a synthetic finite-reservoir / electronic-network device configuration.")
    coupled_template.add_argument("--output", required=True)
    coupled_template.add_argument("--top-voltage", type=float, default=.5)
    coupled_template.add_argument("--bottom-voltage", type=float, default=0.)
    coupled_template.add_argument("--repeats", nargs=3, type=_positive_int, default=(2, 2, 4))
    device_run = subcommands.add_parser("device-run", help="Run finite oxygen exchange and configuration-dependent effective electronic transport.")
    device_run.add_argument("--catalog", required=True)
    device_run.add_argument("--config", required=True, help="Complete coupled device JSON, generated with coupled-template.")
    device_run.add_argument("--allow-illustrative", action="store_true")
    device_run.add_argument("--nvac", type=_nonnegative_int, default=0)
    device_run.add_argument("--vacancy-indices", nargs="*", type=_nonnegative_int)
    device_run.add_argument("--seed", type=_nonnegative_int, help="Optional explicit override of the seed in the configuration.")
    timing = device_run.add_mutually_exclusive_group(required=True)
    timing.add_argument("--duration", type=float, help="Constant-bias physical duration in seconds.")
    timing.add_argument("--pulse", nargs=2, type=float, action="append", metavar=("BIAS_V", "DURATION_S"),
                        help="Repeat for piecewise-constant requested top-minus-bottom voltage segments.")
    device_run.add_argument("--max-events", type=_positive_int, default=20000)
    device_run.add_argument("--sample-every", type=_positive_int, default=10)
    device_run.add_argument("--plot", action="store_true")
    device_run.add_argument("--output", required=True)
    device_template = subcommands.add_parser("device-template", help="Write an ideal blocking electrode/oxide/electrode boundary JSON.")
    device_template.add_argument("--output", required=True, help="New device JSON path.")
    device_template.add_argument("--bottom-voltage", type=float, default=0.0, metavar="V")
    device_template.add_argument("--top-voltage", type=float, default=1.0, metavar="V")
    device_template.add_argument("--bottom-contact-depth", type=float, default=2.0, metavar="ANGSTROM")
    device_template.add_argument("--top-contact-depth", type=float, default=2.0, metavar="ANGSTROM")
    template = subcommands.add_parser("template", help="Write an illustrative catalogue or an unresolved DFT worklist.")
    template.add_argument("--status", choices=("pending", "illustrative"), default="pending")
    template.add_argument("--charges", nargs="+", type=int, default=(0, 2))
    template.add_argument("--reference", help="JSON conventional-cell reference; defaults to the inherited unrelaxed monoclinic reference.")
    template.add_argument("--output", required=True, help="New catalogue JSON path.")
    _geometry_options(template)
    audit = subcommands.add_parser("audit", help="Inspect coverage, provenance, values and scientific calibration without running KMC.")
    audit.add_argument("catalog_path", nargs="?", help="Catalogue JSON path (or use --catalog).")
    audit.add_argument("--catalog", dest="catalog_option")
    audit.add_argument("--output", help="Optionally save the complete audit report as JSON.")
    _run_options(subcommands.add_parser("run", help="Run one trajectory and export states, event-level parameters and provenance."))
    _run_options(subcommands.add_parser("study", help="Compare independent fixed-charge/environment replicas at a common physical time."), study=True)
    neb = subcommands.add_parser("export-neb", help="Export unrelaxed atom-matched POSCAR endpoints; submits no jobs.")
    source = neb.add_mutually_exclusive_group()
    source.add_argument("--catalog", help="Use this catalogue's geometry and branch labels; unresolved energies are allowed for preparation.")
    source.add_argument("--reference", help="Use a supplied conventional-cell reference JSON.")
    neb.add_argument("--charges", nargs="+", type=int, help="Fixed charge branches (default: catalogue branches, otherwise 0 and 2).")
    neb.add_argument("--channels", nargs="+", help="Optional known canonical channel IDs; otherwise export every candidate channel.")
    neb.add_argument("--output", required=True, help="New or empty endpoint/worklist directory.")
    _geometry_options(neb, optional_cutoffs=True)
    return parser


def _print(value: object) -> None:
    print(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False))


def _config(args: argparse.Namespace) -> SimulationConfig:
    device = DeviceConfig.load(args.device) if args.device else None
    if device is not None and args.field is not None:
        raise ValueError("--device and --field are mutually exclusive, including an explicit zero field")
    if device is not None and args.boundary == "bulk":
        raise ValueError("--device requires a slab; --boundary bulk would connect the two electrodes periodically")
    boundary = args.boundary or ("slab" if device is not None else "bulk")
    return SimulationConfig(repeats=tuple(args.repeats), charge=getattr(args, "charge", 2),
                            temperature_K=args.temperature_K, seed=args.seed, boundary=boundary,
                            field_V_nm=tuple(args.field) if args.field is not None else (0., 0., 0.),
                            pair_scale=args.pair_scale, environment_scale=args.environment_scale,
                            device=device)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "coupled-template":
            config = CoupledConfig(repeats=tuple(args.repeats), top_potential_V=args.top_voltage,
                                   bottom_potential_V=args.bottom_voltage)
            config.write(args.output)
            _print({"output": args.output, "calibration_status": "illustrative",
                    "notice": "Finite oxygen inventory and effective resistor transport; frozen-field ionic kinetics, no Poisson or heat solver. Use --allow-illustrative to run."})
            return 0
        if args.command == "device-run":
            from dataclasses import replace
            catalog = Catalog.load(args.catalog, allow_illustrative=args.allow_illustrative)
            config = CoupledConfig.load(args.config)
            if args.seed is not None:
                config = replace(config, seed=args.seed)
            result = run_coupled(catalog, config, duration_s=args.duration, pulses=args.pulse,
                                 n_vacancies=args.nvac, vacancy_indices=args.vacancy_indices,
                                 max_events=args.max_events, sample_every=args.sample_every,
                                 output=args.output, plot=args.plot, allow_illustrative=args.allow_illustrative)
            _print({"output": args.output, "termination": result["termination"],
                    "calibration_status": result["calibration_status"], "conservation": result["conservation"],
                    "first_model_formation_time_s": result["formation"]["first_criterion_time_s"],
                    "final_conductance_S": result["final"]["equivalent_conductance_S"],
                    "final_vacancies": result["final"]["vacancies"]})
            return 0 if result["termination"]["target_reached"] else 3
        if args.command == "device-template":
            device = DeviceConfig(
                bottom=ElectrodeConfig(potential_V=args.bottom_voltage, contact_depth_A=args.bottom_contact_depth),
                top=ElectrodeConfig(potential_V=args.top_voltage, contact_depth_A=args.top_contact_depth),
            )
            destination = Path(args.output)
            if destination.exists():
                raise FileExistsError(f"Device path already exists: {destination}")
            device.write(destination)
            _print({"output": str(destination), "device": device.to_dict(),
                    "notice": "Ideal planar blocking electrodes; linear dielectric potential; no electronic transport, space charge or interface oxygen exchange."})
            return 0
        if args.command == "template":
            destination = Path(args.output)
            if destination.exists():
                raise FileExistsError(f"Catalogue path already exists: {destination}")
            lattice = Lattice(repeats=tuple(args.repeats), boundary="bulk", migration_cutoff_A=args.migration_cutoff,
                              interaction_cutoff_A=args.interaction_cutoff, reference=_reference(args.reference))
            catalogue = make_catalog(lattice, status=args.status, charges=tuple(args.charges))
            catalogue.write(destination)
            _print({"output": str(destination), "calibration_status": args.status, "structure_id": lattice.structure_id,
                    "n_channels": len(lattice.channel_templates()), "charges": list(catalogue.charges),
                    "notice": "Synthetic non-literature parameters; use --allow-illustrative to run." if args.status == "illustrative"
                              else "Unresolved DFT worklist. Missing data deliberately block simulation."})
            return 0
        if args.command == "audit":
            if bool(args.catalog_path) == bool(args.catalog_option):
                raise ValueError("Supply exactly one catalogue path, positional or with --catalog.")
            catalogue = Catalog.load(args.catalog_path or args.catalog_option, validate=False)
            report = catalogue.audit()
            if args.output:
                destination = Path(args.output)
                if destination.exists():
                    raise FileExistsError(f"Audit output already exists: {destination}")
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
            _print(report)
            return 0 if report["ready"] else 1
        if args.command in ("run", "study"):
            catalogue = Catalog.load(args.catalog, allow_illustrative=args.allow_illustrative)
            config = _config(args)
            common = dict(n_vacancies=args.n_vacancies, duration_s=args.duration_s, max_events=args.max_events,
                          sample_every=args.sample_every, output=args.output, plot=args.plot,
                          allow_illustrative=args.allow_illustrative)
            if args.command == "run":
                result = run_simulation(catalogue, config, vacancy_indices=args.vacancy_indices, **common)
                _print({"output": args.output, "status": result["status"], "target_reached": result["target_reached"],
                        "time_s": result["time_s"], "requested_time_s": args.duration_s, "events": result["events"],
                        "boundary": config.boundary, "effective_field_V_nm": result["effective_field_V_nm"],
                        "device_boundary": result.get("device_boundary"),
                        "calibration_status": catalogue.data["metadata"]["calibration_status"]})
                return 0 if result["target_reached"] else 3
            result = run_study(catalogue, config, replicas=args.replicas, charges=args.charges, ablations=args.ablations,
                               confidence=args.confidence, bootstrap_samples=args.bootstrap_samples, **common)
            complete = all(group["statistics"]["all_replicas_complete"] for group in result["groups"])
            _print({"output": args.output, "groups": len(result["groups"]), "all_replicas_complete": complete,
                    "requested_time_s": args.duration_s, "calibration_status": result["calibration_status"],
                    "warnings": result["warnings"]})
            return 0 if complete else 3
        if args.command == "export-neb":
            if args.catalog:
                catalogue = Catalog.load(args.catalog, validate=False)
                reference = catalogue.structure
                charges = args.charges if args.charges is not None else catalogue.charges
                migration_cutoff = args.migration_cutoff
                if migration_cutoff is None:
                    migration_cutoff = catalogue.data.get("metadata", {}).get("migration_cutoff_A")
                if migration_cutoff is None:
                    migration_cutoff = max(float(channel["length_A"]) for channel in catalogue.data["channels"]) + 1e-8
                pair_cutoff = args.interaction_cutoff
                if pair_cutoff is None:
                    pair_cutoff = max(float(catalogue.charge_model(q)["pair_cutoff_A"]) for q in charges)
            else:
                reference = _reference(args.reference)
                charges = args.charges if args.charges is not None else (0, 2)
                migration_cutoff = args.migration_cutoff if args.migration_cutoff is not None else 3.2
                pair_cutoff = args.interaction_cutoff if args.interaction_cutoff is not None else 3.6
            lattice = Lattice(repeats=tuple(args.repeats), boundary="bulk", migration_cutoff_A=migration_cutoff,
                              interaction_cutoff_A=pair_cutoff, reference=reference)
            if args.catalog and lattice.structure_id != catalogue.structure_id:
                raise ValueError("Catalogue structure identity does not match its reference geometry.")
            result = export_neb(lattice, output=args.output, charges=charges, channel_ids=args.channels,
                                source_catalog=catalogue if args.catalog else None)
            _print({"output": args.output, **result})
            return 0
    except (ValueError, OSError, KeyError, TypeError, OverflowError, RuntimeError) as exc:
        print(f"dftkmc: {exc}", file=sys.stderr)
        return 2
    raise RuntimeError("Unreachable command")
