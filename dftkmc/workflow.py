"""Reproducible run, controlled-study, and unrelaxed NEB worklist workflows."""
from __future__ import annotations

import csv
from dataclasses import asdict, is_dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .analysis import compare_groups, reached_target, summarize_replicas
from .catalog import Catalog, make_catalog
from .device import DeviceBoundary
from .geometry import Lattice
from .model import SimulationConfig, Simulator


ABLATIONS = {
    "noninteracting": (0.0, 0.0),
    "pair-only": (1.0, 0.0),
    "pair-environment": (1.0, 1.0),
}


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def _write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True, default=_json_default, allow_nan=False) + "\n", encoding="utf-8")


def _config_dict(config: SimulationConfig) -> dict[str, Any]:
    return asdict(config) if is_dataclass(config) else dict(vars(config))


def catalog_digest(catalog: Catalog) -> str:
    canonical = json.dumps(catalog.data, sort_keys=True, separators=(",", ":"), default=_json_default, allow_nan=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _graph_coverage(catalog: Catalog) -> dict[str, Any]:
    """Carry the declared migration domain into every standalone result."""
    metadata = catalog.data.get("metadata", {})
    exclusions = metadata.get("channel_exclusions", {})
    if not isinstance(exclusions, dict) or any(
        not isinstance(key, str) or not key.strip() or not isinstance(reason, str) or not reason.strip()
        for key, reason in exclusions.items()
    ):
        raise ValueError("channel_exclusions must map channel IDs to nonempty truncation reasons.")
    cutoff = metadata.get("migration_cutoff_A")
    cutoff_source = "declared in catalogue metadata"
    if cutoff is None:
        cutoff = max(float(channel["length_A"]) for channel in catalog.data["channels"]) + 1e-8
        cutoff_source = "inferred through longest supplied channel plus numerical tolerance"
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("Migration coverage requires a positive finite cutoff.")
    warning = (
        "Dynamics and tracer estimates are conditional on the declared truncated migration graph. "
        "Excluded channels are not simulated; these results cannot establish full diffusion in the material."
        if exclusions else ""
    )
    return {"migration_cutoff_A": float(cutoff), "migration_cutoff_source": cutoff_source,
            "channel_exclusions": dict(sorted(exclusions.items())),
            "n_excluded_channels": len(exclusions), "conditional_on_truncated_graph": bool(exclusions),
            "coverage_warning": warning}


def _new_output(path: str | Path) -> Path:
    output = Path(path)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError(f"Output must be a new or empty directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    return output


def _csv_value(value: Any) -> Any:
    if isinstance(value, (list, tuple, dict, np.ndarray)):
        return json.dumps(value, default=_json_default, sort_keys=True, allow_nan=False, separators=(",", ":"))
    if isinstance(value, np.generic):
        return value.item()
    return value


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], preferred: Sequence[str] = ()) -> None:
    fields = list(dict.fromkeys(preferred))
    for row in rows:
        fields.extend(key for key in row if key not in fields)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def _event_rows(events: Sequence[Mapping[str, Any]], catalog: Catalog, config: SimulationConfig,
                initial_labels: np.ndarray, device: DeviceBoundary | None = None) -> list[dict[str, Any]]:
    result = []
    labels = np.asarray(initial_labels, dtype=int).copy()
    previous_time = 0.0
    digest = catalog_digest(catalog)
    for index, raw in enumerate(events, start=1):
        row = dict(raw)
        row.setdefault("event_index", index)
        source, target = int(row["source"]), int(row["target"])
        if labels[source] < 0 or labels[target] >= 0:
            raise ValueError("Engine event trace violates vacancy identity/conservation.")
        row["vacancy_id"] = int(labels[source])
        labels[target], labels[source] = labels[source], -1
        if "time_s" not in row:
            raise ValueError("Engine event output is missing its physical event time.")
        row["time_before_s"] = previous_time
        row["dt_s"] = float(row["time_s"]) - previous_time
        if row["dt_s"] <= 0:
            raise ValueError("Event times must be strictly increasing.")
        previous_time = float(row["time_s"])
        row["catalog_sha256"] = digest
        row["structure_id"] = catalog.structure_id
        row["seed"] = int(config.seed)
        row["temperature_K"] = float(config.temperature_K)
        row["field_V_nm"] = device.field_V_nm.tolist() if device is not None else list(config.field_V_nm)
        if device is not None:
            row.update(
                source_potential_V=float(device.site_potential_V[source]),
                target_potential_V=float(device.site_potential_V[target]),
                source_region=str(device.region[source]), target_region=str(device.region[target]),
                applied_voltage_V=device.config.top.potential_V-device.config.bottom.potential_V,
                potential_model=device.config.potential_model,
            )
        row["pair_scale"] = float(config.pair_scale)
        row["environment_scale"] = float(config.environment_scale)
        channel_id = row.get("channel_id")
        if channel_id is None:
            raise ValueError("Engine event output is missing channel_id; cannot produce a traceable event record.")
        charge = int(row.get("q", config.charge))
        datum = catalog.channel(str(channel_id), charge)
        row["channel_status"] = datum["status"]
        row["channel_provenance"] = datum["provenance"]
        row["catalog_saddle_eV"] = datum["saddle_eV"]
        row["catalog_attempt_Hz"] = datum["attempt_Hz"]
        row["ts_fraction"] = datum["ts_fraction"]
        if "time_before_s" not in row and "time_s" in row and "dt_s" in row:
            row["time_before_s"] = float(row["time_s"]) - float(row["dt_s"])
        result.append(row)
    return result


def _plot_overview(records: Sequence[Mapping[str, Any]], output: Path, calibration_status: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    times = [r["time_s"] for r in records]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for axis, field, label in zip(axes.flat,
                                ("energy_eV", "msd_A2", "pair_count", "largest_cluster"),
                                ("Configurational energy (eV)", "Raw tracer MSD (angstrom squared)",
                                 "Vacancy pair count", "Largest interaction-graph cluster")):
        axis.plot(times, [r[field] for r in records], lw=1.5)
        axis.set_xlabel("Physical time (s)")
        axis.set_ylabel(label)
        axis.grid(alpha=0.2)
    calibration_label = ("ILLUSTRATIVE PARAMETERS / NOT DFT CALIBRATED"
                         if calibration_status == "illustrative"
                         else f"Catalogue calibration_status: {calibration_status}")
    fig.suptitle(calibration_label + "\nFixed-charge KMC trajectory — sampled states, not independent replicas")
    fig.savefig(output / "overview.png", dpi=160)
    plt.close(fig)


def _plot_device_boundary(simulator: Simulator, trajectory: Mapping[str, Any], output: Path,
                          calibration_status: str) -> None:
    """Plot actual geometry, prescribed voltage and sampled vacancy statistics."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    device = simulator.device
    if device is None:
        raise ValueError("A device boundary is required for a device plot")
    records = trajectory["records"]
    initial, final = trajectory["states"][0], trajectory["states"][-1]
    depth_nm, thickness_nm = device.depth_A/10., device.thickness_A/10.
    in_plane = simulator.lattice.cell_A[0]/np.linalg.norm(simulator.lattice.cell_A[0])
    lateral_nm = simulator.lattice.oxygen_A@in_plane/10.
    times_us = np.asarray([row["time_s"] for row in records])*1e6
    bottom_color, top_color, interior_color = "#16877d", "#8061a8", "#8b969f"
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.3))
    for axis in axes.flat:
        axis.grid(alpha=0.18)
        axis.spines[["top", "right"]].set_visible(False)

    axis = axes[0, 0]
    axis.axhspan(0., device.config.bottom.contact_depth_A/10., color=bottom_color, alpha=.14)
    axis.axhspan(thickness_nm-device.config.top.contact_depth_A/10., thickness_nm, color=top_color, alpha=.14)
    axis.scatter(lateral_nm, depth_nm, s=5, color="#bdc7cc", alpha=.55, label="O lattice sites")
    axis.scatter(lateral_nm[initial], depth_nm[initial], s=45, color="#236f9c", marker="x", label="Initial vacancies", zorder=4)
    axis.scatter(lateral_nm[final], depth_nm[final], s=50, facecolors="none", edgecolors="#cc652f",
                 linewidths=1.6, label="Final vacancies", zorder=5)
    axis.axhline(0., color=bottom_color, lw=2.4)
    axis.axhline(thickness_nm, color=top_color, lw=2.4)
    centre = .5*(lateral_nm.min()+lateral_nm.max())
    axis.text(centre, -.035*thickness_nm, f"Bottom electrode: {device.config.bottom.potential_V:g} V",
              ha="center", va="top", color=bottom_color, fontsize=9)
    axis.text(centre, 1.035*thickness_nm, f"Top electrode: {device.config.top.potential_V:g} V",
              ha="center", va="bottom", color=top_color, fontsize=9)
    axis.set_ylim(-.12*thickness_nm, 1.12*thickness_nm)
    axis.set_xlabel("Projection along the first in-plane cell vector (nm)")
    axis.set_ylabel("Depth along the plane normal (nm)")
    axis.set_title("A  |  Finite oxide slab and contact regions", loc="left", fontsize=11)
    axis.legend(loc="center right", fontsize=8, framealpha=.9)

    axis = axes[0, 1]
    depths = np.linspace(0., device.thickness_A, 101)
    potentials = device.potential_V(depths[:, None]*device.normal)
    axis.plot(potentials, depths/10., color="#273f54", lw=2.2)
    axis.scatter([device.config.bottom.potential_V, device.config.top.potential_V], [0., thickness_nm],
                 color=[bottom_color, top_color], s=45, zorder=4)
    axis.axhspan(0., device.config.bottom.contact_depth_A/10., color=bottom_color, alpha=.1)
    axis.axhspan(thickness_nm-device.config.top.contact_depth_A/10., thickness_nm, color=top_color, alpha=.1)
    axis.set_xlabel("Prescribed electrostatic potential (V)")
    axis.set_ylabel("Depth along the plane normal (nm)")
    axis.set_ylim(-.03*thickness_nm, 1.03*thickness_nm)
    axis.set_title("B  |  Linear Dirichlet potential", loc="left", fontsize=11)
    axis.text(.04, .94, f"Normal field: {float(device.field_V_nm@device.normal):+.4f} V/nm\n"
              "No space-charge or current feedback", transform=axis.transAxes,
              fontsize=9, va="top", bbox=dict(facecolor="white", edgecolor="none", alpha=.9))

    axis = axes[1, 0]
    minimum = np.asarray([row["min_vacancy_depth_A"] for row in records])/10.
    maximum = np.asarray([row["max_vacancy_depth_A"] for row in records])/10.
    mean = np.asarray([row["mean_vacancy_depth_A"] for row in records])/10.
    axis.fill_between(times_us, minimum, maximum, color="#236f9c", alpha=.12, label="Min-max depth")
    axis.plot(times_us, mean, color="#236f9c", lw=1.7, label="Mean vacancy depth")
    axis.axhline(0., color=bottom_color, lw=1.1)
    axis.axhline(thickness_nm, color=top_color, lw=1.1)
    axis.set_ylim(-.03*thickness_nm, 1.03*thickness_nm)
    axis.set_xlabel("Physical time (microseconds)")
    axis.set_ylabel("Vacancy depth (nm)")
    axis.set_title("C  |  Vacancy motion inside blocking boundaries", loc="left", fontsize=11)
    axis.legend(loc="best", fontsize=8)

    axis = axes[1, 1]
    for key, label, color in (("bottom_contact_vacancies", "Bottom contact region", bottom_color),
                              ("interior_vacancies", "Oxide interior", interior_color),
                              ("top_contact_vacancies", "Top contact region", top_color)):
        axis.step(times_us, [row[key] for row in records], where="post", label=label, color=color, lw=1.5)
    axis.set_xlabel("Physical time (microseconds)")
    axis.set_ylabel("Number of vacancies")
    axis.set_ylim(-.25, int(initial.sum())+.5)
    axis.set_title("D  |  Contact occupancy; total vacancies conserved", loc="left", fontsize=11)
    axis.legend(loc="best", fontsize=8)

    title = "Electrode / HfO2 / electrode boundary"
    fig.suptitle(title, fontsize=17, weight="bold", y=.985)
    fig.text(.5, .945,
             f"Oxide thickness {thickness_nm:.4f} nm  |  fixed q={simulator.cfg.charge}  |  "
             f"{simulator.cfg.temperature_K:g} K  |  {int(initial.sum())} vacancies  |  "
             f"{simulator.event_count} events  |  {trajectory['status']}",
             ha="center", fontsize=10, color="#44535f")
    label = "ILLUSTRATIVE PARAMETERS / NOT DFT CALIBRATED" if calibration_status == "illustrative" else f"Catalogue status: {calibration_status}"
    fig.text(.5, .015, label + "  |  Sampled ionic states; no electronic conduction or filament claim.",
             ha="center", fontsize=9, color="#65513c")
    fig.tight_layout(rect=(.01, .04, .99, .925), h_pad=2.1, w_pad=2.)
    fig.savefig(output / "device_boundary.png", dpi=180)
    plt.close(fig)


def run_simulation(
    catalog: Catalog, config: SimulationConfig, *, n_vacancies: int = 1,
    duration_s: float, max_events: int = 100000, sample_every: int = 100,
    output: str | Path | None = None, plot: bool = False,
    vacancy_indices: Sequence[int] | None = None, allow_illustrative: bool = False,
) -> dict[str, Any]:
    """Execute one run and optionally save the complete numerical provenance.

    Return a lightweight replica record with final displacements for subsequent
    independent-replica analysis.  The detailed trajectories remain in files.
    """
    if not np.isfinite(duration_s) or duration_s <= 0:
        raise ValueError("duration_s must be finite and positive.")
    if max_events < 1 or sample_every < 1:
        raise ValueError("max_events and sample_every must be positive integers.")
    catalog.validate(allow_illustrative=allow_illustrative)
    coverage = _graph_coverage(catalog)
    simulator = Simulator(catalog, config)
    simulator.initialize(n_vacancies=n_vacancies, vacancy_indices=vacancy_indices)
    directory = _new_output(output) if output is not None else None
    initial_indices = np.flatnonzero(simulator.vacancy).tolist()
    initial_labels = simulator.labels.copy()
    initial_sites_by_id = sorted(initial_indices, key=lambda site: int(initial_labels[site]))
    trajectory = simulator.run(duration_s=duration_s, max_events=max_events, sample_every=sample_every)
    records = trajectory["records"]
    if not records:
        raise ValueError("Engine returned no initial/final snapshot.")
    if len(records) != len(trajectory["states"]) or len(records) != len(trajectory["displacements"]):
        raise ValueError("Engine state, displacement, and snapshot samples must have identical lengths.")
    final = dict(records[-1])
    final_time = float(final["time_s"])
    complete = reached_target(final_time, duration_s)
    status = str(trajectory["status"])
    n_actual = int(np.count_nonzero(simulator.vacancy))
    replica: dict[str, Any] = {
        "seed": int(config.seed), "status": status, "time_s": final_time,
        "target_reached": complete, "censored": not complete,
        "requested_time_s": float(duration_s), "events": int(simulator.event_count),
        "n_vacancies": n_actual, "initial_vacancy_indices": initial_indices,
        "final": final, "displacements_A": np.asarray(simulator.displacements_A).copy(),
        "output": str(directory) if directory is not None else None,
        "effective_field_V_nm": simulator.field_V_nm.tolist(),
        **coverage,
    }
    if simulator.device is not None:
        replica["device_boundary"] = simulator.device.summary()
    if directory is not None:
        events = _event_rows(trajectory["events"], catalog, config, initial_labels, simulator.device)
        _write_csv(directory / "trace.csv", records,
                   ("time_s", "events", "vacancies", "energy_eV", "mean_displacement_A", "msd_A2", "msd_tensor_A2", "pair_count", "largest_cluster"))
        _write_csv(directory / "events.csv", events,
                   ("event_index", "time_s", "time_before_s", "dt_s", "vacancy_id", "source", "target", "channel_id", "direction", "q", "rate_Hz", "barrier_eV", "delta_eV", "env_n3", "env_n4", "displacement_A"))
        device_arrays = {}
        if simulator.device is not None:
            device = simulator.device
            device_arrays = dict(
                cell_A=simulator.lattice.cell_A, oxygen_A=simulator.lattice.oxygen_A,
                hf_A=simulator.lattice.hf_A, pbc=np.asarray(simulator.lattice.pbc),
                potential_V=device.site_potential_V,
                potential_relative_bottom_V=device.relative_potential_V,
                device_depth_A=device.depth_A, bottom_contact=device.bottom_mask, top_contact=device.top_mask,
            )
            _write_json(directory / "device_boundary.json", device.summary())
            _write_csv(directory / "device_sites.csv", [
                dict(site_id=i, kind=str(simulator.lattice.kind[i]),
                     x_A=float(position[0]), y_A=float(position[1]), z_A=float(position[2]),
                     depth_A=float(device.depth_A[i]), potential_V=float(device.site_potential_V[i]),
                     potential_relative_bottom_V=float(device.relative_potential_V[i]),
                     region=str(device.region[i]), bottom_contact=bool(device.bottom_mask[i]),
                     top_contact=bool(device.top_mask[i]))
                for i, position in enumerate(simulator.lattice.oxygen_A)
            ])
        np.savez_compressed(
            directory / "states.npz",
            time_s=np.asarray([r["time_s"] for r in records], dtype=float),
            events=np.asarray([r["events"] for r in records], dtype=np.int64),
            vacancy=np.asarray(trajectory["states"], dtype=bool),
            displacements_A=np.asarray(trajectory["displacements"], dtype=float),
            initial_vacancy_indices=np.asarray(initial_indices, dtype=np.int64),
            initial_sites_by_vacancy_id=np.asarray(initial_sites_by_id, dtype=np.int64),
            initial_labels=np.asarray(initial_labels, dtype=np.int64),
            final_labels=np.asarray(simulator.labels, dtype=np.int64),
            **device_arrays,
        )
        catalog.write(directory / "catalog.json")
        metadata = catalog.data.get("metadata", {})
        summary = {
            "schema_version": "dftkmc-run-1", "created_utc": datetime.now(timezone.utc).isoformat(),
            "calibration_status": metadata.get("calibration_status"),
            "catalog_sha256": catalog_digest(catalog), "structure_id": catalog.structure_id,
            "catalog_provenance": metadata.get("provenance"),
            "structure_provenance": catalog.structure.get("provenance"),
            **coverage,
            "charge_model": catalog.charge_model(config.charge),
            "config": _config_dict(config),
            "effective_field_V_nm": simulator.field_V_nm.tolist(),
            "field_source": "electrode_voltage_linear_dielectric" if simulator.device is not None else "prescribed_uniform_field",
            "run_controls": {"n_vacancies": n_actual, "duration_s": duration_s, "max_events": max_events, "sample_every": sample_every},
            "termination": {"status": status, "target_reached": complete, "censored": not complete,
                            "actual_time_s": final_time, "requested_time_s": duration_s,
                            "event_count": int(simulator.event_count)},
            "initial_vacancy_indices": initial_indices,
            "initial_sites_by_vacancy_id": initial_sites_by_id,
            "energy_reference": "Catalogue energy zero is specific to this charge branch; do not compare absolute energies across charges.",
            "final": final, "event_counts": dict(simulator.counts), "channel_counts": dict(simulator.channel_counts),
            "diffusion": {"available": False, "reason": "One trajectory does not provide an independent-replica uncertainty estimate. Use study at a common physical time in zero-field bulk."},
            "assumptions": [
                "Fixed vacancy charge branch; isothermal; prescribed spatially uniform field; vacancy number conserved.",
                "Classical energies and equal forward/reverse channel attempt prefactors; vibrational free-energy effects are not included.",
                "Pair/environment closure, when enabled, is the declared catalogue approximation, not environment-resolved NEB data.",
                "Raw trajectory MSD is not by itself a converged diffusion coefficient.",
            ],
            "warnings": ([] if complete else ["The run ended before the target physical time. An event cap is computational censoring, not physical equilibration."])
                        + ([coverage["coverage_warning"]] if coverage["conditional_on_truncated_graph"] else []),
            "files": ["catalog.json", "summary.json", "trace.csv", "events.csv", "states.npz"] + (["overview.png"] if plot else []),
        }
        if simulator.device is not None:
            summary["device_boundary"] = simulator.device.summary()
            summary["energy_reference"] += (
                " energy_eV remains chemical configurational energy. electrostatic_energy_relative_bottom_eV "
                "is q times the prescribed potential relative to the bottom plane, counted once; "
                "total_energy_relative_bottom_eV is their sum at fixed electrode bias."
            )
            summary["assumptions"].extend([
                "Two ideal Dirichlet electrode planes; periodic along the first two lattice vectors and blocking oxygen at both normal boundaries.",
                "Homogeneous linear dielectric background without space-charge or electronic transport feedback; charged vacancies are test charges in the prescribed potential.",
                "Contact masks label oxide-side regions only; contact-region sites are not pinned to the metal potential.",
                "Bulk catalogue channels and site energies are used on a truncated slab; electrode-specific termination, surface energetics and interface chemistry are uncalibrated.",
                "A finite closed slab at fixed voltage has a conservative ionic potential; vacancy accumulation is not sustained ionic current or conductive filament formation.",
            ])
            summary["files"] += ["device_boundary.json", "device_sites.csv"] + (["device_boundary.png"] if plot else [])
        _write_json(directory / "summary.json", summary)
        if plot:
            _plot_overview(records, directory, metadata["calibration_status"])
            if simulator.device is not None:
                _plot_device_boundary(simulator, trajectory, directory, metadata["calibration_status"])
    return replica


def run_study(
    catalog: Catalog, config: SimulationConfig, *, n_vacancies: int = 1,
    duration_s: float, replicas: int = 8, charges: Sequence[int] | None = None,
    ablations: Sequence[str] = ("noninteracting", "pair-only", "pair-environment"),
    max_events: int = 100000, sample_every: int = 100,
    output: str | Path, confidence: float = 0.95, bootstrap_samples: int = 2000,
    plot: bool = False, allow_illustrative: bool = False,
) -> dict[str, Any]:
    """Run independent charge x environment controls at one physical time."""
    if replicas < 2:
        raise ValueError("A study requires at least two independent replicas per group.")
    if not np.isfinite(duration_s) or duration_s <= 0:
        raise ValueError("duration_s must be finite and positive.")
    if not 0 < confidence < 1 or bootstrap_samples < 100:
        raise ValueError("Use 0 < confidence < 1 and at least 100 bootstrap samples.")
    selected_charges = tuple(int(q) for q in (charges if charges is not None else catalog.charges))
    if not selected_charges or len(set(selected_charges)) != len(selected_charges):
        raise ValueError("Specify a nonempty set of unique charges.")
    if any(q not in catalog.charges for q in selected_charges):
        raise ValueError("Every selected charge must be present in the catalogue.")
    if not ablations or len(set(ablations)) != len(ablations) or any(name not in ABLATIONS for name in ablations):
        raise ValueError(f"Select distinct ablations from {tuple(ABLATIONS)}.")
    catalog.validate(allow_illustrative=allow_illustrative)
    coverage = _graph_coverage(catalog)
    device = (Simulator(catalog, replace(config, charge=selected_charges[0])).device
              if config.device is not None else None)
    effective_field = device.field_V_nm if device is not None else np.asarray(config.field_V_nm)
    directory = _new_output(output)
    catalog.write(directory / "catalog.json")
    group_count = len(selected_charges) * len(ablations)
    spawned = np.random.SeedSequence(int(config.seed)).spawn(group_count * (replicas + 1))
    seeds = [int(child.generate_state(1, dtype=np.uint64)[0]) for child in spawned]
    if len(set(seeds)) != len(seeds):
        raise RuntimeError("Unexpected seed collision; choose a different base seed.")
    seed_index = 0
    group_rows, replica_rows, groups = [], [], []
    raw_groups: dict[tuple[int, str], list[dict[str, Any]]] = {}
    for charge in selected_charges:
        for ablation in ablations:
            pair_multiplier, environment_multiplier = ABLATIONS[ablation]
            group_name = f"q{charge}_{ablation}"
            runs = []
            group_config = replace(config, charge=charge,
                                   pair_scale=float(config.pair_scale) * pair_multiplier,
                                   environment_scale=float(config.environment_scale) * environment_multiplier)
            for replica_index in range(replicas):
                seed = seeds[seed_index]
                seed_index += 1
                run_config = replace(group_config, seed=seed)
                run_dir = directory / group_name / f"replica_{replica_index:03d}"
                try:
                    record = run_simulation(catalog, run_config, n_vacancies=n_vacancies,
                                            duration_s=duration_s, max_events=max_events,
                                            sample_every=sample_every, output=run_dir, plot=plot,
                                            allow_illustrative=allow_illustrative)
                except (ValueError, RuntimeError) as exc:
                    # Preserve the attempted seed and explicit failure. Never turn
                    # a rejected model domain into a silently adjusted trajectory.
                    run_dir.mkdir(parents=True, exist_ok=True)
                    record = {"seed": seed, "status": "failed", "time_s": None,
                              "target_reached": False, "censored": True, "error": str(exc),
                              "final": {}, "events": None, "output": str(run_dir)}
                    _write_json(run_dir / "failure.json", {"config": _config_dict(run_config),
                                                           "catalog_sha256": catalog_digest(catalog),
                                                           "error": str(exc), "status": "failed"})
                record["replica_index"] = replica_index
                runs.append(record)
                replica_rows.append({"group": group_name, "charge": charge, "ablation": ablation,
                                     "replica_index": replica_index, "seed": seed,
                                     "status": record["status"], "target_reached": record["target_reached"],
                                     "time_s": record["time_s"], "requested_time_s": duration_s,
                                     "events": record["events"], "error": record.get("error", ""),
                                     "output": str(run_dir.relative_to(directory)),
                                     **{key: record.get("final", {}).get(key) for key in ("energy_eV", "msd_A2", "pair_count", "largest_cluster", "mean_displacement_A")}})
            bootstrap_seed = seeds[seed_index]
            seed_index += 1
            stats = summarize_replicas(runs, duration_s=duration_s, boundary=config.boundary,
                                       field_V_nm=effective_field, confidence=confidence,
                                       bootstrap_samples=bootstrap_samples, bootstrap_seed=bootstrap_seed)
            if coverage["conditional_on_truncated_graph"]:
                stats["warnings"].append(coverage["coverage_warning"])
                stats["diffusion"].update(
                    label="finite-time tracer estimate conditional on the declared truncated graph",
                    **coverage,
                )
                stats["diffusion"].setdefault("limitations", []).append(coverage["coverage_warning"])
            group = {"group": group_name, "charge": charge, "ablation": ablation,
                     "pair_scale": group_config.pair_scale, "environment_scale": group_config.environment_scale,
                     "effective_field_V_nm": effective_field.tolist(),
                     "statistics": stats, **coverage}
            groups.append(group)
            raw_groups[(charge, ablation)] = runs
            row: dict[str, Any] = {"group": group_name, "charge": charge, "ablation": ablation,
                                   "complete_replicas": stats["n_complete_replicas"],
                                   "incomplete_replicas": stats["n_incomplete_replicas"],
                                   "all_replicas_complete": stats["all_replicas_complete"], **coverage}
            for metric, value in stats["scalar_estimates"].items():
                row[f"{metric}_mean"] = value["mean"]
                row[f"{metric}_ci"] = value["ci"]
            diff = stats["diffusion"]
            row["finite_time_tracer_D_trace_over_3_m2_s"] = diff.get("trace_over_3_m2_s")
            row["finite_time_tracer_D_ci_m2_s"] = diff.get("trace_over_3_ci_m2_s")
            row["diffusion_label"] = diff.get("label", "")
            row["diffusion_unavailable_reason"] = diff.get("reason", "")
            group_rows.append(row)
            _write_csv(directory / "replicas.csv", replica_rows)
            _write_csv(directory / "groups.csv", group_rows)
            _write_json(directory / "study.progress.json", {"groups_completed": len(groups), "groups_total": group_count,
                                                            "groups": groups})
    contrasts = []
    if "noninteracting" in ablations:
        for charge in selected_charges:
            for ablation in ablations:
                if ablation != "noninteracting":
                    contrasts.append({"first": f"q{charge}_noninteracting", "second": f"q{charge}_{ablation}",
                                      "statistics": compare_groups(raw_groups[(charge, "noninteracting")], raw_groups[(charge, ablation)],
                                                                   duration_s=duration_s, confidence=confidence)})
    if "pair-only" in ablations and "pair-environment" in ablations:
        for charge in selected_charges:
            contrasts.append({"first": f"q{charge}_pair-only", "second": f"q{charge}_pair-environment",
                              "purpose": "Isolate the additional symmetric environment saddle correction at fixed pair-energy scale.",
                              "statistics": compare_groups(raw_groups[(charge, "pair-only")], raw_groups[(charge, "pair-environment")],
                                                           duration_s=duration_s, confidence=confidence)})
    if len(selected_charges) > 1:
        first_charge = selected_charges[0]
        for charge in selected_charges[1:]:
            for ablation in ablations:
                contrasts.append({"first": f"q{first_charge}_{ablation}", "second": f"q{charge}_{ablation}",
                                  "statistics": compare_groups(raw_groups[(first_charge, ablation)], raw_groups[(charge, ablation)],
                                                               duration_s=duration_s, confidence=confidence,
                                                               metrics=("msd_A2", "pair_count", "largest_cluster", "events")),
                                  "energy_comparison": "Unavailable: different charge branches have independent arbitrary energy zeros; no cross-charge stability inference."})
    warnings = [coverage["coverage_warning"]] if coverage["conditional_on_truncated_graph"] else []
    if n_vacancies == 1:
        warnings.append("One vacancy cannot identify vacancy-vacancy pair or environment effects; these ablations should agree within independent-replica noise.")
    if any(not group["statistics"]["all_replicas_complete"] for group in groups):
        warnings.append("At least one group is incomplete. Do not compare event-cap endpoints as equal-time physics; rerun using a larger cap or a shorter common duration.")
    result = {
        "schema_version": "dftkmc-study-1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "calibration_status": catalog.data["metadata"]["calibration_status"],
        "catalog_sha256": catalog_digest(catalog), "structure_id": catalog.structure_id,
        "catalog_provenance": catalog.data["metadata"].get("provenance"),
        "structure_provenance": catalog.structure.get("provenance"),
        **coverage,
        "initial_placement": "Uniform sample of vacancy sites without replacement, independently seeded for each replica.",
        "config": _config_dict(config), "duration_s": duration_s, "n_vacancies": n_vacancies,
        "effective_field_V_nm": effective_field.tolist(),
        "replicas_per_group": replicas, "base_seed": int(config.seed),
        "seed_scheme": "NumPy SeedSequence child streams; every trajectory and every bootstrap has a distinct seed.",
        "study_controls": {"max_events": max_events, "sample_every": sample_every,
                           "confidence": confidence, "bootstrap_samples": bootstrap_samples},
        "groups": groups, "contrasts": contrasts, "replica_outcomes": replica_rows, "warnings": warnings,
        "energy_reference": "Energy zeros are charge-branch specific. Energy contrasts are reported only within a fixed charge branch.",
        "interpretation": [
            "All comparisons use one requested physical time, not equal event counts.",
            "The uncertainty unit is a whole independent replica; time samples and vacancies are not independent replicas.",
            "Finite-time tracer estimates do not establish an Arrhenius law, a converged long-time coefficient, or experimental device behavior.",
            "Ablation scales modify only the catalogue's explicitly declared pair/environment closure; they do not create new DFT data.",
        ],
    }
    if device is not None:
        result["device_boundary"] = device.summary()
        result["interpretation"].append(
            "Finite blocking device slab with electrode-prescribed linear potential; no bulk diffusion inference, "
            "electronic current, interface oxygen exchange or self-consistent space-charge solution."
        )
    _write_json(directory / "study.json", result)
    (directory / "study.progress.json").unlink(missing_ok=True)
    return result


def _write_poscar(path: Path, cell: np.ndarray, hf: np.ndarray, oxygen: np.ndarray, title: str) -> None:
    lines = [title, "1.0"]
    lines += ["  " + " ".join(f"{x:.14f}" for x in row) for row in cell]
    lines += ["Hf O", f"{len(hf)} {len(oxygen)}", "Cartesian"]
    lines += ["  " + " ".join(f"{x:.14f}" for x in row) for row in np.vstack((hf, oxygen))]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_neb(
    lattice: Lattice, *, output: str | Path, charges: Sequence[int] = (0, 2),
    channel_ids: Sequence[str] | None = None,
    source_catalog: Catalog | None = None,
) -> dict[str, Any]:
    """Export *unrelaxed* ordered endpoints and a pending calculation worklist.

    The vacancy moves source -> target, so the target oxygen moves by minus
    the actual hop displacement.  Wrapped target-source coordinates must not
    be used to construct a path crossing a periodic boundary.
    """
    if lattice.pbc is not True and not np.all(lattice.pbc):
        raise ValueError("NEB endpoint export requires a bulk periodic lattice.")
    if source_catalog is not None and source_catalog.structure_id != lattice.structure_id:
        raise ValueError("Source catalogue and NEB lattice must have the same structure identity.")
    source_coverage = _graph_coverage(source_catalog) if source_catalog is not None else None
    exclusion_handling = (
        "Source catalogue channel exclusions are preserved as provenance only and were NOT automatically applied. "
        "The pending catalogue regenerates all geometric candidate channels within the export cutoff; "
        "an explicit channel selection restricts endpoint tasks only."
    )
    source_provenance = {
        "source_catalog_sha256": catalog_digest(source_catalog) if source_catalog is not None else None,
        "source_catalog_provenance": source_catalog.data.get("metadata", {}).get("provenance") if source_catalog is not None else None,
        "source_migration_cutoff_A": source_coverage["migration_cutoff_A"] if source_coverage else None,
        "source_channel_exclusions": source_coverage["channel_exclusions"] if source_coverage else {},
        "source_channel_exclusions_applied": False,
        "source_channel_exclusions_handling": exclusion_handling,
    }
    qs = tuple(int(q) for q in charges)
    if not qs or len(set(qs)) != len(qs) or any(float(q) != int(q) for q in charges):
        raise ValueError("NEB charges must be distinct integers.")
    templates = {item["id"]: item for item in lattice.channel_templates()}
    selected = list(channel_ids) if channel_ids is not None else list(templates)
    if len(set(selected)) != len(selected) or not selected or any(key not in templates for key in selected):
        raise ValueError("Select distinct known channel IDs.")
    representatives = {}
    for edge in lattice.migration_edges:
        representatives.setdefault(edge.channel_id, edge)
    missing = [key for key in selected if key not in representatives]
    if missing:
        raise ValueError(f"No representative edge in this supercell for channels: {missing}")
    directory = _new_output(output)
    pending = make_catalog(lattice, status="pending", charges=qs)
    pending.data["metadata"].update(source_provenance)
    pending.write(directory / "catalog.pending.json")
    tasks = []
    for channel_id in selected:
        edge = representatives[channel_id]
        source, target = int(edge.source), int(edge.target)
        if source == target:
            raise ValueError("A NEB path cannot use a periodic self-edge; enlarge the supercell.")
        displacement = np.asarray(edge.displacement_A, dtype=float)
        oxygen_indices = [index for index in range(lattice.n) if index != source]
        moving_row = oxygen_indices.index(target)
        initial_oxygen = np.asarray(lattice.oxygen_A[oxygen_indices], dtype=float).copy()
        final_oxygen = initial_oxygen.copy()
        final_oxygen[moving_row] = initial_oxygen[moving_row] - displacement
        atoms = ([{"poscar_index_1based": index + 1, "species": "Hf", "lattice_index": index}
                  for index in range(len(lattice.hf_A))] +
                 [{"poscar_index_1based": len(lattice.hf_A) + row + 1, "species": "O", "lattice_index": index}
                  for row, index in enumerate(oxygen_indices)])
        moving_poscar = len(lattice.hf_A) + moving_row + 1
        for charge in qs:
            task_dir = directory / str(channel_id) / f"q{charge}"
            task_dir.mkdir(parents=True)
            title = f"UNRELAXED dftkmc {channel_id} q={charge}; one oxygen vacancy"
            _write_poscar(task_dir / "POSCAR.initial", np.asarray(lattice.cell_A), np.asarray(lattice.hf_A), initial_oxygen, title)
            _write_poscar(task_dir / "POSCAR.final", np.asarray(lattice.cell_A), np.asarray(lattice.hf_A), final_oxygen, title)
            metadata = {
                "status": "unrelaxed_endpoints_require_charge_consistent_DFT_relaxation",
                "structure_id": lattice.structure_id, "channel_id": channel_id,
                "template": templates[channel_id], "charge": charge,
                "source_channel_exclusion_reason": source_provenance["source_channel_exclusions"].get(channel_id),
                "source_channel_exclusions_applied": False,
                "charge_definition": "Net charge of the entire defect supercell containing one vacancy; not a per-atom charge assignment.",
                "vacancy_source_site": source, "vacancy_target_site": target,
                "vacancy_displacement_A": displacement.tolist(),
                "oxygen_displacement_A": (-displacement).tolist(),
                "moving_oxygen_poscar_index_1based": moving_poscar,
                "moving_oxygen_lattice_index": target,
                "moving_oxygen_initial_A": initial_oxygen[moving_row].tolist(),
                "moving_oxygen_final_A": final_oxygen[moving_row].tolist(),
                "vacancy_source_wrapped_A": np.asarray(lattice.oxygen_A[source]).tolist(),
                "vacancy_target_wrapped_A": np.asarray(lattice.oxygen_A[target]).tolist(),
                "periodic_image_rule": "Final moving oxygen = initial target oxygen minus the actual canonical vacancy hop displacement; Cartesian coordinates can lie outside the cell.",
                "atom_identity_order": atoms,
                "atom_counts": {"Hf": len(lattice.hf_A), "O": len(oxygen_indices)},
                "electron_count_rule": "NELECT = sum(species_count * POTCAR_ZVAL) - charge; positive charge removes electrons. Use identical atom count, POTCAR and NELECT along the path.",
                "relaxation_required": True,
                "calculation_requirements": [
                    "Select and converge a consistent exchange-correlation functional, pseudopotentials, k mesh, cutoffs, charge correction and spin treatment.",
                    "Relax both endpoints at fixed cell with identical charge and preserve/reconcile atom identity before interpolating NEB images.",
                    "Verify both relaxed endpoints remain the intended metastable vacancy sites; a geometric hop is only a candidate transition.",
                    "Keep the chosen periodic image of the moving oxygen continuous when constructing images; do not blindly wrap all images independently.",
                    "Record the initial, final and transition-state energies with a common energy reference; add barriers and prefactors only after completed calculations.",
                ],
            }
            _write_json(task_dir / "mapping.json", metadata)
            tasks.append({"channel_id": channel_id, "charge": charge, "structure_id": lattice.structure_id,
                          "source_site": source, "target_site": target,
                          "from_basis": templates[channel_id]["from_basis"], "to_basis": templates[channel_id]["to_basis"],
                          "cell_shift": templates[channel_id]["cell_shift"],
                          "vacancy_displacement_A": displacement.tolist(), "length_A": float(np.linalg.norm(displacement)),
                          "moving_oxygen_poscar_index_1based": moving_poscar,
                          "directory": str(task_dir.relative_to(directory)), "status": "pending_unrelaxed",
                          "source_channel_exclusion_reason": source_provenance["source_channel_exclusions"].get(channel_id, ""),
                          "source_channel_exclusions_applied": False,
                          "initial_energy_eV": None, "final_energy_eV": None, "saddle_eV": None,
                          "attempt_Hz": None, "provenance": ""})
    _write_csv(directory / "tasks.csv", tasks)
    summary = {"schema_version": "dftkmc-neb-worklist-1", "structure_id": lattice.structure_id,
               "migration_cutoff_A": float(lattice.migration_cutoff_A), **source_provenance,
               "n_channels": len(selected), "charges": list(qs), "n_tasks": len(tasks),
               "status": "pending_unrelaxed", "relaxation_required": True,
               "symmetry_reduction": "Translation and path reversal equivalence only; no full space-group reduction is claimed.",
               "files": ["tasks.csv", "catalog.pending.json", "README.txt", "worklist.json"],
               "no_jobs_submitted": True, "potcar_included": False}
    _write_json(directory / "worklist.json", summary)
    (directory / "README.txt").write_text(
        "DFT/NEB WORKLIST — UNRELAXED ENDPOINTS\n\n"
        "Each task contains POSCAR.initial, POSCAR.final and mapping.json. Atom identities and ordering are identical in both endpoint files. The target oxygen fills the source vacancy, moving by minus the actual periodic vacancy-hop displacement. Cartesian positions outside the cell preserve the intended path image.\n\n"
        "These are candidate structures, not completed DFT or NEB calculations. Relax both endpoints at the same fixed cell and charge, verify metastability and atom mapping, then construct and converge NEB images. The conventional-cell geometry and supercell size also need convergence/relaxation assessment.\n\n"
        "Set NELECT from the valence electron count of your own POTCAR minus q. Positive q removes electrons. Use identical POTCAR, atom count, charge, spin treatment and numerical settings throughout each path. No POTCAR is distributed and no calculation has been submitted.\n\n"
        "catalog.pending.json deliberately contains missing numerical data and cannot run. Fill it only with completed, traceable calculations on the same structure and charge branches. Site and saddle energies require a common reference. An activation barrier alone is not an absolute saddle energy; saddle = initial site energy + forward barrier, and verify the reverse relation. Attempt prefactors require a declared calculation or approximation.\n\n"
        "The channel list removes translation and reverse duplicates only; it is not a complete space-group symmetry reduction. A requested subset exports only those tasks, while the pending catalogue retains the complete channel worklist so uncovered rates remain blocked.\n\n"
        + exclusion_handling + " Source IDs/reasons and catalogue identity are retained in worklist.json and catalog.pending.json metadata.\n",
        encoding="utf-8")
    return summary
