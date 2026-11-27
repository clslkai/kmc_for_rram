"""Plots and a dependency-free offline viewer of saved coupled KMC states.

Only recorded states and electronic graph displacements are drawn.  This module
does not run KMC, interpolate defect motion, infer missing paths, or change any
scientific output.  Periodic bonds are split at A/B cell faces; C is never
wrapped.  Electronic-network results remain explicitly illustrative/isothermal.
"""
from __future__ import annotations

import base64
import csv
import json
from pathlib import Path

import numpy as np


_NUMERIC_TRACE = (
    "time_s", "events", "vacancies", "oxide_oxygen", "bottom_reservoir_oxygen",
    "top_reservoir_oxygen", "total_oxygen", "mass_residual", "reaction_charge_residual_e",
    "equivalent_conductance_S", "terminal_current_A", "requested_voltage_V",
    "applied_voltage_V", "read_current_A",
)


def _boolean(value):
    if value in (True, "True", "true", "1"):
        return True
    if value in (False, "False", "false", "0"):
        return False
    raise ValueError(f"unrecognized saved boolean: {value!r}")


def _load(output):
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    with (output / "trace.csv").open(newline="", encoding="utf-8") as handle:
        raw_records = list(csv.DictReader(handle))
    records = []
    for raw in raw_records:
        record = {name: float(raw[name]) for name in _NUMERIC_TRACE}
        for name in ("vacancy_graph_spanning", "model_formation_criterion"):
            record[name] = _boolean(raw[name])
        record["conducting_path_site_ids"] = json.loads(raw["conducting_path_site_ids"])
        records.append(record)
    if not records:
        raise ValueError("coupled plot requires at least one saved snapshot")
    with np.load(output / "states.npz", allow_pickle=False) as archive:
        states = {name: archive[name].copy() for name in archive.files}
    frames = len(records)
    if states["vacancy"].ndim != 2:
        raise ValueError("states vacancy must have shape (saved frames, sites)")
    count = states["vacancy"].shape[1]
    edges = len(states["electronic_pairs"])
    required_shapes = {
        "time_s": (frames,), "events": (frames,), "vacancy": (frames, count),
        "oxygen_ids": (frames, count), "vacancy_ids": (frames, count),
        "potential_V": (frames, count), "oxygen_A": (count, 3), "cell_A": (3, 3),
        "depth_A": (count,), "bottom_contact": (count,), "top_contact": (count,),
        "electronic_pairs": (edges, 2), "electronic_displacements_A": (edges, 3),
        "edge_current_A": (frames, edges), "edge_field_V_nm": (frames, edges),
    }
    for name, shape in required_shapes.items():
        if states[name].shape != shape or not np.isfinite(states[name]).all():
            raise ValueError(f"saved {name} must be finite with shape {shape}")
    if tuple(states["pbc"].tolist()) != (True, True, False):
        raise ValueError("coupled viewer supports only A/B-periodic, C-finite slabs")
    times = np.array([record["time_s"] for record in records])
    if (not np.array_equal(times, states["time_s"])
            or not np.array_equal([record["events"] for record in records], states["events"])):
        raise ValueError("trace.csv and states.npz must align one-to-one; no frame interpolation is performed")
    if np.any(np.diff(times) < 0):
        raise ValueError("saved frame times must be nondecreasing (pulse-boundary duplicates are allowed)")
    if not all(np.isfinite(record[name]) for record in records for name in _NUMERIC_TRACE):
        raise ValueError("saved trace values must be finite")
    pair_index = {tuple(sorted(map(int, pair))): index for index, pair in enumerate(states["electronic_pairs"])}
    path_edges = []
    for frame, record in enumerate(records):
        if record["vacancies"] != int(np.count_nonzero(states["vacancy"][frame])):
            raise ValueError("saved vacancy counts disagree between trace.csv and states.npz")
        path = record["conducting_path_site_ids"]
        if (not isinstance(path, list) or any(isinstance(i, bool) or not isinstance(i, int)
                                             or not 0 <= i < count for i in path)):
            raise ValueError("saved conducting path must contain valid integer site IDs")
        if bool(path) != record["vacancy_graph_spanning"]:
            raise ValueError("saved path and vacancy_graph_spanning flag disagree")
        if path and (not states["bottom_contact"][path[0]] or not states["top_contact"][path[-1]]
                     or not all(states["vacancy"][frame, i] for i in path)):
            raise ValueError("saved path must join the contact regions through recorded vacancy sites")
        try:
            path_edges.append([pair_index[tuple(sorted((i, j)))] for i, j in zip(path[:-1], path[1:])])
        except KeyError as exc:
            raise ValueError("saved conducting path contains a non-electronic edge") from exc
    segments = _periodic_edge_segments(states["oxygen_A"], states["cell_A"],
                                       states["electronic_pairs"], states["electronic_displacements_A"])
    return summary, records, states, segments, path_edges


def _periodic_edge_segments(positions, cell, pairs, displacements):
    """Split real short edge vectors into pieces within the A/B cell faces."""
    inverse = np.linalg.inv(cell)
    result = []
    for (source, target), displacement in zip(pairs, displacements):
        start = positions[source] @ inverse
        step = displacement @ inverse
        endpoint_error = (positions[source] + displacement - positions[target]) @ inverse
        if abs(endpoint_error[2]) > 1e-8 or np.max(np.abs(endpoint_error[:2] - np.rint(endpoint_error[:2]))) > 1e-8:
            raise ValueError("electronic displacement is inconsistent with lateral-only periodic endpoints")
        cuts = [0., 1.]
        for axis in (0, 1):
            if abs(step[axis]) < 1e-14:
                continue
            lo, hi = sorted((start[axis], start[axis] + step[axis]))
            for face in range(int(np.floor(lo)) + 1, int(np.ceil(hi))):
                t = (face - start[axis]) / step[axis]
                if 1e-12 < t < 1 - 1e-12:
                    cuts.append(float(t))
        cuts = sorted(set(cuts))
        pieces = []
        for left, right in zip(cuts[:-1], cuts[1:]):
            middle = start + (left + right) / 2 * step
            image = np.array([np.floor(middle[0]), np.floor(middle[1]), 0.])
            pieces.append(np.stack(((start + left * step - image) @ cell,
                                    (start + right * step - image) @ cell)))
        result.append(pieces)
    return result


def _time_scale(times):
    maximum = float(np.max(np.abs(times), initial=0.0))
    for factor, label in ((1., "s"), (1e-3, "ms"), (1e-6, "µs"), (1e-9, "ns"), (1e-12, "ps")):
        if maximum >= factor or factor == 1e-12:
            return factor, label
    return 1., "s"


def _current_scale(currents):
    maximum = float(np.max(np.abs(currents), initial=0.0))
    for factor, label in ((1., "A"), (1e-3, "mA"), (1e-6, "µA"), (1e-9, "nA"), (1e-12, "pA")):
        if maximum >= factor or factor == 1e-12:
            return factor, label
    return 1., "A"


def _figure(output, summary, records, states, segments, path_edges):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    cell, positions = states["cell_A"], states["oxygen_A"]
    normal = np.cross(cell[0], cell[1]); normal /= np.linalg.norm(normal)
    horizontal = cell[0] / np.linalg.norm(cell[0])
    projection = np.stack((horizontal, normal), axis=1) / 10.
    xy = positions @ projection
    thickness_nm = float(cell[2] @ normal) / 10.
    corners = np.array([a * cell[0] + b * cell[1] + c * cell[2]
                        for a in (0, 1) for b in (0, 1) for c in (0, 1)]) @ projection
    potential_min, potential_max = float(states["potential_V"].min()), float(states["potential_V"].max())
    if potential_min == potential_max:
        color_norm = Normalize(potential_min - .5, potential_max + .5)
    else:
        color_norm = Normalize(potential_min, potential_max)
    times = states["time_s"]
    time_factor, time_unit = _time_scale(times)
    t = times / time_factor
    column = lambda key: np.asarray([record[key] for record in records])
    conductance = column("equivalent_conductance_S")
    criterion_frames = np.flatnonzero(column("model_formation_criterion").astype(bool))
    selected_frame = int(criterion_frames[0]) if len(criterion_frames) else len(records) - 1
    selected_title = "First joint-criterion state" if len(criterion_frames) else "Final recorded state"
    current = column("terminal_current_A")
    current_factor, current_unit = _current_scale(current)
    cfg = summary["config"]
    termination = summary.get("termination", {})
    status = termination.get("status", "saved run") if isinstance(termination, dict) else str(termination)
    with plt.rc_context({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
                         "legend.fontsize": 8.5, "xtick.labelsize": 9, "ytick.labelsize": 9,
                         "axes.spines.top": False, "axes.spines.right": False}):
        fig, axes = plt.subplots(3, 2, figsize=(13.4, 13.4), layout="constrained")
        fig.set_layout_engine("constrained", rect=(0., .075, 1., .925),
                              w_pad=.12, h_pad=.12, hspace=.06, wspace=.06)
        fig.suptitle("Coupled oxygen exchange + effective electronic network\n"
                     f"ILLUSTRATIVE · q = {cfg['charge']:+d}e · {cfg['temperature_K']:g} K · "
                     f"{len(records)} saved states · {status.replace('_', ' ')}", fontsize=15, fontweight="bold")
        for ax, frame, letter, title in ((axes[0, 0], 0, "A", "Initial recorded configuration"),
                                          (axes[0, 1], selected_frame, "B", selected_title)):
            vacancy = states["vacancy"][frame]
            scatter = ax.scatter(xy[:, 0], xy[:, 1], c=states["potential_V"][frame], norm=color_norm,
                                 cmap="viridis", s=15, alpha=.45, edgecolors="none", zorder=2)
            ax.scatter(xy[vacancy, 0], xy[vacancy, 1], c=states["potential_V"][frame, vacancy],
                       norm=color_norm, cmap="viridis", s=57, edgecolors="#e67320", linewidths=1.5, zorder=4)
            pieces = [piece @ projection for edge in path_edges[frame] for piece in segments[edge]]
            if pieces:
                ax.add_collection(LineCollection(pieces, colors="#d74032", linewidths=2.5, zorder=3))
            for face, offset, label in ((0, np.zeros(3), "Bottom electrode"), (1, cell[2], "Top electrode")):
                endpoints = np.array([offset, offset + cell[0], offset + cell[1], offset + cell[0] + cell[1]]) @ projection
                ax.hlines(face * thickness_nm, endpoints[:, 0].min(), endpoints[:, 0].max(),
                          color="#4b5563", linewidth=2.)
            ax.set_xlim(corners[:, 0].min() - .15, corners[:, 0].max() + .15)
            ax.set_ylim(-.10, thickness_nm + .10)
            ax.set_xlabel("Lateral projection along A (nm)")
            ax.set_ylabel("Depth normal to electrodes (nm)")
            ax.set_title(f"{letter}  {title}\n"
                         f"t = {times[frame] / time_factor:.4g} {time_unit}, "
                         f"Nv = {int(vacancy.sum())}, vacancy path = {'yes' if path_edges[frame] else 'no'}", loc="left")
            ax.grid(alpha=.12)
            ax.set_aspect("equal", adjustable="box")
        legend = [Line2D([], [], marker='o', linestyle='none', markerfacecolor='#7c8a76',
                         markeredgecolor='#e67320', markeredgewidth=1.5, markersize=7, label='Vacancy site'),
                  Line2D([], [], color='#d74032', linewidth=2.5, label='Recorded vacancy path')]
        axes[0, 0].legend(handles=legend, loc="upper left", framealpha=.92)
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=color_norm, cmap="viridis"),
                               ax=list(axes[0]), fraction=.035, pad=.015, shrink=.83)
        colorbar.set_label("Solved site potential (V)")
        if potential_min == potential_max:
            colorbar.set_ticks([potential_min])

        ax = axes[1, 0]
        positive = conductance > 0
        if np.any(positive):
            ax.plot(t, np.where(positive, conductance, np.nan), color="#1c7890", linewidth=1.7,
                    marker='.', markersize=3, label="Frozen-state conductance")
            ax.set_yscale('log')
        else:
            ax.plot(t, conductance, color="#1c7890", linewidth=1.7, label="Frozen-state conductance")
            ax.text(.04, .8, "No terminal-connected electronic network: G = 0", transform=ax.transAxes, fontsize=9)
        threshold = max(float(cfg["formation_min_conductance_S"]),
                        float(cfg["formation_min_gain"]) * conductance[0])
        if np.any(positive):
            ax.axhline(threshold, color="#a86419", linestyle='--', linewidth=1., label="Declared G criterion")
        met = column("model_formation_criterion").astype(bool)
        if np.any(met & positive):
            ax.scatter(t[met & positive], conductance[met & positive], s=22, facecolors='none',
                       edgecolors='#ad3a24', linewidths=1., label="Full model criterion met", zorder=4)
        if np.any(~positive) and np.any(positive):
            ax.text(.03, .03, f"{int(np.count_nonzero(~positive))} zero-G states omitted on log axis", transform=ax.transAxes, fontsize=8)
        ax.set_title("C  Electrical evolution of saved defect states", loc="left")
        ax.set_ylabel("Equivalent conductance G (S)")
        ax.legend(loc="best", framealpha=.9)

        ax = axes[1, 1]
        ax.plot(t, current / current_factor, color="#8050a3", linewidth=1.8, label="Top-terminal current")
        ax.set_ylabel(f"Signed terminal current ({current_unit})", color="#8050a3")
        ax.tick_params(axis='y', colors="#8050a3")
        bias_ax = ax.twinx()
        bias_ax.spines['right'].set_visible(True)
        bias_ax.plot(t, column("requested_voltage_V"), color="#9da5ad", linestyle='--', linewidth=1.2, label="Requested bias")
        bias_ax.plot(t, column("applied_voltage_V"), color="#264c65", linewidth=1.4, label="Actual bias")
        bias_ax.set_ylabel("Top − bottom voltage (V)", color="#264c65")
        lines = ax.get_lines() + bias_ax.get_lines()
        ax.legend(lines, [line.get_label() for line in lines], loc='best', framealpha=.9)
        ax.set_title("D  Write current and compliance-limited voltage", loc="left")

        ax = axes[2, 0]
        for key, color, label in (("vacancies", "#dd7d27", "Oxide vacancies"),
                                  ("bottom_reservoir_oxygen", "#276b9a", "Bottom reservoir O"),
                                  ("top_reservoir_oxygen", "#31816a", "Top reservoir O")):
            ax.plot(t, column(key), linewidth=2.1 if key == "vacancies" else 1.4,
                    linestyle='--' if key == "top_reservoir_oxygen" else '-', color=color, label=label)
        ax.set_title("E  Defect population and finite oxygen reservoirs", loc="left")
        ax.set_ylabel("Count (sites or O atoms)")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc="best", framealpha=.9)
        ax.text(.02, .03, f"Oxide O = {len(positions)} − Nv; total O = "
                f"{int(records[0]['total_oxygen'])} → {int(records[-1]['total_oxygen'])}",
                transform=ax.transAxes, fontsize=8.5, bbox=dict(facecolor='white', alpha=.8, edgecolor='none'))

        ax = axes[2, 1]
        mass, charge = column("mass_residual"), column("reaction_charge_residual_e")
        ax.plot(t, mass, color='#237d9a', linewidth=2., label='Oxygen balance (O atoms)')
        ax.plot(t, charge, color='#c27031', linestyle='--', linewidth=1.6, label='Reaction charge balance (e)')
        if max(np.max(np.abs(mass)), np.max(np.abs(charge))) == 0:
            ax.set_ylim(-1., 1.)
        ax.set_title("F  Recorded oxygen and reaction-charge ledgers", loc="left")
        ax.set_ylabel("Residual (O atoms / e)")
        ax.legend(loc='upper left', framealpha=.9)
        ax.text(.03, .13, f"Max |O residual| = {np.max(np.abs(mass)):g}\n"
                f"Max |reaction-charge residual| = {np.max(np.abs(charge)):g}",
                transform=ax.transAxes, fontsize=9)
        for ax in axes[1:].ravel():
            ax.set_xlabel(f"Time ({time_unit})")
            ax.grid(alpha=.2, which='major')
            if t[-1] > t[0]:
                ax.set_xlim(t[0], t[-1])
        fig.text(.5, .016,
                 "ILLUSTRATIVE effective-network result · Isothermal · No space-charge Poisson / quantum transport\n"
                 "A/B are 2D projections; overlap is not a bond. Paths use saved graph edges and short periodic images.\n"
                 "B shows the first saved joint-criterion state when available; final states remain in the time traces and viewer.\n"
                 "Lines connect saved snapshots. The declared model criterion is not experimental forming validation.",
                 fontsize=8.7, ha="center", va="bottom")
        fig.savefig(output / "coupled_overview.png", dpi=180, facecolor='white')
        plt.close(fig)


def _viewer(output, summary, records, states, segments, path_edges):
    cell = states["cell_A"]
    normal = np.cross(cell[0], cell[1]); normal /= np.linalg.norm(normal)
    x_axis = cell[0] / np.linalg.norm(cell[0])
    y_axis = np.cross(normal, x_axis)
    basis = np.stack((x_axis, y_axis, normal), axis=1)
    centers = (cell[0] + cell[1] + cell[2]) / 2
    corners = np.array([a * cell[0] + b * cell[1] + c * cell[2]
                        for a in (0, 1) for b in (0, 1) for c in (0, 1)])
    convert = lambda value: np.round((np.asarray(value) - centers) @ basis, 8).tolist()
    # Compact display-only data.  Full-precision NPZ/CSV files are untouched.
    # The overlay needs relative magnitudes, so normalize before float32
    # encoding to avoid loss of physically small current units (e.g. pA).
    currents = states["edge_current_A"]
    peak = np.max(np.abs(currents), axis=1, initial=0.0)
    relative_currents = np.divide(currents, peak[:, None], out=np.zeros_like(currents), where=peak[:, None] > 0)
    encoded_currents = base64.b64encode(relative_currents.astype("<f4").tobytes()).decode("ascii")
    displayed_potential = [[float(f"{value:.6g}") for value in row] for row in states["potential_V"]]
    data = {
        "title": "Saved oxygen-exchange / resistor-network trajectory",
        "summary": {key: summary.get(key) for key in ("config", "calibration_status", "initial_spanning", "termination")},
        "records": records,
        "positions": convert(states["oxygen_A"]), "hf": convert(states["hf_A"]),
        "corners": convert(corners), "segments": [[convert(piece) for piece in pieces] for pieces in segments],
        "pairs": states["electronic_pairs"].tolist(), "pathEdges": path_edges,
        "vacancy": states["vacancy"].astype(int).tolist(), "potential": displayed_potential,
        "oxygenIds": states["oxygen_ids"].tolist(), "vacancyIds": states["vacancy_ids"].tolist(),
        "edgeCurrentRelativeF32LE": encoded_currents,
        "bottom": states["bottom_contact"].astype(int).tolist(), "top": states["top_contact"].astype(int).tolist(),
        "potentialMin": float(states["potential_V"].min()), "potentialMax": float(states["potential_V"].max()),
    }
    payload = json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(',', ':')).replace('<', '\\u003c')
    (output / "viewer.html").write_text(_HTML.replace("__PAYLOAD__", payload), encoding="utf-8")


def plot_coupled_run(output: Path) -> None:
    """Create ``coupled_overview.png`` and offline ``viewer.html`` from outputs."""
    output = Path(output)
    summary, records, states, segments, path_edges = _load(output)
    _figure(output, summary, records, states, segments, path_edges)
    _viewer(output, summary, records, states, segments, path_edges)


_HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Saved KMC device trajectory</title>
<style>
:root{color-scheme:light;font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#183246;background:#edf2f5}
*{box-sizing:border-box}body{margin:0}header{padding:22px 28px 17px;background:#102d3f;color:white}h1{margin:0 0 7px;font-size:23px;font-weight:650}header p{margin:0;color:#bcd2df;font-size:13px;line-height:1.55}.badge{display:inline-block;padding:4px 8px;margin-right:8px;border-radius:5px;font-size:10px;letter-spacing:.08em;background:#f5c66a;color:#493208;font-weight:750}
main{max-width:1550px;margin:0 auto;padding:18px;display:grid;grid-template-columns:300px minmax(0,1fr);gap:17px}.card{border:1px solid #d7e1e8;border-radius:11px;background:white;overflow:hidden}.side{padding:17px}.side h2{font-size:14px;margin:0 0 13px}.kv{display:grid;grid-template-columns:1fr auto;gap:8px;font-size:12px;padding:10px 0;border-bottom:1px solid #e7edf1}.kv strong{font-variant-numeric:tabular-nums;font-weight:650;text-align:right}.small{font-size:11px;line-height:1.6;color:#637989}.caption{font-size:12px;line-height:1.6;color:#5b7080;margin:13px 0 0}.check{display:flex;align-items:center;gap:7px;font-size:12px;margin:10px 0}.status{font-size:12px;padding:9px;border-radius:6px;background:#eff6f7;margin:13px 0}.good{color:#167357}.warn{color:#9b5d19}.controls{padding:14px 16px;border-bottom:1px solid #dce6ec;display:flex;align-items:center;gap:11px;flex-wrap:wrap}button,select{font:inherit;font-size:12px;border:1px solid #c8d7df;border-radius:6px;padding:7px 10px;background:#f7fafc;color:#183246;cursor:pointer}button:hover{background:#e8f2f7}#play{background:#176580;color:white;border-color:#176580;min-width:64px}input[type=range]{flex:1;min-width:150px;accent-color:#15728e}#frameLabel{font-size:12px;font-variant-numeric:tabular-nums;min-width:85px}.canvas-wrap{position:relative;background:#f8fbfd}canvas{display:block;width:100%;height:560px;touch-action:none;cursor:grab}canvas:active{cursor:grabbing}.hint{position:absolute;top:12px;left:14px;color:#6b8190;font-size:11px;pointer-events:none}.legend{padding:12px 17px;display:flex;align-items:center;gap:12px;flex-wrap:wrap;font-size:11px;border-top:1px solid #e2eaf0}.gradient{height:10px;width:150px;border-radius:5px;background:linear-gradient(90deg,#440154,#3b528b,#21918c,#5ec962,#fde725)}.dot{display:inline-block;width:10px;height:10px;background:#80aaa3;border:2px solid #ed8d36;border-radius:50%;vertical-align:middle;margin-right:4px}.path{display:inline-block;width:21px;height:3px;background:#d74032;vertical-align:middle;margin-right:4px}.bottom-note{padding:13px 17px;background:#f6f9fb;border-top:1px solid #e2eaf0;font-size:11px;line-height:1.6;color:#536c7c}.inspector{background:#f7fafb;padding:10px;border-radius:6px;font-size:11px;line-height:1.7;margin-top:14px}.strip{padding:0 16px 13px}#spark{height:80px;cursor:default}footer{text-align:center;padding:0 22px 20px;color:#607887;font-size:11px;line-height:1.6}@media(max-width:900px){main{grid-template-columns:1fr}.side{display:grid;grid-template-columns:1fr 1fr;gap:0 20px}.side h2,.side .status,.side .inspector,.side .caption{grid-column:1/-1}canvas{height:470px}}@media(max-width:540px){header{padding:18px}h1{font-size:19px}main{padding:10px}.side{display:block}canvas{height:410px}.controls{gap:7px}}
</style></head><body>
<header><h1>Oxygen exchange and electronic-network evolution</h1><p><span class="badge">ILLUSTRATIVE MODEL</span>Actual saved configurations · finite oxygen reservoirs · classical resistor-network potentials</p></header>
<main><aside class="card side"><h2>Selected saved state</h2><div id="stats"></div><div id="criterion" class="status"></div>
<label class="check"><input type="checkbox" id="showO" checked> Show oxygen-site backdrop</label>
<label class="check"><input type="checkbox" id="showHf"> Show fixed Hf reference atoms</label>
<label class="check"><input type="checkbox" id="showPath" checked> Show recorded vacancy path</label>
<label class="check"><input type="checkbox" id="showCurrent"> Show electronic edge currents</label>
<div id="inspect" class="inspector">Click a site to inspect its saved identity and local potential.</div>
<p class="caption">Drag to rotate; scroll to zoom. The current overlay uses the saved edge currents. Color and width indicate relative magnitude; it is not a separate conduction threshold.</p>
</aside><section class="card"><div class="controls"><button id="play">Play</button><button id="prev" aria-label="Previous saved frame">◀</button><button id="next" aria-label="Next saved frame">▶</button><input id="frame" type="range" min="0" value="0" aria-label="Saved frame"><span id="frameLabel"></span><select id="speed" aria-label="Playback speed"><option value="2">2 frames/s</option><option value="5" selected>5 frames/s</option><option value="10">10 frames/s</option></select><button id="formed">First criterion</button><button id="final">Final</button><button id="reset">Reset view</button></div>
<div class="canvas-wrap"><canvas id="view" aria-label="Rotatable view of actual saved oxide coordinates"></canvas><div class="hint">Orthographic view · coordinates in Å · A/B periodic, electrode-normal direction finite</div></div>
<div class="legend"><span class="dot"></span><span>Vacancy site</span><span class="path"></span><span>Recorded vacancy-only graph path</span><span class="gradient"></span><span id="potentialLegend"></span></div>
<div class="strip"><canvas id="spark" aria-label="Saved conductance history"></canvas></div>
<div class="bottom-note">Frames are recorded states, with no interpolated ionic motion. Playback advances a fixed number of frames per second, not uniform physical time. A/B boundary-crossing edges are split using their saved short periodic displacements. Display values are rounded; full precision remains in the source files. The oxide contact sites retain solved potentials; the metal planes are separate ideal terminals.</div>
</section></main><footer id="footer"></footer>
<script>
"use strict";
const DATA=__PAYLOAD__;
const bytes=Uint8Array.from(atob(DATA.edgeCurrentRelativeF32LE),c=>c.charCodeAt(0)),binary=new DataView(bytes.buffer),relativeCurrent=new Float32Array(bytes.length/4);for(let i=0;i<relativeCurrent.length;i++)relativeCurrent[i]=binary.getFloat32(4*i,true);DATA.edgeCurrent=DATA.records.map((_,i)=>relativeCurrent.subarray(i*DATA.pairs.length,(i+1)*DATA.pairs.length));delete DATA.edgeCurrentRelativeF32LE;
const $=id=>document.getElementById(id), canvas=$('view'), ctx=canvas.getContext('2d'), spark=$('spark'), sctx=spark.getContext('2d');
const frame=$('frame'); frame.max=DATA.records.length-1;
let index=0,yaw=.65,elevation=.2,zoom=1,selected=null,points=[],timer=null,drag=null;
let width=0,height=0;
const radius=Math.max(...DATA.corners.map(p=>Math.hypot(...p)),1);
const cfg=DATA.summary.config;
const fmt=(x,unit='',digits=4)=>{if(!Number.isFinite(x))return '—';let [m,e]=Number(x).toPrecision(digits).split('e');if(m.includes('.'))m=m.replace(/0+$/,'').replace(/\.$/,'');return m+(e?'e'+e:'')+(unit?' '+unit:'')};
function engineering(x,unit){if(x===0)return '0 '+unit;const e=Math.max(-12,Math.min(3,3*Math.floor(Math.log10(Math.abs(x))/3)));const p={'-12':'p','-9':'n','-6':'µ','-3':'m','0':'','3':'k'};return fmt(x/Math.pow(10,e),(p[e]||'')+unit,4)}
function color(v){const stops=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];let f=DATA.potentialMax===DATA.potentialMin?.5:(v-DATA.potentialMin)/(DATA.potentialMax-DATA.potentialMin);f=Math.max(0,Math.min(1,f))*4;const i=Math.min(3,Math.floor(f)),u=f-i;return 'rgb('+stops[i].map((c,j)=>Math.round(c*(1-u)+stops[i+1][j]*u)).join(',')+')'}
function project(p){const x=Math.cos(yaw)*p[0]-Math.sin(yaw)*p[1],d=Math.sin(yaw)*p[0]+Math.cos(yaw)*p[1];const v=Math.cos(elevation)*p[2]-Math.sin(elevation)*d;const depth=Math.sin(elevation)*p[2]+Math.cos(elevation)*d;const scale=Math.min(width*.83,height*.83)/(2*radius)*zoom;return [width/2+x*scale,height/2-v*scale,depth]}
function line(a,b,c,w=1,alpha=1){const p=project(a),q=project(b);ctx.globalAlpha=alpha;ctx.strokeStyle=c;ctx.lineWidth=w;ctx.beginPath();ctx.moveTo(p[0],p[1]);ctx.lineTo(q[0],q[1]);ctx.stroke();ctx.globalAlpha=1}
function polygon(ids,fill){ctx.fillStyle=fill;ctx.beginPath();ids.forEach((id,k)=>{const p=project(DATA.corners[id]);k?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1])});ctx.closePath();ctx.fill()}
function draw(){ctx.clearRect(0,0,width,height);polygon([0,2,6,4],'rgba(78,113,136,.10)');polygon([1,3,7,5],'rgba(78,113,136,.16)');for(let a=0;a<8;a++)for(let b=a+1;b<8;b++)if([1,2,4].includes(a^b))line(DATA.corners[a],DATA.corners[b],'#b0c1cb',1);
if($('showCurrent').checked){const currents=DATA.edgeCurrent[index],maximum=Math.max(...currents.map(Math.abs),0);if(maximum>0){const order=currents.map((x,i)=>i).sort((a,b)=>Math.abs(currents[a])-Math.abs(currents[b]));for(const edge of order){const f=Math.abs(currents[edge])/maximum;if(f===0)continue;for(const seg of DATA.segments[edge])line(seg[0],seg[1],'#408eaf',.35+2*Math.sqrt(f),.10+.60*Math.sqrt(f));}}}
if($('showPath').checked)for(const edge of DATA.pathEdges[index])for(const seg of DATA.segments[edge])line(seg[0],seg[1],'#d74032',3.1,.95);
const items=[];if($('showHf').checked)DATA.hf.forEach((p,i)=>items.push({p:project(p),hf:true}));points=DATA.positions.map(project);points.forEach((p,i)=>{if(DATA.vacancy[index][i]||$('showO').checked)items.push({p,i,hf:false})});items.sort((a,b)=>a.p[2]-b.p[2]);for(const item of items){const [x,y]=item.p;if(item.hf){ctx.fillStyle='rgba(140,155,173,.30)';ctx.beginPath();ctx.arc(x,y,1.8,0,Math.PI*2);ctx.fill();continue}const i=item.i,v=DATA.vacancy[index][i];ctx.globalAlpha=v?1:.37;ctx.fillStyle=color(DATA.potential[index][i]);ctx.beginPath();ctx.arc(x,y,v?5.2:2.25,0,Math.PI*2);ctx.fill();ctx.globalAlpha=1;if(v){ctx.strokeStyle='#ed8d36';ctx.lineWidth=1.7;ctx.stroke()}if(i===selected){ctx.strokeStyle='#123b53';ctx.lineWidth=2;ctx.beginPath();ctx.arc(x,y,9,0,Math.PI*2);ctx.stroke()}}
ctx.fillStyle='#38576c';ctx.font='12px system-ui';for(const [label,ids]of[['Bottom electrode',[0,2,6,4]],['Top electrode',[1,3,7,5]]]){const projected=ids.map(i=>project(DATA.corners[i]));const anchor=projected.reduce((a,b)=>a[0]>b[0]?a:b);ctx.fillText(label,Math.min(width-105,anchor[0]+7),anchor[1]+4)}drawSpark()}
function drawSpark(){const w=spark.clientWidth,h=spark.clientHeight;sctx.clearRect(0,0,w,h);const rows=DATA.records,positive=rows.map(r=>r.equivalent_conductance_S).filter(g=>g>0);sctx.fillStyle='#617b8c';sctx.font='11px system-ui';sctx.fillText('Saved G(t), logarithmic scale',4,13);if(!positive.length){sctx.fillText('G = 0 in all saved states',4,40);return}const lo=Math.log10(Math.min(...positive)),hi=Math.log10(Math.max(...positive)),t0=rows[0].time_s,t1=rows[rows.length-1].time_s;const px=r=>40+(w-50)*(t1===t0?0:(r.time_s-t0)/(t1-t0)),py=r=>h-10-(h-31)*(hi===lo?.5:(Math.log10(r.equivalent_conductance_S)-lo)/(hi-lo));sctx.strokeStyle='#237f99';sctx.lineWidth=1.3;sctx.beginPath();let started=false;for(const r of rows){if(r.equivalent_conductance_S<=0){started=false;continue}if(started)sctx.lineTo(px(r),py(r));else{sctx.moveTo(px(r),py(r));started=true}}sctx.stroke();const r=rows[index];sctx.strokeStyle='#d57632';sctx.beginPath();sctx.moveTo(px(r),19);sctx.lineTo(px(r),h-7);sctx.stroke();if(r.equivalent_conductance_S>0){sctx.fillStyle='#d57632';sctx.beginPath();sctx.arc(px(r),py(r),3.5,0,Math.PI*2);sctx.fill()}}
function inspect(){if(selected===null)return;const r=DATA.records[index],i=selected,v=DATA.vacancy[index][i],region=DATA.bottom[i]?'bottom contact':DATA.top[i]?'top contact':'interior';$('inspect').innerHTML='<strong>Site '+i+' · '+region+'</strong><br>'+(v?'Vacancy label '+DATA.vacancyIds[index][i]:'Oxygen atom label '+DATA.oxygenIds[index][i])+'<br>Local potential: '+fmt(DATA.potential[index][i],'V',6)+'<br>Saved state '+(index+1)+'; event '+r.events}
function update(){const r=DATA.records[index];frame.value=index;$('frameLabel').textContent=(index+1)+' / '+DATA.records.length;const entries=[['Physical time',engineering(r.time_s,'s')],['KMC event count',String(r.events)],['Vacancies / oxide O',r.vacancies+' / '+r.oxide_oxygen],['Reservoir O: bottom / top',r.bottom_reservoir_oxygen+' / '+r.top_reservoir_oxygen],['Total oxygen atoms',String(r.total_oxygen)],['O / reaction charge residual',r.mass_residual+' / '+r.reaction_charge_residual_e+' e'],['Conductance G',engineering(r.equivalent_conductance_S,'S')],['Signed write current',engineering(r.terminal_current_A,'A')],['Actual / requested voltage',fmt(r.applied_voltage_V,'V')+' / '+fmt(r.requested_voltage_V,'V')],['Frozen-state read current',engineering(r.read_current_A,'A')]];$('stats').innerHTML=entries.map(([k,v])=>'<div class="kv"><span>'+k+'</span><strong>'+v+'</strong></div>').join('');$('criterion').innerHTML='<strong>Vacancy graph path: '+(r.vacancy_graph_spanning?'present':'absent')+'</strong><br>Declared model criterion: '+(r.model_formation_criterion?'met':'not met');$('criterion').className='status '+(r.model_formation_criterion?'good':'');$('potentialLegend').textContent=fmt(DATA.potentialMin,'V',3)+' to '+fmt(DATA.potentialMax,'V',3)+' (site potential)';inspect();draw()}
function stop(){if(timer!==null)clearInterval(timer);timer=null;$('play').textContent='Play'}
function play(){stop();if(index===DATA.records.length-1)index=0;$('play').textContent='Pause';timer=setInterval(()=>{if(index>=DATA.records.length-1){stop();return}index++;update()},1000/Number($('speed').value));update()}
function resize(){const ratio=window.devicePixelRatio||1;width=canvas.clientWidth;height=canvas.clientHeight;canvas.width=Math.round(width*ratio);canvas.height=Math.round(height*ratio);ctx.setTransform(ratio,0,0,ratio,0,0);spark.width=Math.round(spark.clientWidth*ratio);spark.height=Math.round(spark.clientHeight*ratio);sctx.setTransform(ratio,0,0,ratio,0,0);draw()}
frame.addEventListener('input',()=>{stop();index=Number(frame.value);update()});$('prev').onclick=()=>{stop();index=Math.max(0,index-1);update()};$('next').onclick=()=>{stop();index=Math.min(DATA.records.length-1,index+1);update()};$('play').onclick=()=>timer===null?play():stop();$('speed').onchange=()=>{if(timer!==null)play()};const firstCriterion=DATA.records.findIndex(r=>r.model_formation_criterion);$('formed').disabled=firstCriterion<0;$('formed').onclick=()=>{if(firstCriterion>=0){stop();index=firstCriterion;update()}};$('final').onclick=()=>{stop();index=DATA.records.length-1;update()};$('reset').onclick=()=>{yaw=.65;elevation=.2;zoom=1;selected=null;$('inspect').textContent='Click a site to inspect its saved identity and local potential.';draw()};for(const id of ['showO','showHf','showPath','showCurrent'])$(id).onchange=draw;
canvas.addEventListener('pointerdown',e=>{canvas.setPointerCapture(e.pointerId);drag={x:e.clientX,y:e.clientY,startX:e.clientX,startY:e.clientY}});canvas.addEventListener('pointermove',e=>{if(!drag)return;yaw+=(e.clientX-drag.x)*.009;elevation=Math.max(-1.3,Math.min(1.3,elevation+(e.clientY-drag.y)*.007));drag.x=e.clientX;drag.y=e.clientY;draw()});canvas.addEventListener('pointerup',e=>{if(drag&&Math.hypot(e.clientX-drag.startX,e.clientY-drag.startY)<4){const rect=canvas.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top;let best=12,bestDepth=-Infinity;selected=null;points.forEach((p,i)=>{const d=Math.hypot(p[0]-x,p[1]-y);if(d<best||(Math.abs(d-best)<1&&p[2]>bestDepth)){best=d;bestDepth=p[2];selected=i}});inspect();draw()}drag=null});canvas.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.4,Math.min(3.5,zoom*Math.exp(-e.deltaY*.001)));draw()},{passive:false});window.addEventListener('keydown',e=>{if(e.key==='ArrowRight'){$('next').click();e.preventDefault()}else if(e.key==='ArrowLeft'){$('prev').click();e.preventDefault()}});
const termination=DATA.summary.termination||{};$('footer').textContent='ILLUSTRATIVE parameters · Fixed q = '+cfg.charge+'e · '+cfg.temperature_K+' K, isothermal · No space-charge Poisson, quantum transport, or thermal feedback. Initial vacancy graph path: '+(DATA.summary.initial_spanning?'present':'absent')+'. Run status: '+String(termination.status||'saved output').replaceAll('_',' ')+'. The model criterion combines the recorded path and declared conductance/gain thresholds; it is not experimental forming validation. Read current is a frozen-configuration linear diagnostic.';
new ResizeObserver(resize).observe(canvas.parentElement);update();resize();
</script></body></html>'''
