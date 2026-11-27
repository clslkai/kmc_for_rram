"""Plot saved replica summaries without rerunning or modifying KMC results.

Usage from the project root:
    python examples/dftkmc/plot_study.py examples/dftkmc/mechanism_study/study.json
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.study.read_text(encoding="utf-8"))
    output = args.output or args.study.with_name("mechanism_comparison.png")
    groups = data["groups"]
    charges = list(dict.fromkeys(g["charge"] for g in groups))
    ablations = list(dict.fromkeys(g["ablation"] for g in groups))
    colors = plt.get_cmap("tab10")
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.9))
    for index, charge in enumerate(charges):
        offset = (index-(len(charges)-1)/2)*0.22
        for group in (g for g in groups if g["charge"] == charge):
            x = ablations.index(group["ablation"])+offset
            stats = group["statistics"]
            if not stats["all_replicas_complete"]:
                continue  # Censored groups must never be plotted as complete inference.
            diffusion = stats["diffusion"]
            values = [
                (diffusion.get("trace_over_3_m2_s"), diffusion.get("trace_over_3_ci_m2_s")),
                (stats["scalar_estimates"]["pair_count"]["mean"], stats["scalar_estimates"]["pair_count"]["ci"]),
                (stats["scalar_estimates"]["largest_cluster"]["mean"], stats["scalar_estimates"]["largest_cluster"]["ci"]),
            ]
            for panel, (axis, (mean, interval)) in enumerate(zip(axes, values)):
                if mean is None:
                    continue
                scale = 1e12 if panel == 0 else 1.
                if interval is not None:
                    low, high = [value*scale for value in interval]
                    axis.vlines(x, low, high, color=colors(index), linewidth=1.6)
                    axis.hlines([low, high], x-.04, x+.04, color=colors(index), linewidth=1.2)
                axis.plot(x, mean*scale, "o", color=colors(index), markersize=6,
                          label=f"q = {charge}" if group["ablation"] == ablations[0] else None)
    labels = {"noninteracting": "No pair / no TS", "pair-only": "Pair only", "pair-environment": "Pair + TS"}
    for axis, title, ylabel in zip(axes,
            ("Finite-time tracer motion", "Local aggregation", "Cluster size"),
            (r"tr[$D(t)$] / 3 ($10^{-12}$ m$^2$/s)", "Mean vacancy pair count", "Mean largest cluster (vacancies)")):
        axis.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=14)
        axis.set_xticks(range(len(ablations)), [labels.get(x,x) for x in ablations], rotation=18)
        axis.set_ylabel(ylabel)
        axis.set_ylim(bottom=0)
        axis.grid(axis="y", alpha=.18)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=9)
    status = data["calibration_status"]
    headline = "SYNTHETIC PARAMETERS — NOT DFT-CALIBRATED" if status == "illustrative" else f"Declared calibration: {status}"
    fig.suptitle(headline, x=.07, ha="left", fontsize=14, fontweight="bold")
    settings = (f"{data['config']['temperature_K']:g} K  |  {data['duration_s']*1e6:g} microseconds  |  "
                f"{data['n_vacancies']} vacancies  |  {data['replicas_per_group']} independent replicas per group")
    fig.text(.07, .89, settings, fontsize=10, color="#454545")
    confidence = data["study_controls"]["confidence"]*100
    fig.text(.07, .035, f"{confidence:g}% intervals: replica bootstrap for tracer motion; Student t for cluster observables.\n"
             "TS = transition-state environment correction. One observation time does not establish long-time diffusion.",
             fontsize=9, color="#454545")
    fig.subplots_adjust(left=.07, right=.985, bottom=.25, top=.77, wspace=.38)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
