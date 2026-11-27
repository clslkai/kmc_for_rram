"""Run a reproducible pulse experiment and export actual states, not stylized filaments."""
import argparse
import csv
import json
import os
from pathlib import Path
import numpy as np
from .model import Config, Simulator


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, help="JSON overrides of Config fields")
    p.add_argument("--shape", nargs=3, type=int, default=None, metavar=("NX", "NY", "NZ"))
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--initial", choices=("pristine", "random", "filament"), default="random")
    p.add_argument("--fraction", type=float, default=0.04)
    p.add_argument("--pulse", nargs=2, type=float, action="append", metavar=("VOLTAGE", "DURATION_S"),
                   help="repeat for piecewise constant bias; e.g. --pulse 2 0.001 --pulse -2 0.001")
    p.add_argument("--max-events", type=int, default=2000)
    p.add_argument("--sample-every", type=int, default=20)
    p.add_argument("--output", type=Path, default=Path("results/demo"))
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args(argv)
    if args.max_events < 1 or args.sample_every < 1: p.error("event limits must be positive")
    data = json.loads(args.config.read_text()) if args.config else {}
    if "shape" in data: data["shape"] = tuple(data["shape"])
    if args.shape: data["shape"] = tuple(args.shape)
    if args.seed is not None: data["seed"] = args.seed
    config = Config(**data)
    pulses = args.pulse or [(2.0, 0.000001), (-2.0, 0.000001)]
    if any(not np.isfinite(v) or not np.isfinite(t) or t <= 0 for v,t in pulses):
        p.error("each pulse needs finite voltage and positive duration")
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    sim = Simulator(config); sim.initialize(args.initial, args.fraction)
    records = []; vacancies = []; interstitials = []; potentials = []; temperatures = []
    def capture(v):
        records.append(sim.snapshot(v))
        vacancies.append(sim.vacancy.copy().reshape(config.shape))
        interstitials.append(sim.interstitial.copy().reshape(config.shape))
        potentials.append(sim.phi.copy().reshape(config.shape))
        temperatures.append(sim.temperature.copy().reshape(config.shape))
    end = 0.0
    status = "completed_waveform"
    event_log = []
    for v, duration in pulses:
        capture(v)
        end += duration
        while sim.time_s < end and sim.event_count < args.max_events:
            event = sim.step(v, end)
            if event is not None:
                event_log.append((sim.time_s, event.kind, event.source, event.target, event.delta_eV, event.rate_Hz))
                if sim.event_count % args.sample_every == 0: capture(v)
        capture(v)
        if sim.event_count >= args.max_events and sim.time_s < end:
            status = "event_limit_before_waveform_end"
            break
    with (out / "trace.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=records[0].keys()); writer.writeheader(); writer.writerows(records)
    with (out / "events.csv").open("w", newline="") as f:
        writer = csv.writer(f); writer.writerow(["time_s", "kind", "source", "target", "delta_eV", "rate_Hz"]); writer.writerows(event_log)
    np.savez_compressed(out / "states.npz", time_s=np.array([r["time_s"] for r in records]),
                        vacancy=vacancies, interstitial=interstitials, potential_V=potentials,
                        temperature_K=temperatures)
    summary = {"config": sim.config_dict(), "initial": args.initial, "fraction": args.fraction,
               "pulses_V_s": pulses, "status": status, "event_counts": sim.counts,
               "final": records[-1], "warning": "Illustrative closed-oxygen coarse-grained model; not calibrated HfO2 predictions."}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    if not args.no_plot: plot_results(out, config, records, vacancies, interstitials)
    print(json.dumps({"output": str(out.resolve()), "status": status, "events": sim.event_count,
                      "time_s": sim.time_s, "counts": sim.counts}, indent=2))


def plot_results(out, config, records, vacancies, interstitials):
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/vcm3d-matplotlib")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp/vcm3d-cache")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(12, 8))
    for n, idx, title in [(1, 0, "Initial state"), (2, -1, "Final state")]:
        ax = fig.add_subplot(2, 2, n, projection="3d")
        for state, color, size, name in [(vacancies[idx], "#de7b20", 24, "Vacancy"),
                                       (interstitials[idx], "#207caf", 10, "O interstitial")]:
            xyz = np.argwhere(state) * config.spacing_m * 1e9
            ax.scatter(*xyz.T, c=color, s=size, alpha=.7, label=name)
        ax.set(xlim=(0,(config.shape[0]-1)*config.spacing_m*1e9),
               ylim=(0,(config.shape[1]-1)*config.spacing_m*1e9),
               zlim=(-.5*config.spacing_m*1e9,(config.shape[2]-.5)*config.spacing_m*1e9),
               xlabel="x (nm)", ylabel="y (nm)", zlabel="z (nm)", title=title)
        ax.legend(fontsize=8)
    times = np.array([r["time_s"] for r in records])
    ax = fig.add_subplot(2,2,3)
    ax.plot(times, [r["vacancies"] for r in records], label="Vacancies")
    ax.plot(times, [r["interstitials"] for r in records], "--", label="Interstitials")
    ax.set(xlabel="Time (s)", ylabel="Defect count"); ax.legend()
    ax = fig.add_subplot(2,2,4)
    ax.semilogy(times, [r["conductance_S"] for r in records], label="Conductance")
    ax.set(xlabel="Time (s)", ylabel="Frozen-state conductance (S)")
    twin = ax.twinx(); twin.plot(times, [r["max_temperature_K"] for r in records], color="red")
    twin.set_ylabel("Maximum temperature (K)", color="red")
    fig.suptitle("3D VCM reference KMC — illustrative, not calibrated")
    fig.tight_layout(); fig.savefig(out/"overview.png", dpi=160); plt.close(fig)

if __name__ == "__main__":
    main()
