"""Predeclared synthetic interface/transport demonstrations, including controls.

No seed is selected by whether it forms a path. Every scheduled outcome,
including model-domain failures and event-cap censoring, is retained.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
import json
from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

from dftkmc.catalog import Catalog
from dftkmc.coupled import CoupledConfig
from dftkmc.coupled_workflow import run_coupled
from dftkmc.exchange import ExchangeConfig, ReservoirConfig
from dftkmc.transport import TransportConfig
from dftkmc.workflow import _write_json, catalog_digest


def first_persistent_onset(transitions, end_time, hold_time):
    start = None
    for row in transitions:
        if row["criterion_met"]:
            start = row["time_s"]
        elif start is not None:
            if row["time_s"]-start >= hold_time:
                return start
            start = None
    return start if start is not None and end_time-start >= hold_time else None


def run_case(task):
    root = Path(task["output"])
    config = CoupledConfig.from_dict(task["config"])
    catalog = Catalog.load(root/"catalog.virtual.json", allow_illustrative=True)
    try:
        summary = run_coupled(catalog, config, pulses=task["pulses"], n_vacancies=0,
            max_events=task["max_events"], sample_every=task["sample_every"],
            output=root/task["name"], plot=task["plot"], allow_illustrative=True)
        final, termination, forming = summary["final"], summary["termination"], summary["formation"]
        return dict(name=task["name"], role=task["role"], seed=config.seed,
            status=termination["status"], error=termination["error"],
            target_reached=termination["target_reached"], time_s=final["time_s"], events=final["events"],
            initial_vacancies=summary["initial"]["vacancies"], initial_spanning=summary["initial_spanning"],
            final_vacancies=final["vacancies"], final_spanning=final["vacancy_graph_spanning"],
            initial_conductance_S=summary["initial"]["equivalent_conductance_S"],
            peak_conductance_S=final["peak_conductance_S"], final_conductance_S=final["equivalent_conductance_S"],
            first_spanning_time_s=forming["first_spanning_time_s"],
            first_model_criterion_time_s=forming["first_criterion_time_s"],
            first_persistent_criterion_onset_s=first_persistent_onset(forming["criterion_transitions"],
                final["time_s"], task["diagnostic_hold_time_s"]),
            longest_observed_criterion_duration_s=forming["longest_observed_criterion_duration_s"],
            final_model_criterion=forming["final_criterion_met"],
            oxygen_out=summary["event_counts"].get("oxygen_out", 0),
            oxygen_in=summary["event_counts"].get("oxygen_in", 0),
            hopping=summary["event_counts"].get("hop", 0),
            conservation=summary["conservation"], numerical_checks=summary["numerical_checks"])
    except Exception as exc:
        # A failed worker is an outcome, never an invitation to resample seeds.
        failed = dict(name=task["name"], role=task["role"], seed=config.seed,
                      status="workflow_error", target_reached=False,
                      error=f"{type(exc).__name__}: {exc}")
        _write_json(root/(task["name"]+".failure.json"), failed)
        return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--allow-illustrative", action="store_true")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    if not args.allow_illustrative:
        parser.error("explicit --allow-illustrative is required")
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("--output must be a new or empty directory")
    catalog = Catalog.load(PROJECT/"examples/virtual_data_v3_1/catalog.virtual.json", allow_illustrative=True)
    args.output.mkdir(parents=True, exist_ok=True)
    catalog.write(args.output/"catalog.virtual.json")
    reservoirs = ExchangeConfig(
        vacancy_formation_offset_eV=.4, reference_basis=0,
        bottom=ReservoirConfig(capacity=96, initial_oxygen=0, mu0_eV=-.15, enabled=False),
        top=ReservoirConfig(capacity=96, initial_oxygen=0, mu0_eV=-.15, enabled=True),
        provenance="Synthetic device-protocol chemistry; chosen before the scheduled runs, not fitted to DFT or a forming curve.")
    base = CoupledConfig(repeats=(2, 2, 4), charge=2, temperature_K=800., seed=17,
        bottom_potential_V=0., top_potential_V=.5, exchange=reservoirs,
        transport=TransportConfig(current_compliance_A=5e-8),
        formation_min_conductance_S=1e-9, formation_min_gain=1000.)
    pulses = [[.5, 6e-6], [-.5, 2e-6]]
    configurations = [(f"biased_seed{seed}", "primary", replace(base, seed=seed), pulses) for seed in (17, 29, 43)]
    configurations += [
        ("zero_bias_seed17", "zero_bias_control", base, [[0., 8e-6]]),
        ("exchange_blocked_seed17", "exchange_blocked_control",
         replace(base, exchange=replace(reservoirs, top=replace(reservoirs.top, enabled=False))), pulses),
        ("constant_conductance_seed17", "constant_conductance_control",
         replace(base, transport=replace(base.transport,
             vacancy_conductance_S=base.transport.background_conductance_S)), pulses),
    ]
    tasks = []
    for name, role, config, waveform in configurations:
        config.write(args.output/(name+".config.json"))
        tasks.append(dict(name=name, role=role, config=config.to_dict(), pulses=waveform,
                          output=str(args.output.resolve()), max_events=30000, sample_every=20,
                          diagnostic_hold_time_s=1e-8, plot=not args.no_plot))
    protocol = dict(schema_version="coupled-demo-protocol-1", calibration_status="illustrative",
        catalog_sha256=catalog_digest(catalog), representative_case="biased_seed17",
        initial_state="No vacancies; initially nonspanning oxide. No seed selection or prescribed filament.",
        requested_common_end_time_s=8e-6, diagnostic_hold_time_s=1e-8, tasks=tasks,
        interpretation=[
            "Protocol is saved before any trajectory. Primary seeds are predeclared and every outcome is retained.",
            "Three primary replicas and one replica per control are software demonstrations, not a converged forming distribution.",
            "The joint geometry/conductance criterion and 10 ns duration diagnostic are predeclared engineering observables, not validated physical memory retention.",
            "The constant-conductance control changes both readout mapping and ionic voltage feedback; it is an engineering ablation, not an isolated measurement of one material mechanism.",
            "A negative pulse is a waveform test; RESET is not assumed to occur.",
        ])
    _write_json(args.output/"protocol.json", protocol)
    outcomes = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_case, task): task["name"] for task in tasks}
        for future in as_completed(futures):
            outcome = future.result()
            outcomes[outcome["name"]] = outcome
            _write_json(args.output/"progress.json", [outcomes[t["name"]] for t in tasks if t["name"] in outcomes])
            print(json.dumps(outcome, allow_nan=False), flush=True)
    ordered = [outcomes[t["name"]] for t in tasks]
    report = dict(code_version="0.4.0", calibration_status="illustrative",
        protocol_file="protocol.json", representative_case="biased_seed17", outcomes=ordered,
        all_complete=all(r["target_reached"] for r in ordered),
        planned_cases=len(tasks), completed_cases=sum(r["target_reached"] for r in ordered),
        interpretation=protocol["interpretation"])
    _write_json(args.output/"demo_report.json", report)
    (args.output/"progress.json").unlink(missing_ok=True)
    return 0 if report["all_complete"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
