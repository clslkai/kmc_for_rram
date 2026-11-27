"""Reproduce three declared synthetic device-boundary trajectories.

Run from the project root:
python examples/device_boundary_v3_2/run_demo.py --allow-illustrative --output results/device_boundary
These are paired demonstration cases, not independent statistical replicas.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

from dftkmc.catalog import Catalog
from dftkmc.device import DeviceConfig, ElectrodeConfig
from dftkmc.model import SimulationConfig
from dftkmc.workflow import run_simulation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-illustrative", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--duration", type=float, default=1e-6)
    parser.add_argument("--max-events", type=int, default=20000)
    args = parser.parse_args()
    if not args.allow_illustrative:
        parser.error("--allow-illustrative is required: these parameters are synthetic, not completed DFT data")
    if not np.isfinite(args.duration) or args.duration <= 0 or args.max_events < 1:
        parser.error("use a positive finite duration and a positive event cap")
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("--output must be a new or empty directory")
    catalog = Catalog.load(PROJECT / "examples/virtual_data_v3_1/catalog.virtual.json", allow_illustrative=True)
    args.output.mkdir(parents=True, exist_ok=True)
    catalog.write(args.output / "catalog.virtual.json")
    cases = []
    for name, top_voltage in (("negative", -1.0), ("zero", 0.0), ("positive", 1.0)):
        device = DeviceConfig(bottom=ElectrodeConfig(potential_V=0.0),
                              top=ElectrodeConfig(potential_V=top_voltage))
        device.write(args.output / f"device_{name}.json")
        config = SimulationConfig(repeats=(3, 3, 8), charge=2, temperature_K=800.,
                                  seed=17, boundary="slab", device=device)
        result = run_simulation(catalog, config, n_vacancies=8, duration_s=args.duration,
                                max_events=args.max_events, sample_every=10,
                                output=args.output/name, plot=True, allow_illustrative=True)
        summary = json.loads((args.output/name/"summary.json").read_text())
        with np.load(args.output/name/"states.npz") as state:
            initial_depth = state["device_depth_A"][state["vacancy"][0]]
            conserved = bool(np.all(state["vacancy"].sum(axis=1) == 8))
        case = dict(name=name, top_voltage_V=top_voltage, bottom_voltage_V=0.,
                    effective_field_V_nm=result["effective_field_V_nm"],
                    thickness_nm=summary["device_boundary"]["thickness_nm"],
                    initial_vacancy_indices=result["initial_vacancy_indices"],
                    initial_mean_depth_A=float(initial_depth.mean()),
                    final_mean_depth_A=result["final"]["mean_vacancy_depth_A"],
                    events=result["events"], time_s=result["time_s"],
                    target_reached=result["target_reached"], status=result["status"],
                    vacancy_count_conserved_in_saved_states=conserved,
                    bottom_contact_vacancies=result["final"]["bottom_contact_vacancies"],
                    interior_vacancies=result["final"]["interior_vacancies"],
                    top_contact_vacancies=result["final"]["top_contact_vacancies"])
        cases.append(case)
        print(json.dumps(case, ensure_ascii=False, allow_nan=False), flush=True)
    report = dict(
        version="0.3.2", calibration_status="illustrative", seed=17,
        repeats=[3, 3, 8], temperature_K=800., fixed_charge=2, n_vacancies=8,
        requested_duration_s=args.duration, max_events=args.max_events, cases=cases,
        all_complete=all(case["target_reached"] for case in cases),
        same_initial_sites=all(case["initial_vacancy_indices"] == cases[0]["initial_vacancy_indices"] for case in cases),
        interpretation=("Three paired synthetic boundary demonstrations at a common physical time. "
                        "Same seed and initial-site rule across voltages; not independent statistical replicas. "
                        "Observed vacancy redistribution is not electronic current, filament formation or a material prediction."),
    )
    (args.output/"demo_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
    return 0 if report["all_complete"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
