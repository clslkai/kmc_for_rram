"""End-to-end virtual-data experiments with independent finite-time references.

This suite only accepts explicitly illustrative data. It does not run DFT.
Run: python -m dftkmc.virtual_suite --allow-illustrative --output results/virtual
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import math
import multiprocessing
from numbers import Integral
from pathlib import Path
import sys
import time

import numpy as np
from scipy.stats import t as student_t

from . import __version__
from .analysis import compare_groups, reached_target, summarize_replicas
from .catalog import Catalog
from .model import KB, SimulationConfig, Simulator
from .workflow import (ABLATIONS, _graph_coverage, _json_default, _new_output,
                       _write_csv, _write_json, catalog_digest, run_simulation)


PROFILES = {
    "standard": dict(temperatures=[700., 800., 900.], benchmark_hops=[50, 200],
                     benchmark_replicas=24, mechanism_replicas=16,
                     mechanism_charges=[0, 2], mechanism_hops=40,
                     field_replicas=16, field_hops=200, bootstrap_samples=2000),
    "smoke": dict(temperatures=[800.], benchmark_hops=[2, 5],
                  benchmark_replicas=2, mechanism_replicas=2,
                  mechanism_charges=[2], mechanism_hops=2,
                  field_replicas=2, field_hops=5, bootstrap_samples=100),
}
_WORKER_CATALOG = None


def _positive_integer(value, name):
    if not isinstance(value, Integral) or isinstance(value, (bool, np.bool_)) or value < 1:
        raise ValueError(name + " must be a positive integer")


def build_suite_plan(catalog, profile="standard", base_seed=20261111, max_events=5000):
    """Freeze physical durations and independent random streams before sampling."""
    from .reference import periodic_tracer_reference, finite_time_tracer_reference
    if profile not in PROFILES:
        raise ValueError("profile must be standard or smoke")
    if not isinstance(base_seed, Integral) or isinstance(base_seed, (bool, np.bool_)) or base_seed < 0:
        raise ValueError("base_seed must be a nonnegative integer")
    _positive_integer(max_events, "max_events")
    catalog.validate(allow_illustrative=True)
    if catalog.data["metadata"]["calibration_status"] != "illustrative":
        raise ValueError("virtual suite accepts illustrative data only")
    if not {0, 2}.issubset(catalog.charges):
        raise ValueError("virtual suite requires both q=0 and q=2 branches")
    profile_data = dict(PROFILES[profile])
    conditions, references = [], {}
    for q in (0, 2):
        for temperature in profile_data["temperatures"]:
            reference = periodic_tracer_reference(catalog, q, temperature)
            references[(q, temperature)] = reference
            rate = reference["mean_jump_rate_Hz"]
            if not math.isfinite(rate) or rate <= 0:
                raise ValueError("the reference must have a positive finite mean jump rate")
            for hops in profile_data["benchmark_hops"]:
                duration = hops/rate
                finite = finite_time_tracer_reference(catalog, q, temperature, duration)
                conditions.append(dict(
                    id=f"diffusion_q{q}_T{temperature:g}_h{hops}",
                    family="diffusion", charge=q, temperature_K=temperature,
                    field_V_nm=[0., 0., 0.], pair_scale=0., environment_scale=0.,
                    n_vacancies=1, initial_mode="equilibrium_single",
                    basis_probabilities=reference["basis_probabilities"],
                    duration_s=duration, target_expected_hops=hops,
                    n_replicas=profile_data["benchmark_replicas"],
                    reference_long_time=reference, reference_finite_time=finite,
                ))
    fast_rate = max(references[(q, 800.)]["mean_jump_rate_Hz"] for q in (0, 2))
    for q in profile_data["mechanism_charges"]:
        for ablation, (pair_scale, environment_scale) in ABLATIONS.items():
            conditions.append(dict(
                id=f"mechanism_q{q}_{ablation}", family="mechanism", ablation=ablation,
                charge=q, temperature_K=800., field_V_nm=[0., 0., 0.],
                pair_scale=pair_scale, environment_scale=environment_scale,
                n_vacancies=8, initial_mode="uniform_random",
                duration_s=profile_data["mechanism_hops"]/fast_rate,
                n_replicas=profile_data["mechanism_replicas"],
            ))
    for q in (0, 2):
        for field in (-0.15, 0., 0.15):
            conditions.append(dict(
                id=f"field_q{q}_F{field:+.2f}", family="field", charge=q,
                temperature_K=800., field_V_nm=[0., 0., field],
                pair_scale=0., environment_scale=0., n_vacancies=1,
                initial_mode="equilibrium_single",
                basis_probabilities=references[(q, 800.)]["basis_probabilities"],
                duration_s=profile_data["field_hops"]/fast_rate,
                n_replicas=profile_data["field_replicas"],
                seed_pairing="common streams across neutral fields" if q == 0 else "independent groups",
            ))
    sequence = np.random.SeedSequence(int(base_seed))
    neutral_streams = None
    for condition in conditions:
        if condition["family"] == "field" and condition["charge"] == 0 and neutral_streams is not None:
            streams = [dict(item) for item in neutral_streams]
        else:
            streams = []
            for replica_index, child in enumerate(sequence.spawn(condition["n_replicas"])):
                initialization, trajectory = child.spawn(2)
                streams.append(dict(
                    replica_index=replica_index,
                    initialization_seed=int(initialization.generate_state(1, dtype=np.uint64)[0]),
                    trajectory_seed=int(trajectory.generate_state(1, dtype=np.uint64)[0]),
                ))
            if condition["family"] == "field" and condition["charge"] == 0:
                neutral_streams = streams
        condition["replica_seeds"] = streams
        condition["bootstrap_seed"] = int(sequence.spawn(1)[0].generate_state(1, dtype=np.uint64)[0])
        condition["max_events"] = int(max_events)
        condition["repeats"] = [3, 3, 3]
    seed_pairs = {(s["initialization_seed"], s["trajectory_seed"])
                  for c in conditions for s in c["replica_seeds"]}
    if any(a == b for a, b in seed_pairs):
        raise RuntimeError("Unexpected initialization/trajectory seed collision")
    return dict(
        schema_version="dftkmc-virtual-plan-1", profile=profile, base_seed=int(base_seed),
        catalog_sha256=catalog_digest(catalog), calibration_status="illustrative",
        conditions=conditions, profile_controls=profile_data, confidence=.95,
        n_conditions=len(conditions),
        n_trajectory_executions=sum(c["n_replicas"] for c in conditions),
        n_distinct_stream_pairs=len(seed_pairs),
        duration_policy="Every physical duration is fixed before sampling. Expected-hop targets select the duration, never the observed event count.",
        initialization_policy="Separate SeedSequence children for initial positions, trajectory randomness and bootstrap. Neutral-field controls intentionally reuse paired initial/trajectory streams.",
        **_graph_coverage(catalog),
    )


def _simulate_replica(catalog, condition, stream):
    config = SimulationConfig(
        repeats=tuple(condition["repeats"]), charge=condition["charge"],
        temperature_K=condition["temperature_K"], seed=stream["trajectory_seed"],
        field_V_nm=tuple(condition["field_V_nm"]), pair_scale=condition["pair_scale"],
        environment_scale=condition["environment_scale"],
    )
    base = dict(condition_id=condition["id"], replica_index=stream["replica_index"],
                initialization_seed=stream["initialization_seed"], seed=config.seed,
                duration_s=condition["duration_s"], charge=config.charge,
                family=condition["family"])
    try:
        simulator = Simulator(catalog, config)
        initial_rng = np.random.default_rng(stream["initialization_seed"])
        probability = None
        if condition["initial_mode"] == "equilibrium_single":
            probability = np.asarray(condition["basis_probabilities"])[simulator.lattice.basis_id]
            probability = probability/probability.sum()
        indices = np.sort(initial_rng.choice(simulator.n, size=condition["n_vacancies"],
                                            replace=False, p=probability))
        simulator.initialize(vacancy_indices=indices)
        trajectory = simulator.run(condition["duration_s"], max_events=condition["max_events"],
                                   sample_every=condition["max_events"]+1)
        final = simulator.snapshot()
        final["z_displacement_A"] = final["mean_displacement_A"][2]
        return dict(
            **base, status=trajectory["status"], time_s=simulator.time_s,
            target_reached=reached_target(simulator.time_s, condition["duration_s"]),
            events=simulator.event_count, initial_vacancy_indices=indices.tolist(),
            final_vacancy_indices=np.flatnonzero(simulator.vacancy).tolist(),
            final=final, displacements_A=simulator.displacements_A.tolist(),
            channel_counts=dict(simulator.channel_counts),
        )
    except (ValueError, RuntimeError, ArithmeticError) as exc:
        return dict(**base, status="failed", time_s=None, target_reached=False,
                    events=None, final={}, error=str(exc))


def _worker_init(catalog_data):
    global _WORKER_CATALOG
    _WORKER_CATALOG = Catalog(catalog_data)
    _WORKER_CATALOG.validate(allow_illustrative=True)


def _worker(job):
    return _simulate_replica(_WORKER_CATALOG, *job)


def _mean_interval(values, confidence=.95):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return dict(mean=None, ci=None)
    mean = float(values.mean())
    if len(values) < 2:
        return dict(mean=mean, ci=None)
    half = float(student_t.ppf((1+confidence)/2, len(values)-1) *
                 values.std(ddof=1)/np.sqrt(len(values)))
    return dict(mean=mean, ci=[mean-half, mean+half])


def summarize_suite(plan, records):
    """Infer only from complete groups; stochastic CI coverage is diagnostic."""
    groups, comparisons, group_records = [], [], {}
    for condition in plan["conditions"]:
        runs = sorted([r for r in records if r["condition_id"] == condition["id"]],
                      key=lambda r: r["replica_index"])
        if len(runs) != condition["n_replicas"]:
            raise ValueError("Missing replica records for " + condition["id"])
        group_records[condition["id"]] = runs
        stats = summarize_replicas(
            runs, duration_s=condition["duration_s"], boundary="bulk",
            field_V_nm=condition["field_V_nm"], confidence=plan["confidence"],
            bootstrap_samples=plan["profile_controls"]["bootstrap_samples"],
            bootstrap_seed=condition["bootstrap_seed"],
        )
        group = dict(condition={key: value for key, value in condition.items()
                                if key != "replica_seeds"}, statistics=stats)
        if stats["all_replicas_complete"]:
            group["z_displacement_A"] = _mean_interval(
                [r["final"]["z_displacement_A"] for r in runs], plan["confidence"])
        else:
            group["z_displacement_A"] = dict(mean=None, ci=None, reason="Incomplete group")
        if condition["family"] == "diffusion":
            finite = condition["reference_finite_time"]["trace_over_3_m2_s"]
            long_time = condition["reference_long_time"]["trace_over_3_m2_s"]
            diffusion = stats["diffusion"]
            estimate = diffusion.get("trace_over_3_m2_s")
            interval = diffusion.get("trace_over_3_ci_m2_s")
            comparison = dict(
                condition_id=condition["id"],
                reference_finite_time_m2_s=finite, reference_long_time_m2_s=long_time,
                estimate_m2_s=estimate, estimate_ci_m2_s=interval,
                relative_error_vs_finite_time=None if estimate is None else estimate/finite-1,
                finite_time_bias_vs_limit=finite/long_time-1,
                ci_contains_finite_time_reference=None if interval is None else interval[0] <= finite <= interval[1],
                interpretation="CI coverage and Monte Carlo error are diagnostics, not individually forced pass/fail gates. Finite-time reference matches the actual observation time.",
            )
            group["reference_comparison"] = comparison
            comparisons.append(comparison)
        groups.append(group)
    contrasts = []
    for q in (0, 2):
        for first, second in (("noninteracting", "pair-only"), ("pair-only", "pair-environment")):
            a, b = f"mechanism_q{q}_{first}", f"mechanism_q{q}_{second}"
            if a in group_records and b in group_records:
                contrasts.append(dict(first=a, second=b, statistics=compare_groups(
                    group_records[a], group_records[b],
                    duration_s=next(c["duration_s"] for c in plan["conditions"] if c["id"] == a),
                    metrics=("msd_A2", "pair_count", "largest_cluster", "events"),
                )))
    neutral = [c for c in plan["conditions"] if c["family"] == "field" and c["charge"] == 0]
    neutral_equal = None
    if neutral:
        baseline = group_records[neutral[0]["id"]]
        neutral_equal = all(
            a.get("target_reached") and b.get("target_reached")
            and a["seed"] == b["seed"]
            and a["initialization_seed"] == b["initialization_seed"]
            and a["initial_vacancy_indices"] == b["initial_vacancy_indices"]
            and a["events"] == b["events"]
            and a["displacements_A"] == b["displacements_A"]
            and a["final_vacancy_indices"] == b["final_vacancy_indices"]
            and a["channel_counts"] == b["channel_counts"]
            for c in neutral[1:] for a, b in zip(baseline, group_records[c["id"]])
        )
    arrhenius = []
    for q in (0, 2):
        selected = [g for g in groups if g["condition"]["family"] == "diffusion"
                    and g["condition"]["charge"] == q
                    and g["condition"]["target_expected_hops"] == max(plan["profile_controls"]["benchmark_hops"])]
        if len(selected) < 3:
            continue
        temperatures = np.asarray([g["condition"]["temperature_K"] for g in selected])
        for label, values in (
            ("exact_long_time_reference", [g["reference_comparison"]["reference_long_time_m2_s"] for g in selected]),
            ("finite_time_KMC_protocol", [g["reference_comparison"]["estimate_m2_s"] for g in selected]),
        ):
            if any(v is None or v <= 0 for v in values):
                continue
            coefficients = np.polyfit(1/temperatures, np.log(values), 1)
            arrhenius.append(dict(
                charge=q, source=label, apparent_activation_eV=float(-KB*coefficients[0]),
                temperatures_K=temperatures.tolist(), n_temperatures=len(temperatures),
                interpretation="Descriptive three-temperature slope of this synthetic network/protocol. It is not an individual NEB barrier, a fitted material value or a statistically converged activation energy.",
            ))
    all_complete = all(g["statistics"]["all_replicas_complete"] for g in groups)
    return dict(
        schema_version="dftkmc-virtual-results-1", calibration_status="illustrative",
        n_conditions=len(groups), n_trajectory_executions=len(records),
        n_distinct_stream_pairs=plan["n_distinct_stream_pairs"],
        total_events=sum(r.get("events") or 0 for r in records),
        all_trajectories_complete=all_complete,
        neutral_field_paired_invariance=neutral_equal,
        deterministic_checks_passed=all_complete and neutral_equal is True,
        groups=groups, reference_comparisons=comparisons,
        mechanism_contrasts=contrasts, apparent_arrhenius_slopes=arrhenius,
        statistical_diagnostics=dict(
            n_diffusion_conditions=len(comparisons),
            n_ci_available=sum(c["ci_contains_finite_time_reference"] is not None for c in comparisons),
            n_ci_cover_reference=sum(c["ci_contains_finite_time_reference"] is True for c in comparisons),
            caution="Nominal 95% intervals can miss the reference. These correlated/multiple diagnostics are not all required to pass; no seeds or data are adjusted to force coverage.",
        ),
        limitations=[
            "Synthetic input only; no DFT/NEB calculations or material claims.",
            "The reference is exact for the declared single-vacancy periodic rate network, not for real hafnia.",
            "The mechanism study starts from uniform random multi-vacancy configurations; final cluster metrics are not equilibrium lifetimes.",
            "Neutral field groups use paired streams; every group has independent replicas internally, but those three groups are not independent of one another.",
            "Fields +/-0.15 V/nm test driven response, not a claimed weak-field Einstein relation.",
        ],
    )


def _write_report(result, output):
    lines = [
        "# 虚拟数据测试实测报告", "",
        "**所有输入能量与势垒均为合成数据；本次没有进行真实 DFT/NEB 计算。**", "",
        f"- 条件数：{result['n_conditions']}；轨迹执行数：{result['n_trajectory_executions']}。",
        f"- 不重复的初始化/动力学子流对：{result['n_distinct_stream_pairs']}；中性场控制明确使用配对流。",
        f"- 实际事件总数：{result['total_events']}；全部完成指定物理时长：{result['all_trajectories_complete']}。",
        f"- 中性空位配对外场控制保持相同路径：{result['neutral_field_paired_invariance']}。",
        f"- 同种子完整事件轨迹重放核对：{result.get('sample_replay_verified')}。",
        "", "## 单空位扩散：与相同时间的确定性参考比较", "",
        "物理时长在抽样前按平衡平均跳率确定，所有重复运行到固定终止时间。初始化随机数与动力学随机数使用独立子流。",
        "矩方程参考给出同一时刻的 D(t)；长时间 D∞ 只用于考察有限时偏差。区间是否覆盖参考值仅作诊断。", "",
        "| 条件 | KMC / 精确 D(t) | 95% 区间覆盖参考 | 精确 D(t)/D∞ |",
        "|---|---:|---|---:|",
    ]
    for row in result["reference_comparisons"]:
        ratio = None if row["relative_error_vs_finite_time"] is None else 1+row["relative_error_vs_finite_time"]
        lines.append(f"| {row['condition_id']} | {'未提供' if ratio is None else f'{ratio:.4f}'} | {row['ci_contains_finite_time_reference']} | {1+row['finite_time_bias_vs_limit']:.5f} |")
    lines += ["", "## 环境消融和外场响应", "",
              "| 条件 | 平均事件数 | 平均空位对数 | 平均最大簇 | 平均 z 位移 / Å |",
              "|---|---:|---:|---:|---:|"]
    def fmt(value):
        return "未提供" if value is None else f"{value:.4g}"
    for group in result["groups"]:
        if group["condition"]["family"] == "diffusion":
            continue
        metrics = group["statistics"]["scalar_estimates"]
        lines.append(f"| {group['condition']['id']} | {fmt(metrics['events']['mean'])} | {fmt(metrics['pair_count']['mean'])} | {fmt(metrics['largest_cluster']['mean'])} | {fmt(group['z_displacement_A']['mean'])} |")
    lines += ["", "环境鞍点项在两个方向对称加入，能改变有限时动力学，但不改变同一构型能量模型的平衡分布。不同模型的瞬时簇大小排序不应预设。", "",
              "## 温度扫描的表观斜率", "",
              "下面是三个温度点的描述性斜率，未作为独立 NEB 势垒或实际材料激活能：", ""]
    for slope in result["apparent_arrhenius_slopes"]:
        lines.append(f"- q={slope['charge']}，{slope['source']}：{slope['apparent_activation_eV']:.5f} eV。")
    lines += ["", "## 文件与复现", "",
              "suite_plan.json 固化参数、时长和所有种子；replicas.jsonl 保存每条轨迹初态、终态、累计位移及通道计数；groups.csv 与 suite_summary.json 保存统计；sample_trajectory/ 保存一个实际条件的完整事件重放。",
              "原始虚拟计算输入为 mock_dft_input.json，导入后目录为 catalog.virtual.json。导入将正反势垒转换到同一鞍点，并检查闭合关系。",
              ""]
    if result.get("plot_requested"):
        lines += ["![虚拟数据测试总览](overview.png)", ""]
    lines += ["## 解释范围", ""]
    lines += ["- "+item for item in result["limitations"]]
    (output/"REPORT.md").write_text("\n".join(lines)+"\n", encoding="utf-8")


def plot_suite(result, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    colors = {0: "#2878ad", 2: "#dc7a24"}
    temperatures = sorted({g["condition"]["temperature_K"] for g in result["groups"] if g["condition"]["family"] == "diffusion"})
    diffusion_groups = [g for g in result["groups"] if g["condition"]["family"] == "diffusion"]
    hop_targets = sorted({g["condition"]["target_expected_hops"] for g in diffusion_groups})
    axis = axes[0, 0]
    for group in diffusion_groups:
        c, comparison = group["condition"], group["reference_comparison"]
        ratio = comparison["estimate_m2_s"]
        if ratio is None:
            continue
        reference = comparison["reference_finite_time_m2_s"]
        x = temperatures.index(c["temperature_K"])+(0 if c["charge"] == 0 else len(temperatures))
        offset = -.12 if c["target_expected_hops"] == hop_targets[0] else .12
        interval = comparison["estimate_ci_m2_s"]
        if interval:
            axis.vlines(x+offset, interval[0]/reference, interval[1]/reference, color=colors[c["charge"]])
        marker = "o" if c["target_expected_hops"] == hop_targets[0] else "s"
        axis.plot(x+offset, ratio/reference, marker, color=colors[c["charge"]], ms=5)
    axis.axhline(1, color="#555555", linestyle="--", linewidth=1)
    axis.set_xticks(range(2*len(temperatures)), [f"q={q}\n{t:g} K" for q in (0, 2) for t in temperatures])
    axis.set_ylabel("KMC D(t) / exact finite-time D(t)")
    axis.set_title("A  Numerical reference comparison", loc="left", fontweight="bold")
    axis.plot([], [], "o", color="#555555", label=f"{hop_targets[0]} expected hops")
    axis.plot([], [], "s", color="#555555", label=f"{hop_targets[-1]} expected hops")
    axis.legend(loc="upper right", frameon=False, fontsize=8)
    axis = axes[0, 1]
    for q in (0, 2):
        selected = sorted([g for g in diffusion_groups if g["condition"]["charge"] == q and
                           g["condition"]["target_expected_hops"] == hop_targets[-1]],
                          key=lambda g:g["condition"]["temperature_K"])
        ts = [g["condition"]["temperature_K"] for g in selected]
        exact = [g["reference_comparison"]["reference_long_time_m2_s"] for g in selected]
        estimates = [g["reference_comparison"]["estimate_m2_s"] for g in selected]
        axis.plot(ts, exact, "-", color=colors[q], label=f"q={q} exact long-time")
        axis.plot(ts, estimates, "s", color=colors[q])
    axis.set_yscale("log")
    axis.set_xlabel("Temperature (K)")
    axis.set_ylabel("Tracer D (m$^2$/s)")
    axis.set_title("B  Synthetic temperature response", loc="left", fontweight="bold")
    axis.plot([], [], "s", color="#555555", label="KMC finite-time estimate")
    axis.legend(loc="upper left", frameon=False, fontsize=8)
    axis = axes[1, 0]
    names = list(ABLATIONS)
    for q in (0, 2):
        for g in result["groups"]:
            c = g["condition"]
            if c["family"] != "mechanism" or c["charge"] != q:
                continue
            x = names.index(c["ablation"])+(-.1 if q == 0 else .1)
            metric = g["statistics"]["scalar_estimates"]["pair_count"]
            if metric["mean"] is None:
                continue
            if metric["ci"]:
                axis.vlines(x, *metric["ci"], color=colors[q])
            axis.plot(x, metric["mean"], "o", color=colors[q], label=f"q={q}" if x < .5 else None)
    axis.set_xticks(range(3), ["No pair / no TS", "Pair only", "Pair + TS"])
    axis.set_ylabel("Final vacancy pair count")
    axis.set_title("C  Equal-time environment controls", loc="left", fontweight="bold")
    axis.legend(loc="upper left", frameon=False, fontsize=8)
    axis = axes[1, 1]
    for q in (0, 2):
        groups = sorted([g for g in result["groups"] if g["condition"]["family"] == "field" and g["condition"]["charge"] == q],
                        key=lambda g:g["condition"]["field_V_nm"][2])
        for g in groups:
            field, metric = g["condition"]["field_V_nm"][2], g["z_displacement_A"]
            if metric["mean"] is not None:
                if metric["ci"]:
                    axis.vlines(field, *metric["ci"], color=colors[q])
                axis.plot(field, metric["mean"], "o", color=colors[q])
        axis.plot([], [], "o", color=colors[q], label=f"q={q}")
    axis.axhline(0, color="#777777", linestyle="--", linewidth=1)
    axis.set_xlabel("Applied field along z (V/nm)")
    axis.set_ylabel("Mean unwrapped z displacement (angstrom)")
    axis.set_title("D  Field polarity and neutral control", loc="left", fontweight="bold")
    axis.legend(frameon=False, fontsize=8)
    for axis in axes.flat:
        axis.grid(axis="y", alpha=.18)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("VIRTUAL DATA — NO DFT/NEB CALCULATIONS", fontweight="bold", fontsize=16)
    fig.text(.07, .035, "95% intervals use independent replicas within each group; neutral-field groups deliberately share paired streams.\n"
             "Finite-time CI coverage is diagnostic. Synthetic rates, cluster metrics and apparent slopes are not material predictions.",
             fontsize=9, color="#444444")
    fig.subplots_adjust(left=.09, right=.97, bottom=.14, top=.89, hspace=.4, wspace=.28)
    fig.savefig(Path(output)/"overview.png", dpi=170, facecolor="white")
    plt.close(fig)


def run_virtual_suite(catalog, output, profile="standard", workers=4, base_seed=20261111,
                      allow_illustrative=False, max_events=5000, mock_bundle=None, plot=True):
    _positive_integer(workers, "workers")
    catalog.validate(allow_illustrative=allow_illustrative)
    plan = build_suite_plan(catalog, profile, base_seed, max_events)
    directory = _new_output(output)
    catalog.write(directory/"catalog.virtual.json")
    _write_json(directory/"catalog.audit.json", catalog.audit())
    if mock_bundle is not None:
        from .mockdata import write_mock_bundle, import_mock_bundle
        if catalog_digest(import_mock_bundle(mock_bundle)) != catalog_digest(catalog):
            raise ValueError("mock bundle does not reproduce the supplied catalogue")
        write_mock_bundle(mock_bundle, directory/"mock_dft_input.json")
    _write_json(directory/"suite_plan.json", plan)
    jobs = [(c, stream) for c in plan["conditions"] for stream in c["replica_seeds"]]
    records = []
    start = time.monotonic()
    with (directory/"replicas.jsonl").open("w", encoding="utf-8") as handle:
        def record(value):
            records.append(value)
            handle.write(json.dumps(value, default=_json_default, allow_nan=False)+"\n")
            handle.flush()
            if len(records) % 16 == 0 or len(records) == len(jobs):
                _write_json(directory/"progress.json", dict(
                    completed_executions=len(records), total_executions=len(jobs),
                    completed_duration=sum(r["target_reached"] for r in records),
                    elapsed_seconds=time.monotonic()-start))
        if workers == 1:
            for c, stream in jobs:
                record(_simulate_replica(catalog, c, stream))
        else:
            with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                                     initializer=_worker_init, initargs=(catalog.data,)) as pool:
                futures = [pool.submit(_worker, job) for job in jobs]
                for future in as_completed(futures):
                    record(future.result())
    records.sort(key=lambda r:(r["condition_id"], r["replica_index"]))
    result = summarize_suite(plan, records)
    result.update(created_utc=datetime.now(timezone.utc).isoformat(), software_version=__version__,
                  elapsed_seconds=time.monotonic()-start, profile=profile,
                  catalog_sha256=catalog_digest(catalog), base_seed=base_seed,
                  mock_source=catalog.data["metadata"].get("mock_source"),
                  plot_requested=bool(plot),
                  **_graph_coverage(catalog))
    candidate = next((r for r in records if r["target_reached"] and
                      r["condition_id"] == "mechanism_q2_pair-environment"), None)
    if candidate is None:
        candidate = next((r for r in records if r["target_reached"]), None)
    result["sample_replay_verified"] = None
    if candidate is not None:
        condition = next(c for c in plan["conditions"] if c["id"] == candidate["condition_id"])
        config = SimulationConfig(
            repeats=tuple(condition["repeats"]), charge=condition["charge"],
            temperature_K=condition["temperature_K"], seed=candidate["seed"],
            field_V_nm=tuple(condition["field_V_nm"]), pair_scale=condition["pair_scale"],
            environment_scale=condition["environment_scale"])
        replay = run_simulation(catalog, config, duration_s=condition["duration_s"],
                                max_events=condition["max_events"], sample_every=10,
                                vacancy_indices=candidate["initial_vacancy_indices"],
                                output=directory/"sample_trajectory", allow_illustrative=True)
        result["sample_replay_verified"] = bool(
            replay["events"] == candidate["events"]
            and np.array_equal(replay["displacements_A"], candidate["displacements_A"]))
        result["sample_replay_record"] = dict(condition_id=candidate["condition_id"],
                                               replica_index=candidate["replica_index"])
    result["deterministic_checks_passed"] &= result["sample_replay_verified"] is True
    rows = []
    for group in result["groups"]:
        c, stats = group["condition"], group["statistics"]
        row = {key:c.get(key) for key in ("id", "family", "charge", "temperature_K", "field_V_nm", "duration_s", "target_expected_hops")}
        row.update(n_complete=stats["n_complete_replicas"], n_incomplete=stats["n_incomplete_replicas"],
                   z_displacement_mean_A=group["z_displacement_A"]["mean"],
                   z_displacement_ci_A=group["z_displacement_A"]["ci"])
        for metric, value in stats["scalar_estimates"].items():
            row[metric+"_mean"], row[metric+"_ci"] = value["mean"], value["ci"]
        row.update(group.get("reference_comparison", {}))
        rows.append(row)
    _write_csv(directory/"groups.csv", rows)
    _write_csv(directory/"replicas.csv", [
        {key:r.get(key) for key in ("condition_id", "replica_index", "seed", "initialization_seed", "status", "time_s", "duration_s", "events", "initial_vacancy_indices", "displacements_A", "error")}
        for r in records])
    _write_json(directory/"suite_summary.json", result)
    _write_report(result, directory)
    if plot:
        plot_suite(result, directory)
    (directory/"progress.json").unlink(missing_ok=True)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--mock-input", help="Existing explicitly synthetic raw calculation bundle; default generates the built-in virtual dataset.")
    parser.add_argument("--allow-illustrative", action="store_true")
    parser.add_argument("--profile", choices=tuple(PROFILES), default="standard")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument("--max-events", type=int, default=5000)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)
    try:
        if not args.allow_illustrative:
            raise ValueError("Virtual data require explicit --allow-illustrative")
        from .mockdata import make_mock_bundle, load_mock_bundle, import_mock_bundle
        bundle = load_mock_bundle(args.mock_input) if args.mock_input else make_mock_bundle()
        result = run_virtual_suite(
            import_mock_bundle(bundle), args.output, profile=args.profile, workers=args.workers,
            base_seed=args.seed, allow_illustrative=True, max_events=args.max_events,
            mock_bundle=bundle, plot=not args.no_plot)
        print(json.dumps({key:result[key] for key in (
            "n_conditions", "n_trajectory_executions", "n_distinct_stream_pairs",
            "total_events", "all_trajectories_complete", "neutral_field_paired_invariance",
            "sample_replay_verified", "deterministic_checks_passed", "statistical_diagnostics",
            "elapsed_seconds")}, indent=2))
        return 0 if result["deterministic_checks_passed"] else 3
    except (ValueError, OSError, RuntimeError) as exc:
        print("virtual_suite: "+str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
