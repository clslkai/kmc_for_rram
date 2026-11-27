"""Replica-level analysis for fixed-time, fixed-charge KMC experiments.

Time samples and vacancies in one trajectory are deliberately *not* treated as
independent observations.  Diffusion output is a finite-time tracer estimate;
its physical long-time convergence must be established in a separate study.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from scipy.stats import t as student_t


SCALAR_METRICS = ("energy_eV", "msd_A2", "pair_count", "largest_cluster", "events")


def reached_target(time_s: float, duration_s: float) -> bool:
    """Use a relative tolerance, without a large absolute tolerance for short runs."""
    return bool(np.isfinite(time_s) and abs(time_s - duration_s) <= duration_s * 1e-10)


def _scalar_summary(values: Sequence[float], confidence: float, infer: bool) -> dict[str, Any]:
    a = np.asarray(values, dtype=float)
    n = len(a)
    if not n:
        return {"n_replicas": 0, "mean": None, "std": None, "ci": None}
    mean = float(a.mean())
    sd = float(a.std(ddof=1)) if n > 1 else None
    interval = None
    if infer and n > 1:
        halfwidth = float(student_t.ppf((1.0 + confidence) / 2.0, n - 1) * sd / np.sqrt(n))
        interval = [mean - halfwidth, mean + halfwidth]
    return {
        "n_replicas": n, "mean": mean, "std": sd, "ci": interval,
        "ci_method": "Student t interval over independent replica means" if interval else None,
    }


def _moments(displacements: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    means, seconds = [], []
    for item in displacements:
        a = np.asarray(item, dtype=float)
        if a.ndim != 2 or a.shape[1] != 3 or not len(a) or not np.isfinite(a).all():
            raise ValueError("Each replica must contain finite unwrapped displacements of shape (n_vacancies, 3).")
        means.append(a.mean(axis=0))
        seconds.append(np.einsum("ni,nj->ij", a, a) / len(a))
    return np.asarray(means), np.asarray(seconds)


def _centered_tracer_tensor(means: np.ndarray, seconds: np.ndarray) -> np.ndarray:
    """Unbiased estimate of E[dd^T] - E[d]E[d]^T, grouped by replica.

    For several correlated vacancy tracers, uncertainty in the ensemble drift
    comes from the replica mean, not the number of vacancies.  The correction
    is Cov(replica mean)/R.  For one vacancy this reduces to sample covariance.
    """
    n = len(means)
    if n < 2:
        raise ValueError("An ensemble-centered tracer tensor needs at least two independent replicas.")
    drift = means.mean(axis=0)
    deviations = means - drift
    drift_covariance = deviations.T @ deviations / (n - 1)
    tensor = seconds.mean(axis=0) - np.outer(drift, drift) + drift_covariance / n
    return (tensor + tensor.T) * 0.5


def finite_time_tracer_estimate(
    displacements: Sequence[np.ndarray], duration_s: float, *,
    confidence: float = 0.95, bootstrap_samples: int = 2000,
    bootstrap_seed: int = 731,
) -> dict[str, Any]:
    """Estimate D_ij(t) = centered <dr_i dr_j> / (2t), in m^2/s.

    The caller must establish bulk periodic geometry, zero field, identical
    physical observation time, and absence of censored replicas.
    """
    if not np.isfinite(duration_s) or duration_s <= 0:
        raise ValueError("duration_s must be finite and positive.")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one.")
    if bootstrap_samples < 100:
        raise ValueError("bootstrap_samples must be at least 100.")
    means, seconds = _moments(displacements)
    n = len(means)
    if n < 2:
        return {"available": False, "reason": "At least two independent replicas are needed to estimate ensemble-centered displacement."}
    covariance = _centered_tracer_tensor(means, seconds)
    conversion = 1e-20 / (2.0 * duration_s)
    tensor = covariance * conversion
    rng = np.random.default_rng(bootstrap_seed)
    samples = np.empty((bootstrap_samples, 3, 3), dtype=float)
    for index in range(bootstrap_samples):
        selected = rng.integers(0, n, size=n)
        samples[index] = _centered_tracer_tensor(means[selected], seconds[selected]) * conversion
    tail = (1.0 - confidence) / 2.0
    bounds = np.quantile(samples, [tail, 1.0 - tail], axis=0)
    isotropic_samples = np.trace(samples, axis1=1, axis2=2) / 3.0
    counts = [int(len(x)) for x in displacements]
    uncertainty_resolved = bool(np.any(np.ptp(means, axis=0) != 0.0) or np.any(np.ptp(seconds, axis=0) != 0.0))
    return {
        "available": True,
        "label": "finite-time tracer estimate",
        "definition": "D_ij(t) = ensemble-centered unwrapped tracer displacement covariance_ij / (2 t)",
        "time_s": float(duration_s), "units": "m^2/s",
        "n_independent_replicas": n, "vacancies_per_replica": counts,
        "tensor_m2_s": tensor.tolist(),
        "tensor_ci_low_m2_s": bounds[0].tolist() if uncertainty_resolved else None,
        "tensor_ci_high_m2_s": bounds[1].tolist() if uncertainty_resolved else None,
        "trace_over_3_m2_s": float(np.trace(tensor) / 3.0),
        "trace_over_3_ci_m2_s": np.quantile(isotropic_samples, [tail, 1.0 - tail]).tolist() if uncertainty_resolved else None,
        "uncertainty_resolved": uncertainty_resolved,
        "mean_displacement_A": means.mean(axis=0).tolist(),
        "centered_displacement_tensor_A2": covariance.tolist(),
        "centering_correction": "Covariance of independent replica mean displacements divided by replica count; accounts for estimated ensemble drift.",
        "confidence": confidence,
        "ci_method": "Percentile bootstrap resampling entire independent replicas; all vacancy tracers in a replica remain together.",
        "bootstrap_samples": int(bootstrap_samples), "bootstrap_seed": int(bootstrap_seed),
        "limitations": [
            "A single observation time does not establish the diffusive long-time regime or a converged diffusion coefficient.",
            "This is tracer displacement, not collective chemical diffusion or electrical conductivity.",
            "Vacancy tracers within a replica may be correlated; they are not separate uncertainty samples.",
        ] + (["Fewer than eight replicas: bootstrap intervals can be unstable; add independent replicas."] if n < 8 else [])
          + ([] if uncertainty_resolved else ["All replicas have identical first and second displacement moments. Empirical resampling cannot establish uncertainty; intervals are withheld. Zero observed motion does not bound unobserved rare-event diffusivity."]),
    }


def summarize_replicas(
    replicas: Sequence[Mapping[str, Any]], *, duration_s: float,
    boundary: str = "bulk", field_V_nm: Sequence[float] = (0.0, 0.0, 0.0),
    confidence: float = 0.95, bootstrap_samples: int = 2000,
    bootstrap_seed: int = 731,
) -> dict[str, Any]:
    """Summarize final states at one requested physical time.

    Each input has ``time_s``, ``status``, ``final`` (snapshot), and
    ``displacements_A``.  Failed runs may omit displacements.  If any intended
    replica is incomplete, headline inference is withheld: selecting only runs
    that stayed below an event cap can introduce kinetic selection bias.
    """
    if not np.isfinite(duration_s) or duration_s <= 0:
        raise ValueError("duration_s must be finite and positive.")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one.")
    if not replicas:
        raise ValueError("At least one replica is required.")
    field = np.asarray(field_V_nm, dtype=float)
    if field.shape != (3,) or not np.isfinite(field).all():
        raise ValueError("field_V_nm must contain three finite numbers.")
    completed = [r for r in replicas if not r.get("error") and reached_target(float(r.get("time_s", 0.0)), duration_s)]
    # Avoid equality comparisons of dictionaries containing NumPy arrays.
    complete_ids = {id(r) for r in completed}
    incomplete = [r for r in replicas if id(r) not in complete_ids]
    all_complete = not incomplete
    status_counts: dict[str, int] = {}
    for r in replicas:
        status = str(r.get("status", "unknown"))
        status_counts[status] = status_counts.get(status, 0) + 1
    scalars = {}
    for metric in SCALAR_METRICS:
        values = [float(r.get("final", {}).get(metric, r.get("events", 0))) for r in completed
                  if metric in r.get("final", {}) or (metric == "events" and "events" in r)]
        scalars[metric] = _scalar_summary(values, confidence, all_complete)
    warnings = []
    if incomplete:
        warnings.append("Some intended replicas did not reach the requested physical time. Completed-only scalar values are descriptive and conditional on completion; diffusion and inferential confidence intervals are withheld to avoid event-cap selection bias.")
    if len(completed) < 2:
        warnings.append("Fewer than two complete independent replicas: uncertainty cannot be estimated.")
    result: dict[str, Any] = {
        "requested_time_s": float(duration_s),
        "n_requested_replicas": len(replicas), "n_complete_replicas": len(completed),
        "n_incomplete_replicas": len(incomplete), "all_replicas_complete": all_complete,
        "termination_counts": status_counts, "confidence": confidence,
        "statistical_unit": "one independent simulation replica, not time samples or individual vacancies",
        "scalar_estimates": scalars,
        "incomplete_replicas": [{"seed": r.get("seed"), "status": r.get("status"),
                                  "time_s": r.get("time_s"), "error": r.get("error")} for r in incomplete],
        "warnings": warnings,
    }
    if not all_complete:
        diffusion = {"available": False, "reason": "Incomplete replica set; rerun with a larger event cap or a shorter common physical duration."}
    elif boundary != "bulk":
        diffusion = {"available": False, "reason": "Tracer diffusion is reported only for bulk periodic geometry."}
    elif np.any(field != 0.0):
        diffusion = {"available": False, "reason": "Tracer diffusion is reported only at exactly zero applied field; finite-field drift is not interpreted as equilibrium diffusion."}
    elif len(completed) < 2:
        diffusion = {"available": False, "reason": "At least two independent replicas are needed for ensemble-centered displacement."}
    else:
        diffusion = finite_time_tracer_estimate(
            [np.asarray(r["displacements_A"], dtype=float) for r in completed], duration_s,
            confidence=confidence, bootstrap_samples=bootstrap_samples, bootstrap_seed=bootstrap_seed,
        )
    result["diffusion"] = diffusion
    return result


def compare_groups(
    first: Sequence[Mapping[str, Any]], second: Sequence[Mapping[str, Any]], *,
    duration_s: float, confidence: float = 0.95, metrics: Sequence[str] = SCALAR_METRICS,
) -> dict[str, Any]:
    """Independent-group differences (second minus first) at a common time."""
    for group in (first, second):
        if not group or any(r.get("error") or not reached_target(float(r.get("time_s", 0)), duration_s) for r in group):
            return {"available": False, "reason": "A contrast requires all intended replicas to reach the common physical time."}
    differences: dict[str, Any] = {}
    for metric in metrics:
        if not all(metric in r.get("final", {}) for r in (*first, *second)):
            continue
        a = np.asarray([r["final"][metric] for r in first], dtype=float)
        b = np.asarray([r["final"][metric] for r in second], dtype=float)
        difference = float(b.mean() - a.mean())
        interval = None
        degrees = None
        if len(a) > 1 and len(b) > 1:
            va = float(a.var(ddof=1) / len(a))
            vb = float(b.var(ddof=1) / len(b))
            variance = va + vb
            if variance == 0.0:
                interval = [difference, difference]
            else:
                degrees = variance ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
                halfwidth = float(student_t.ppf((1.0 + confidence) / 2.0, degrees) * np.sqrt(variance))
                interval = [difference - halfwidth, difference + halfwidth]
        differences[metric] = {"mean_difference": difference, "ci": interval, "degrees_of_freedom": degrees}
    return {"available": True, "orientation": "second minus first", "confidence": confidence,
            "ci_method": "Welch Student t interval for independent replica groups; no multiple-comparison adjustment",
            "differences": differences}
