"""Deterministic references for one vacancy in a declared periodic hop graph.

These solve the *supplied model*, including synthetic parameters and any stated
channel exclusions.  They do not validate the model as a material prediction.
Only zero-field, isothermal, noninteracting single-vacancy motion is covered.

The internal state is the eight-site crystallographic basis.  A hop also carries
an unwrapped Cartesian displacement.  Thus a hop to another periodic image of
the same basis state changes position and jump count even though it leaves the
basis-state generator unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from numbers import Integral, Real
from typing import Any

import numpy as np
from scipy.linalg import expm

from .catalog import Catalog
from .model import KB


_BALANCE_TOL = 1e-10
_SOLVE_TOL = 1e-10
_MOMENT_TOL = 1e-9
_CONDITION_LIMIT = 1e12
_SECOND_COMPONENTS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


@dataclass
class _PeriodicModel:
    generator: np.ndarray  # Row generator, divided by the common rate scale.
    first_derivatives: np.ndarray  # Three tilted-generator derivatives / scale.
    second_derivatives: np.ndarray  # 3 x 3 tilted-generator derivatives / scale.
    probabilities: np.ndarray
    sources: np.ndarray
    targets: np.ndarray
    displacements_A: np.ndarray
    scaled_rates: np.ndarray
    rate_scale_Hz: float
    local_drift_scaled: np.ndarray
    mean_jump_rate_Hz: float
    drift_A_s: np.ndarray
    length_scale_A: float
    diagnostics: dict[str, Any]
    metadata: dict[str, Any]


def _components(n: int, sources: np.ndarray, targets: np.ndarray) -> list[list[int]]:
    adjacency = [set() for _ in range(n)]
    for source, target in zip(sources, targets):
        if source != target:
            adjacency[int(source)].add(int(target))
            adjacency[int(target)].add(int(source))
    remaining = set(range(n))
    components = []
    while remaining:
        stack = [min(remaining)]
        remaining.remove(stack[0])
        component = []
        while stack:
            node = stack.pop()
            component.append(node)
            neighbors = adjacency[node] & remaining
            remaining -= neighbors
            stack.extend(neighbors)
        components.append(sorted(component))
    return components


def _prepare_model(catalog: Catalog, charge: int, temperature_K: float) -> _PeriodicModel:
    if not isinstance(charge, Integral) or isinstance(charge, (bool, np.bool_)):
        raise ValueError("charge must be an integer fixed-charge branch.")
    if (not isinstance(temperature_K, Real) or isinstance(temperature_K, (bool, np.bool_))
            or not math.isfinite(float(temperature_K)) or temperature_K <= 0):
        raise ValueError("temperature_K must be finite and positive.")
    # A synthetic catalogue is legitimate input to an analytic software benchmark;
    # its declared status is retained in every result. Pending/inconsistent data fail.
    catalog.validate(allow_illustrative=True)
    energies = np.asarray(catalog.site_energies(charge), dtype=float)
    n = len(energies)
    if n != 8:
        raise ValueError("This reference currently requires exactly eight basis states.")
    thermal_energy = KB * float(temperature_K)
    log_weights = -(energies - energies.min()) / thermal_energy
    weights = np.exp(log_weights)
    if not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Equilibrium basis probabilities exceed the floating-point dynamic range.")
    probabilities = weights / weights.sum()
    log_probabilities = log_weights - math.log(float(weights.sum()))
    sources, targets, displacements, log_rates = [], [], [], []
    for channel in catalog.data["channels"]:
        source, target = int(channel["from_basis"]), int(channel["to_basis"])
        displacement = np.asarray(channel["displacement_A"], dtype=float)
        parameter = catalog.channel(channel["id"], charge)
        saddle, attempt = float(parameter["saddle_eV"]), float(parameter["attempt_Hz"])
        for start, finish, vector in ((source, target, displacement), (target, source, -displacement)):
            sources.append(start)
            targets.append(finish)
            displacements.append(vector)
            log_rates.append(math.log(attempt) - (saddle - energies[start]) / thermal_energy)
    sources = np.asarray(sources, dtype=int)
    targets = np.asarray(targets, dtype=int)
    displacements = np.asarray(displacements, dtype=float)
    log_rates = np.asarray(log_rates, dtype=float)
    if not len(log_rates) or not np.isfinite(log_rates).all():
        raise ValueError("The reference requires finite positive directed-event rates.")
    log_scale = float(log_rates.max())
    try:
        rate_scale = math.exp(log_scale)
    except OverflowError as exc:
        raise ValueError("Event rates exceed the floating-point dynamic range.") from exc
    scaled_rates = np.exp(log_rates - log_scale)
    if rate_scale <= 0 or not math.isfinite(rate_scale) or np.any(scaled_rates <= 0):
        raise ValueError("Event-rate underflow would change the graph; this numerical reference cannot resolve it.")
    components = _components(n, sources, targets)
    if len(components) != 1:
        raise ValueError(f"Disconnected basis transition graph: components {components}. This reference requires a unique equilibrium basis chain.")

    generator = np.zeros((n, n))
    first = np.zeros((3, n, n))
    second = np.zeros((3, 3, n, n))
    outgoing = np.zeros(n)
    for source, target, displacement, rate in zip(sources, targets, displacements, scaled_rates):
        # Self-basis translations cancel in Q(0), but NOT in its derivatives.
        if source != target:
            generator[source, target] += rate
            generator[source, source] -= rate
        first[:, source, target] += rate * displacement
        second[:, :, source, target] += rate * np.outer(displacement, displacement)
        outgoing[source] += rate
    local_drift = first.sum(axis=2).T
    mean_rate_scaled = float(probabilities @ outgoing)
    drift_scaled = probabilities @ local_drift
    mean_jump_rate = mean_rate_scaled * rate_scale
    drift = drift_scaled * rate_scale
    if not math.isfinite(mean_jump_rate) or not np.isfinite(drift).all():
        raise ValueError("Physical rate/drift output exceeds the floating-point dynamic range.")
    length_scale = float(np.max(np.linalg.norm(displacements, axis=1)))
    generator_norm = float(np.linalg.norm(generator, ord=np.inf))
    row_sum_error = float(np.max(np.abs(generator.sum(axis=1)))) / generator_norm
    stationary_error = float(np.max(np.abs(probabilities @ generator))) / generator_norm
    # Every canonical path contributes two adjacent directed entries.
    log_fluxes = log_probabilities[sources] + log_rates
    log_flux_ratio = log_fluxes[0::2] - log_fluxes[1::2]
    balance_error = float(np.max(np.abs(log_flux_ratio)))
    drift_error = float(np.linalg.norm(drift_scaled)) / (mean_rate_scaled * length_scale)
    if max(row_sum_error, stationary_error, balance_error, drift_error) > _BALANCE_TOL:
        raise ValueError("Zero-field detailed balance, stationarity or zero-drift check failed; no reference is reported.")

    original_metadata = catalog.data["metadata"]
    exclusions = dict(original_metadata.get("channel_exclusions", {}))
    cutoff = original_metadata.get("migration_cutoff_A")
    if cutoff is None:
        cutoff = max(float(channel["length_A"]) for channel in catalog.data["channels"]) + 1e-8
    coverage_warning = (
        "This reference is conditional on the declared truncated migration graph; excluded paths do not contribute and full material diffusion is not established."
        if exclusions else ""
    )
    canonical = json.dumps(catalog.data, sort_keys=True, separators=(",", ":"), allow_nan=False)
    warnings = [coverage_warning] if exclusions else []
    if original_metadata["calibration_status"] == "illustrative":
        warnings.append("Synthetic illustrative parameters: this analytic model result is not DFT-calibrated material data.")
    metadata = {
        "schema_version": "dftkmc-reference-1", "charge": int(charge),
        "temperature_K": float(temperature_K), "units": "m^2/s",
        "structure_id": catalog.structure_id,
        "catalog_sha256": hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        "calibration_status": original_metadata["calibration_status"],
        "catalog_provenance": original_metadata.get("provenance"),
        "migration_cutoff_A": float(cutoff), "channel_exclusions": exclusions,
        "conditional_on_truncated_graph": bool(exclusions), "coverage_warning": coverage_warning,
        "basis_probabilities": probabilities.tolist(), "mean_jump_rate_Hz": mean_jump_rate,
        "drift_A_s": drift.tolist(), "warnings": warnings,
        "assumptions": [
            "One vacancy, a fixed charge branch, constant temperature and exactly zero applied field.",
            "Infinite periodic continuation of the declared geometric hop catalogue; no vacancy-vacancy interactions or environmental saddle shifts.",
            "Each canonical path is included in both directions with its actual unwrapped Cartesian displacement and one shared attempt prefactor.",
            "Classical rates k_ij = attempt_Hz * exp(-(saddle_eV - site_energy_i)/(k_B T)); equilibrium basis probabilities are proportional to exp(-site_energy_i/(k_B T)).",
            "The eight-state basis graph is connected. This alone does not guarantee unbounded transport in all three spatial directions.",
            "An analytic solution of supplied parameters tests KMC implementation; it does not establish those parameters as physically calibrated.",
        ],
    }
    diagnostics = {
        "basis_state_count": n, "basis_connected_components": components,
        "canonical_channel_count": len(catalog.data["channels"]), "directed_event_count": len(sources),
        "self_basis_directed_event_count": int(np.count_nonzero(sources == targets)),
        "rate_scale_Hz": rate_scale,
        "log_rate_dynamic_range": float(log_rates.max() - log_rates.min()),
        "generator_row_sum_relative_residual": row_sum_error,
        "stationary_probability_relative_residual": stationary_error,
        "detailed_balance_max_abs_log_flux_ratio": balance_error,
        "zero_field_drift_relative_residual": drift_error,
        "balance_tolerance": _BALANCE_TOL,
    }
    return _PeriodicModel(generator, first, second, probabilities, sources, targets,
                          displacements, scaled_rates, rate_scale, local_drift,
                          mean_jump_rate, drift, length_scale, diagnostics, metadata)


def _check_tensor(tensor: np.ndarray, name: str, tolerance: float) -> dict[str, float]:
    if not np.isfinite(tensor).all():
        raise ValueError(f"{name} contains non-finite values.")
    scale = float(np.linalg.norm(tensor, ord=2))
    symmetry_error = float(np.max(np.abs(tensor - tensor.T))) / scale if scale else 0.0
    minimum = float(np.linalg.eigvalsh((tensor + tensor.T) / 2.0).min())
    negative_relative = max(0.0, -minimum / scale) if scale else 0.0
    if symmetry_error > tolerance or negative_relative > tolerance:
        raise ValueError(f"{name} failed symmetry/positive-semidefinite numerical checks.")
    return {"tensor_symmetry_relative_residual": symmetry_error,
            "tensor_minimum_eigenvalue": minimum,
            "tensor_negative_eigenvalue_relative_residual": negative_relative}


def periodic_tracer_reference(catalog: Catalog, charge: int, temperature_K: float) -> dict[str, Any]:
    """Long-time D of the declared periodic, noninteracting single-vacancy model.

    Solve Q chi = -b with pi chi = 0, then average corrected hop outer
    products.  All rates are scaled by a common factor for the linear solve.
    Disconnected or numerically ill-conditioned basis chains fail explicitly.
    Synthetic catalogues are permitted and remain labelled as synthetic.
    """
    model = _prepare_model(catalog, charge, temperature_K)
    n = len(model.probabilities)
    q_norm = float(np.linalg.norm(model.generator, ord=np.inf))
    normalized_q = model.generator / q_norm
    normalized_b = model.local_drift_scaled / q_norm
    bordered = np.zeros((n + 1, n + 1))
    bordered[:n, :n] = normalized_q
    bordered[:n, n] = 1.0
    bordered[n, :n] = model.probabilities
    condition = float(np.linalg.cond(bordered))
    if not math.isfinite(condition) or condition > _CONDITION_LIMIT:
        raise ValueError(f"The connected basis corrector is numerically singular or ill-conditioned (condition number {condition:.3g}); no reliable double-precision reference is reported.")
    rhs = np.zeros((n + 1, 3))
    rhs[:n] = -normalized_b
    try:
        solution = np.linalg.solve(bordered, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("The basis corrector linear system is singular.") from exc
    corrector = solution[:n]
    equation_error = normalized_q @ corrector + normalized_b
    denominator = float(np.max(np.abs(corrector))) + float(np.max(np.abs(normalized_b)))
    equation_relative = float(np.max(np.abs(equation_error))) / (denominator or model.length_scale_A)
    gauge_A = float(np.max(np.abs(model.probabilities @ corrector)))
    gauge_relative = gauge_A / max(model.length_scale_A, float(np.max(np.abs(corrector))))
    if max(equation_relative, gauge_relative) > _SOLVE_TOL:
        raise ValueError("Cell-corrector equation or equilibrium-weighted gauge residual exceeds tolerance.")
    corrected = model.displacements_A + corrector[model.targets] - corrector[model.sources]
    flux = model.probabilities[model.sources] * model.scaled_rates
    gram = 0.5 * np.einsum("n,ni,nj->ij", flux, corrected, corrected)
    tensor = gram * (model.rate_scale_Hz * 1e-20)
    tensor = (tensor + tensor.T) / 2.0
    diagnostics = dict(model.diagnostics)
    diagnostics.update({
        "corrector_bordered_condition_number": condition, "corrector_condition_limit": _CONDITION_LIMIT,
        "corrector_equation_relative_residual": equation_relative,
        "corrector_equation_max_abs_residual_A_s": float(np.max(np.abs(equation_error))) * q_norm * model.rate_scale_Hz,
        "corrector_gauge_max_abs_residual_A": gauge_A,
        "corrector_gauge_relative_residual": gauge_relative,
        "corrector_solvability_multiplier_scaled": solution[n].tolist(),
        "corrector_solve_tolerance": _SOLVE_TOL,
        **_check_tensor(tensor, "Long-time diffusion tensor", _SOLVE_TOL),
    })
    label = "long-time periodic single-vacancy tracer reference"
    if model.metadata["conditional_on_truncated_graph"]:
        label += " conditional on the declared truncated graph"
    return {**model.metadata, "reference_kind": "periodic_long_time_single_vacancy", "label": label,
            "limit": "t -> infinity", "tensor_m2_s": tensor.tolist(),
            "trace_over_3_m2_s": float(np.trace(tensor) / 3.0),
            "corrector_A": corrector.tolist(), "residual_diagnostics": diagnostics,
            "method": "Eight-state row generator and equilibrium-gauged cell corrector; D = (1/2) sum_i,j pi_i k_ij (d_ij + chi_j - chi_i)(...)^T.",
            "finite_time_warning": "Equilibrium initial basis does not imply that covariance/(2t) already equals this long-time limit. Use finite_time_tracer_reference at the observation time."}


def finite_time_tracer_reference(
    catalog: Catalog, charge: int, temperature_K: float, duration_s: float,
) -> dict[str, Any]:
    """Exact finite-time displacement moments from an equilibrium basis start.

    Propagate 80 linear variables: 8 state probabilities, 3 x 8 first moments,
    and 6 x 8 independent second moments.  For row-vector moments,

        p' = p Q
        m_a' = m_a Q + p Q_a
        M_ab' = M_ab Q + m_a Q_b + m_b Q_a + p Q_ab.

    Tilted-generator derivatives Q_a/Q_ab include self-basis translations.
    Initial displacement is zero; initial basis probabilities are Boltzmann.
    The returned tensor is the centered covariance/(2t), not an assumed limit.
    """
    if (not isinstance(duration_s, Real) or isinstance(duration_s, (bool, np.bool_))
            or not math.isfinite(float(duration_s)) or duration_s <= 0):
        raise ValueError("duration_s must be finite and positive.")
    model = _prepare_model(catalog, charge, temperature_K)
    n = len(model.probabilities)
    n_blocks = 1 + 3 + len(_SECOND_COMPONENTS)
    matrix = np.zeros((n * n_blocks, n * n_blocks))

    def block(index: int) -> slice:
        return slice(index * n, (index + 1) * n)

    for index in range(n_blocks):
        matrix[block(index), block(index)] = model.generator
    for axis in range(3):
        matrix[block(0), block(1 + axis)] = model.first_derivatives[axis]
    for component, (first_axis, second_axis) in enumerate(_SECOND_COMPONENTS, start=4):
        matrix[block(0), block(component)] = model.second_derivatives[first_axis, second_axis]
        matrix[block(1 + first_axis), block(component)] += model.first_derivatives[second_axis]
        matrix[block(1 + second_axis), block(component)] += model.first_derivatives[first_axis]
    scaled_time = float(duration_s) * model.rate_scale_Hz
    if not math.isfinite(scaled_time) or scaled_time <= 0:
        raise ValueError("Time-rate product exceeds the floating-point range for finite-time propagation.")
    initial = np.zeros(n * n_blocks)
    initial[block(0)] = model.probabilities
    with np.errstate(over="ignore", invalid="ignore"):
        final = initial @ expm(matrix * scaled_time)
    if not np.isfinite(final).all():
        raise ValueError("Finite-time matrix exponential produced non-finite moments.")
    final_probabilities = final[block(0)]
    mean = np.asarray([final[block(1 + axis)].sum() for axis in range(3)])
    raw_second = np.zeros((3, 3))
    for component, (first_axis, second_axis) in enumerate(_SECOND_COMPONENTS, start=4):
        raw_second[first_axis, second_axis] = raw_second[second_axis, first_axis] = final[block(component)].sum()
    covariance = raw_second - np.outer(mean, mean)
    covariance = (covariance + covariance.T) / 2.0
    tensor = covariance * (1e-20 / (2.0 * float(duration_s)))
    probability_error = float(np.max(np.abs(final_probabilities - model.probabilities)))
    mass_error = abs(float(final_probabilities.sum()) - 1.0)
    displacement_scale = max(model.length_scale_A, math.sqrt(max(float(np.trace(raw_second)), 0.0)))
    mean_error = float(np.linalg.norm(mean - model.drift_A_s * float(duration_s))) / displacement_scale
    if max(probability_error, mass_error, mean_error) > _MOMENT_TOL or final_probabilities.min() < -_MOMENT_TOL:
        raise ValueError("Finite-time probability, stationarity or zero-field mean-moment residual exceeds tolerance.")
    diagnostics = dict(model.diagnostics)
    diagnostics.update({
        "moment_system_dimension": n * n_blocks,
        "moment_block_order": ["p", "m_x", "m_y", "m_z", "M_xx", "M_yy", "M_zz", "M_xy", "M_xz", "M_yz"],
        "probability_mass_residual": mass_error,
        "equilibrium_probability_max_abs_residual": probability_error,
        "minimum_final_basis_probability": float(final_probabilities.min()),
        "mean_displacement_relative_residual": mean_error,
        "moment_check_tolerance": _MOMENT_TOL,
        **_check_tensor(covariance, "Finite-time displacement covariance", _MOMENT_TOL),
    })
    label = "finite-time equilibrium single-vacancy tracer reference"
    if model.metadata["conditional_on_truncated_graph"]:
        label += " conditional on the declared truncated graph"
    return {**model.metadata, "reference_kind": "finite_time_equilibrium_single_vacancy", "label": label,
            "duration_s": float(duration_s), "tensor_m2_s": tensor.tolist(),
            "trace_over_3_m2_s": float(np.trace(tensor) / 3.0),
            "mean_displacement_A": mean.tolist(), "raw_second_moment_A2": raw_second.tolist(),
            "covariance_A2": covariance.tolist(), "final_basis_probabilities": final_probabilities.tolist(),
            "expected_hop_count": model.mean_jump_rate_Hz * float(duration_s),
            "initial_condition": "Equilibrium basis probabilities, zero initial unwrapped displacement.",
            "method": "Exact linear Markov-additive probability/first/second moment system evaluated with scipy.linalg.expm (80 variables); tensor = centered covariance/(2t).",
            "residual_diagnostics": diagnostics,
            "finite_time_warning": "This is the model reference at the specified observation time. It need not equal the long-time periodic tracer tensor."}
