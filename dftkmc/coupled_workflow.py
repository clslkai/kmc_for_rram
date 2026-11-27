"""Auditable file outputs for exchange/transport device trajectories."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import numpy as np

from .coupled import CoupledConfig, CoupledDeviceSimulator, FIELD_CLOSURE
from .workflow import _graph_coverage, _new_output, _write_csv, _write_json, catalog_digest


def run_coupled(catalog, config: CoupledConfig, *, duration_s=None, pulses=None,
                n_vacancies=0, vacancy_indices=None, max_events=20000, sample_every=10,
                output=None, plot=False, allow_illustrative=False):
    """Run a declared model, preserving incomplete and model-domain outcomes.

    Exactly one of duration_s and pulses is required. Pulse voltage means the
    requested top-minus-bottom bias; the bottom reference remains fixed.
    """
    if (duration_s is None) == (pulses is None):
        raise ValueError("provide exactly one of duration_s and pulses")
    catalog.validate(allow_illustrative=allow_illustrative)
    if not allow_illustrative:
        raise ValueError("exchange and effective transport are illustrative; explicit --allow-illustrative is required")
    simulator = CoupledDeviceSimulator(catalog, config)
    simulator.initialize(n_vacancies=n_vacancies, vacancy_indices=vacancy_indices)
    initial = simulator.snapshot()
    initial_oxygen_ids = simulator.oxygen_ids.copy()
    initial_vacancy_ids = simulator.vacancy_ids.copy()
    initial_pool_ids = {name: list(pool) for name, pool in simulator.reservoir_oxygen_ids.items()}
    directory = _new_output(output) if output is not None else None
    if pulses is None:
        trajectory = simulator.run(duration_s, max_events=max_events, sample_every=sample_every)
    else:
        trajectory = simulator.run_waveform(pulses, max_events=max_events, sample_every=sample_every)
    records, events = trajectory["records"], trajectory["events"]
    initial = dict(records[0])  # includes the first segment's actual/compliance-limited bias
    final = simulator.snapshot()
    coverage = _graph_coverage(catalog)
    digest = catalog_digest(catalog)
    geometry = simulator.geometry_boundary
    summary = dict(
        schema_version="dftkmc-coupled-run-1", code_version="0.4.0",
        created_utc=datetime.now(timezone.utc).isoformat(),
        calibration_status="illustrative", migration_calibration_status=catalog.data["metadata"]["calibration_status"],
        catalog_sha256=digest, structure_id=catalog.structure_id,
        config=config.to_dict(), **coverage,
        kinetic_field_closure=FIELD_CLOSURE, finite_bias_global_detailed_balance="not_assumed",
        potential_interpretation="effective_network_voltage_used_as_ionic_electrostatic_proxy",
        geometry=dict(cell_A=simulator.lattice.cell_A, normal=geometry.normal, thickness_nm=geometry.thickness_A/10.,
                      oxide_sites=simulator.n, pbc=[True, True, False],
                      bottom_contact_sites=np.flatnonzero(geometry.bottom_mask),
                      top_contact_sites=np.flatnonzero(geometry.top_mask),
                      explicit_metal_atoms=False),
        oxygen_boundaries={name: ("finite_reservoir" if pool.config.enabled else "blocking")
                           for name, pool in simulator.reservoirs.items()},
        electronic_network=simulator.network.summary(),
        run_controls=dict(max_events=max_events, sample_every=sample_every, pulses=trajectory["pulses"]),
        termination=dict(status=trajectory["status"], error=trajectory["error"],
                         target_reached=trajectory["target_reached"], censored=not trajectory["target_reached"],
                         actual_time_s=simulator.time_s, requested_end_time_s=trajectory["requested_end_time_s"],
                         event_count=simulator.event_count),
        initial=initial, final=final, initial_spanning=simulator.initial_spanning,
        initial_vacancy_indices=simulator.initial_vacancy_indices,
        initial_oxygen_atom_ids=initial_oxygen_ids, initial_vacancy_ids=initial_vacancy_ids,
        initial_reservoir_oxygen_atom_ids=initial_pool_ids,
        final_reservoir_oxygen_atom_ids=simulator.reservoir_oxygen_ids,
        event_counts=dict(simulator.counts),
        reservoir_ledger={name: dict(initial=pool.initial_count, final=pool.count,
                           oxygen_out=pool.out_count, oxygen_in=pool.in_count,
                           capacity=pool.config.capacity, reaction_electrons=simulator.electrons_to_electrode[name])
                          for name, pool in simulator.reservoirs.items()},
        conservation=dict(initial_total_oxygen=simulator.initial_total_oxygen,
                          final_total_oxygen=final["total_oxygen"],
                          oxygen_balance_residual=final["mass_residual"],
                          reaction_charge_residual_e=final["reaction_charge_residual_e"],
                          oxygen_atom_identity_check="performed after every event",
                          max_saved_mass_residual=max(abs(r["mass_residual"]) for r in records),
                          max_saved_reaction_charge_residual_e=max(abs(r["reaction_charge_residual_e"]) for r in records)),
        numerical_checks=dict(
            max_saved_relative_kcl_residual=max(r["relative_kcl_residual"] for r in records),
            max_saved_relative_terminal_imbalance=max(r["relative_terminal_current_imbalance"] for r in records),
            max_saved_relative_power_balance_residual=max(r["relative_power_balance_residual"] for r in records),
            transport_residual_checks="also enforced by each electronic solve, including between saved samples"),
        formation=dict(
            criterion="vacancy_graph_spanning AND G>=minimum_conductance AND G/G_initial>=minimum_gain",
            minimum_conductance_S=config.formation_min_conductance_S, minimum_gain=config.formation_min_gain,
            initial_spanning=simulator.initial_spanning, first_spanning_time_s=simulator.first_spanning_time_s,
            first_criterion_time_s=simulator.first_formation_time_s,
            final_criterion_met=final["model_formation_criterion"],
            longest_observed_criterion_duration_s=final["longest_criterion_duration_s"],
            criterion_transitions=simulator.formation_observations,
            interpretation="A declared effective-network criterion, not verified real-material filament formation or retention."),
        diffusion=dict(available=False, reason="Vacancies are born and removed; no fixed-cohort tracer MSD or bulk diffusion inference is made."),
        assumptions=[
            "Fixed oxygen-site lattice, fixed vacancy charge branch, isothermal; no interstitial oxygen or bulk Frenkel creation.",
            "Only enabled contact regions exchange neutral stored oxygen; finite storage free energy includes binomial entropy once when requested.",
            "Formation offset is independent of migration-energy zero; epsilon_i-epsilon_reference_basis supplies relative site chemistry.",
            "Reaction electron counts compensate changes in vacancy charge relative to the initialized state; they are not the continuous terminal electron current or a Poisson charge solution.",
            "Shared saddles satisfy local detailed balance at one frozen site/electrode potential; fields are recomputed after events. Global finite-bias detailed balance is not assumed.",
            "Positive effective resistor graph is independent of ionic/interaction graphs; vacancy-marked conductance is illustrative, not quantum or DFT-calibrated transport.",
            "Actual compliance-limited voltage drives kinetics. Electrical power is diagnostic; no thermal feedback is applied.",
            "Geometric vacancy spanning, electrical conductance threshold, and observed persistence are recorded separately; no experimental SET/RESET/retention claim.",
            "Reservoir atom labels use a deterministic last-in bookkeeping convention; no microscopic electrode diffusion or reservoir tracer dynamics is resolved.",
        ],
        files=["coupled_config.json", "catalog.json", "summary.json", "trace.csv", "events.csv", "oxygen_ledger.csv",
               "states.npz", "electronic_edges.csv", "contact_sites.csv", "final_transport.npz"]+
              (["coupled_overview.png", "viewer.html"] if plot else []),
    )
    if directory is not None:
        config.write(directory/"coupled_config.json")
        catalog.write(directory/"catalog.json")
        for row in events:
            row.update(catalog_sha256=digest, kinetic_field_closure=FIELD_CLOSURE,
                       energy_increment_interpretation="frozen-field local chemical free-energy change plus electrical work",
                       temperature_K=config.temperature_K, seed=config.seed,
                       exchange_provenance=config.exchange.provenance)
            if row["kind"] == "hop":
                datum = catalog.channel(row["channel_id"], config.charge)
                row["channel_provenance"] = datum["provenance"]
                row["catalog_saddle_eV"] = datum["saddle_eV"]
        _write_csv(directory/"trace.csv", records, ("time_s", "events", "pulse_index", "vacancies", "oxide_oxygen",
            "bottom_reservoir_oxygen", "top_reservoir_oxygen", "total_oxygen", "mass_residual", "reaction_charge_residual_e",
            "requested_voltage_V", "applied_voltage_V", "terminal_current_A", "equivalent_conductance_S"))
        _write_csv(directory/"events.csv", events, ("event_index", "time_s", "dt_s", "kind", "source", "target",
            "reservoir", "oxygen_id", "oxygen_from", "oxygen_to", "vacancy_id", "vacancy_identity_action"))
        _write_csv(directory/"oxygen_ledger.csv", [dict(
            event_index=row["event_index"], time_s=row["time_s"], kind=row["kind"], reservoir=row["reservoir"],
            oxygen_id=row["oxygen_id"], oxygen_from=row["oxygen_from"], oxygen_to=row["oxygen_to"],
            vacancies_before=row["vacancies_before"], vacancies_after=row["vacancies_after"],
            reservoir_oxygen_before=row["reservoir_oxygen_before"], reservoir_oxygen_after=row["reservoir_oxygen_after"],
            mass_residual=row["mass_residual"], reaction_charge_residual_e=row["reaction_charge_residual_e"]
        ) for row in events], ("event_index", "time_s", "kind", "oxygen_id", "oxygen_from", "oxygen_to"))
        network = simulator.network
        np.savez_compressed(directory/"states.npz", time_s=np.asarray([r["time_s"] for r in records]),
            events=np.asarray([r["events"] for r in records], dtype=np.int64),
            vacancy=np.asarray(trajectory["states"], dtype=bool),
            oxygen_ids=np.asarray(trajectory["oxygen_ids"], dtype=np.int64),
            vacancy_ids=np.asarray(trajectory["vacancy_ids"], dtype=np.int64),
            potential_V=np.asarray(trajectory["potentials_V"]),
            edge_current_A=np.asarray(trajectory["edge_current_A"]),
            edge_field_V_nm=np.asarray(trajectory["edge_field_V_nm"]),
            oxygen_A=simulator.lattice.oxygen_A, hf_A=simulator.lattice.hf_A,
            cell_A=simulator.lattice.cell_A, pbc=np.asarray(simulator.lattice.pbc),
            depth_A=geometry.depth_A, bottom_contact=geometry.bottom_mask, top_contact=geometry.top_mask,
            electronic_pairs=network.pairs, electronic_displacements_A=network.displacements_A)
        _write_csv(directory/"electronic_edges.csv", [dict(edge_id=index, source=int(pair[0]), target=int(pair[1]),
            displacement_A=network.displacements_A[index], distance_A=float(network.distances_A[index]))
            for index, pair in enumerate(network.pairs)])
        _write_csv(directory/"contact_sites.csv", [dict(site_id=i, x_A=float(p[0]), y_A=float(p[1]), z_A=float(p[2]),
            depth_A=float(geometry.depth_A[i]), region=str(geometry.region[i]),
            bottom_contact=bool(geometry.bottom_mask[i]), top_contact=bool(geometry.top_mask[i]))
            for i, p in enumerate(simulator.lattice.oxygen_A)])
        final_transport = simulator.fields.to_dict(include_arrays=True)
        np.savez_compressed(directory/"final_transport.npz", **{key: value for key, value in final_transport.items()
                           if isinstance(value, (list, tuple, np.ndarray, float, int, bool)) and value is not None})
        _write_json(directory/"summary.json", summary)
        if plot:
            from .coupled_plot import plot_coupled_run
            plot_coupled_run(directory)
    return summary
