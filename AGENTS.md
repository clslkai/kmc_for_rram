# Research development instructions

The current development package is `dftkmc/`. `atomistic/`, `vcm3d/`, their
original tests, and the original C sources preserve historical baselines.

The v3 scope is fixed-lattice, fixed-charge, isothermal, vacancy-conserving KMC
with a prescribed uniform field. State clearly when a change expands this scope.

- Never fabricate completed DFT/NEB calculations or provenance. Unknown required
  parameters stay `null`; synthetic runnable examples are `illustrative` and
  require explicit CLI permission.
- Use the same state energy and common transition state for forward/reverse
  rates. Do not add a NEB activation barrier to the legacy uphill penalty.
- Keep migration and interaction graphs independent. Geometric clusters are
  not evidence of electronic conduction or a switching filament.
- Preserve atom identity, full-cell periodic image displacement, charge, and
  atom count in every exported NEB path. Endpoints require relaxation.
- Fixed charge branches do not implement electron capture/emission. Do not
  compare absolute energies between branches with independent energy zeros.
- Compare independent replicas at a common physical time. Report event-cap
  censoring; no inferential diffusion estimate from an incomplete replica set.
- A single observation time gives a finite-time tracer estimate, not a
  converged long-time diffusion coefficient or an experimental prediction.
- For physical/numerical changes, run `python -m unittest discover -s tests -v`.
  Add meaningful invariant or analytic tests for new physics, not tests that
  simply repeat the implementation. Update `DFT_KMC_V3.md` for changed assumptions.

See `VALIDATION_V3.md` for the baseline, validation evidence and reproducible demos.
