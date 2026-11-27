# Model definition, assumptions and validation boundary

## 1. What is represented

The original `kmc-ecm.c` is left unchanged. It describes metal oxidation/reduction and diffusion on a 2D lattice; it is **not** an oxygen-vacancy VCM model. The new Python code is a separate model, not a verified port of that code.

Each Cartesian voxel has two binary, distinct coarse-grained sites: one oxygen lattice site (vacancy flag `v`) and one interstitial site (`o`). Both may be occupied in the same voxel; they are not the same atomic site. Positions at voxel centers are an approximation for the unresolved sublattices. A same-cell Frenkel pair has zero resolved dipole. Pair creation across a cell or a nearest-neighbor cell and its exact reverse are both enumerated. Defects cannot leave any boundary. The model represents a **closed oxygen inventory**, not an electrode oxygen reservoir.

The bookkeeping invariant is `sum(v) == sum(o)` for the supplied neutral initial conditions. Pair creation removes an oxygen from its lattice site and places it on an interstitial site; recombination reverses that. Vacancy motion is oxygen motion in the opposite direction. Initial interstitials are an effective representation of oxygen displaced by pre-existing vacancies. Initial states with missing oxygen in an external reservoir are not implemented.

All six nearest-neighbor directions are present in the volume. There is no periodic wrapping; lateral electrical/thermal boundaries are insulating, and material/defect boundaries are reflecting. Electrode faces are half a voxel beyond first/last centers in z. Dimensions are `nx*a × ny*a × nz*a`.

## 2. Configurable phenomenological energetics

For *frozen* potential, an effective configurational energy in eV is

`U = (F/2)*(Nv + Ni) - J*Bvv + q*sum_i[(v_i-o_i)*phi_i]`,

where `F` is pair formation energy, `J` is the attractive nearest-neighbor vacancy bond energy, `Bvv` counts each vacancy-vacancy bond once, `q` is an effective charge in units of the elementary charge, and phi is in volts. Thus `q*delta_phi` is numerically in eV; do not multiply again by the elementary charge when used with eV barriers.

Event energy changes:

- pair generation at oxygen site i and interstitial site j: `dU = F - J*n_v(i) + q*(phi_i-phi_j)`;
- reverse recombination: negative of that value, using the existing vacancy's neighbors;
- vacancy hop i→j: `dU = q*(phi_j-phi_i) - J*(n_v(j)-1-n_v(i))`;
- interstitial hop i→j: `dU = -q*(phi_j-phi_i)`.

The `-1` in vacancy-hop neighbor bookkeeping excludes the vacancy being moved from i. Charge-state transitions, defect-defect Coulomb repulsion, strain, entropic degeneracies, structural phase transformations and oxygen molecules are not resolved. A fixed charge of +2/-2 is an effective approximation, especially questionable in a dense conducting filament.

For each reversible mechanism,

`rate = nu * exp[-(E0 + max(0,dU))/(kB*T_local)]`,

with `T_local=(T_i+T_j)/2`. E0 is a nonnegative intrinsic barrier. This Metropolis transition-state choice ensures a nonnegative barrier, rate at most nu, and forward/reverse ratio `exp(-dU/kBT)` for frozen fields and equal local temperature. It is one defensible coarse-grained choice, **not a DFT-derived saddle point**. At nonuniform temperature it is an effective bond-temperature approximation.

The electrical potential is recomputed after each state change. The above energy differences use the pre-event field and omit the full electrode-work/field-relaxation free-energy change. Therefore the fixed-field detailed-balance tests must **not** be read as proof of global equilibrium detailed balance of a self-consistent electrostatic many-body model. With applied voltage the system is driven and generally nonequilibrium.

## 3. Electrical network, not Poisson electrostatics

We solve `div(sigma*grad(phi))=0`, discretized as a 3D resistor network. A vacancy-containing cell has sigma_on; all others have sigma_off. An internal edge has harmonic-mean conductivity times `a` in siemens; an electrode contact has `2*sigma*a` (half-cell distance). Every current path is represented, unlike independent column resistance estimates.

This is a **quasi-neutral current-continuity approximation**, not a solution of `div(epsilon*grad(phi))=-rho`. The model does not calculate long-range defect space charge or screening. Nor does it implement tunneling or hopping electronic conduction. The harmonic edge law makes a vacant/nonvacant interface bottleneck, which is a modeling choice, not measured HfO2 transport.

The unit-bias linear solve gives terminal conductance G. Ideal instantaneous current compliance uses

`V_device = sign(V_applied)*min(abs(V_applied), I_comp/G)`.

The limited voltage is used for every potential, Joule source and event rate in that state. This is not a transistor load line or a finite-bandwidth circuit. A finite nonzero sigma_off keeps the network connected and numerically well-posed. A top-to-bottom connected vacancy component is recorded as a **geometric** filament criterion; electrical conductance remains finite with or without it.

## 4. Heat

For each internal resistor, dissipated power is `g*(phi_i-phi_j)^2`; half is assigned to each adjacent cell. For each electrode contact, half the contact dissipation is assigned to the cell and half to the ideal electrode. Total electrical dissipation, including full contact power, is checked against `I*V_device`.

A separate linear heat-conduction network solves `K*(T-T_ambient)=P_cells`, with constant thermal conductivity, isothermal top/bottom electrodes and insulated lateral faces. Heat flux into electrodes equals deposited **cell** power. This corresponds to quasistatic relaxation between events, not a transient thermal equation. It has no temperature-dependent conductivity, interface thermal resistance or reaction enthalpy. A configurable maximum temperature aborts the calculation rather than silently clamping temperature and altering kinetics.

## 5. KMC time and pulses

All legal events are explicitly enumerated with rates r_i. Select an event with probability r_i/R, R=sum(r_i), and advance time by `-ln(u)/R`. Random uniforms for event selection and waiting time are independent. This direct stochastic simulation algorithm is exact for the stated state-dependent rate model within a constant-voltage interval, **not** for all underlying materials physics.

If the sampled event lies beyond the next voltage boundary, advance only to that boundary, do not fire the event, and rebuild fields/rates for the new voltage. Memorylessness justifies a new waiting-time draw. A cap on event count is a runtime safety limit, not a physical stopping rule; output states explicitly whether the prescribed waveform finished. Diagnostics and stored potentials/temperatures are recomputed for the **post-event** state.

## 6. Default illustrative values and units

| Quantity | Value | Unit / meaning |
|---|---:|---|
| voxel side | 0.5 | nm; coarse graining, not an explicit HfO2 lattice |
| ambient | 500 | K; elevated demonstration temperature |
| attempt frequency | 1e12 | s^-1 |
| effective vacancy/interstitial charges | +2 / -2 | elementary charge units |
| pair energy F | 0.8 | eV; effective placeholder, not isolated HfO2 formation energy |
| pair reversible barrier E0 | 0.45 | eV |
| vacancy / interstitial hop E0 | 0.65 / 0.70 | eV |
| vacancy bond J | 0.08 | eV per nearest-neighbor pair |
| sigma_off / sigma_on | 1e-3 / 1e3 | S/m |
| thermal conductivity | 1 | W/(m K) |
| compliance | 10 | microampere |
| temperature validity limit | 1500 | K; stop, not clip |

These numbers support reproducible numerical examples only. In particular the effective formation energy is deliberately not presented as a real HfO2 defect formation energy. No switching voltage, time, on/off ratio or endurance derived from them is a validated prediction. At 300 K the dynamics may be extremely slow; increasing temperature just to get activity changes the physical problem.

## 7. Verification and limitations

Tests cover homogeneous analytic conductance and linear potential, signed polarity, six-neighbor topology, KCL, electrical and thermal power balance, current compliance, forward/reverse energy/rate consistency, exact pulse cutoff, oxygen conservation, reproducibility and geometric connectivity. Numerical consistency does not establish material realism.

The supplied random-state trajectory has 16 saved samples with a geometrically spanning vacancy cluster, followed by disconnection; its final state is nonspanning. This is transient behavior within the illustrative coarse-grained conductance model, not validation of a stable physical filament or real SET. The seeded example starts with an explicitly imposed conducting chain and may break because random displaced oxygen is initially nearby; this is not independent evidence for voltage-controlled RESET. Demonstrating robust reversible forming/SET/RESET requires further physical ingredients and calibration, not just tuning these demo values until attractive plots appear.

This implementation intentionally enumerates all events and solves sparse systems after every event for clarity. It is suitable for small reference grids. Production-scale work requires validated local event updates, sparse-factor reuse where justified, event trees, and a convergence study; never simply multiply all rates without rescaling time and understanding its consequences.

## 8. Numerical-method sources

- D. T. Gillespie, *Exact stochastic simulation of coupled chemical reactions*, J. Phys. Chem. 81, 2340–2361 (1977), https://doi.org/10.1021/j100540a008 — direct stochastic simulation method.
- See the project assessment for material-specific VCM literature. The constants here are intentionally not attributed to a particular fitted paper or specific fabricated stack.

## 9. Frozen-state read diagnostic

`conductance_S` and `resistance_ohm` measure the same instantaneous resistor network independent of write-voltage polarity. `frozen_read_0p1V_current_A = G*0.1 V` is a hypothetical ohmic readout without defect evolution or a read-compliance circuit; it is not a separately simulated read pulse. Since conductivity is temperature-independent here, it is also independent of read heating. This simplification must be replaced if transport acquires temperature or field dependence.
