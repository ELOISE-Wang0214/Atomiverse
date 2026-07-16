---
title: CPU breadth 03 — charged variable-cell NPT with Ewald/PME
slug: cpu-breadth-03-charged-npt-kspace
status: approved
execution: cpu_complete
last_progress_update: 2026-07-16
science_required: true
created: 2026-07-15
acceptance: cpu-breadth-03-charged-npt-kspace.acceptance.md
depends_on: cpu-breadth-01-acceptance-foundation
---

# CPU breadth 03 — charged variable-cell NPT with Ewald/PME

## Summary

Compose existing CPU subsystems into the first accepted charged,
variable-cell application: `CoulLong + KSpaceEwald/KSpacePME +
NoseHooverNPT`. The feature is a minimal isotropic/orthorhombic NPT run on a
neutral ionic crystal, not a new electrostatics method.

This probe asks whether Box mutation, stress capability, long-range parameter
pairing, KSpace cache invalidation, FFT workspace lifetime, MPI reduction, and
the MD Plan remain correct when exercised together. Each piece has isolated
tests today; their composition does not.

## Domain basis

Real and reciprocal Ewald terms share one splitting parameter `alpha` [1/Å].
The reference static system is a neutral NaCl lattice. Its Coulomb energy per
ion pair follows the Madelung expression `E = -M k_e q^2 / r0`, with the NaCl
Madelung constant and lattice convention documented in `REFERENCE.md`.
Uniformly scaling all coordinates and the cell by `s` scales the pure Coulomb
energy by `1/s`; this is the cache-invalidation oracle at two volumes.

The dynamic path uses Martyna-Tobias-Klein/Nose-Hoover NPT equations already
implemented by `NoseHooverNPT`. Acceptance checks the extended conserved
quantity drift for a short deterministic isotropic run and validates pressure,
energy, force, and stress against LAMMPS at selected frames. PME references
Essmann et al., JCP 103, 8577 (1995); Ewald/Madelung and MTK references are
listed with exact equations in the regression/validation references.

Units: distance [Å], energy [eV], force [eV/Å], stress/pressure in the engine's
documented conversion, `alpha` [1/Å], temperature [K], time [fs].

## Progress / current evidence (2026-07-16)

**Execution: CPU-complete — 9/9 tasks; 9/9 active acceptance criteria
verified. GPU KSpace remains deferred (ac-010).**

- `LongRangeParticipant` lets ordinary components declare real/reciprocal role
  and `alpha`. CPU Plan validation now rejects missing, duplicate, and
  mismatched halves before launch; actual CoulLong+Ewald and CoulLong+PME Plan
  tests pin the contract.
- `Box` has monotonic mutation revision and cell caches key on `(identity,
  revision)`, preventing both stale in-place NPT data and accidental reuse
  across equal-revision replica Frames.
- PME spline workspace is component-owned and capacity-persistent. Lifecycle
  counters distinguish build, cell refresh, and capacity/grid rebuild; FFTW
  plan ownership is move-only so refresh cannot duplicate a raw plan handle.
- **PME reciprocal virial is implemented** and advertises `Capability::Stress`.
  Unit gates pin Ewald parity, finite-strain diagonal oracle
  (`V_αα = −∂E/∂ε_α`), and NPT+PME Plan acceptance.
- **Public regression** `regressions/md/npt_pme_nacl/` delivers Madelung +
  `1/s` scaling (pure electrostatic Evaluate) and short MTK NPT with soft LJ
  core, cell motion, and conserved-H bound.
- **LAMMPS parity** (`bench_pme_nacl_lammps`): three static variable-cell
  frames vs pinned coul/long+pppm reference (energy, 6-virial, pressure,
  volume, force floor). PASS.
- **1/2-rank MPI parity** (`test_mpi_pme_nacl_parity`): static Madelung +
  tag-gather + short NPT rank-identical H/V on n=1 and n=2. PASS.
- **ac-008/ac-009 closed:** stable-grid lifecycle/memory gate +
  `docs/developer/kspace-long-range.md`.

## Architecture questions

1. Is Box mutation observable to every cell-dependent component without
   letting each method invent an unsafe cache policy?
2. Does `Capability::Stress` mean the component actually contributes the
   complete tensor needed by the active barostat?
3. Can real-space and reciprocal-space halves declare and validate one coupled
   electrostatics contract while remaining ordinary components?
4. Can a changing cell refresh coefficients without turning every MD step
   into an allocation or FFT-plan rebuild?

## Design

### Fail-closed long-range pairing

Add a narrow opt-in long-range electrostatics protocol implemented by
`CoulLong`, `KSpaceEwald`, and `KSpacePME`. It exposes only the information the
engine must validate: real/reciprocal role, splitting `alpha`, and required
cutoff/accuracy metadata. MD validation rejects:

- a real-space long-range component without one reciprocal partner;
- a reciprocal component without the real-space partner;
- duplicate partners;
- unequal `alpha` beyond an explicitly documented floating comparison;
- an NPT Plan whose selected reciprocal component lacks real Stress output.

The protocol is not a new evaluator hierarchy and performs no computation.

### Separate cell refresh from capacity rebuild

KSpace components distinguish three events:

1. **build**: allocate atom/grid scratch and create FFT plans;
2. **cell refresh**: update reciprocal basis, volume factors, influence data,
   and spline moduli when the cell changes but grid shape/capacity does not;
3. **rebuild**: grow atom scratch or change grid dimensions/FFT plans, exposed
   and counted as a cold event.

PME moves per-atom `base/theta/dtheta` scratch to component-owned storage sized
at build/rebuild. A same-shape cell refresh overwrites existing storage and
uses `fill`, not `assign` that reallocates. A grid-dimension change may replan,
but the counter and timing make it visible. Ewald follows the same semantic
event split for reciprocal-vector tables.

### Complete stress, then couple NPT

PME implements the reciprocal virial/stress tensor from the same influence
function used for energy/force, including self/net-charge terms as required by
the selected convention. Finite strain derivatives and Ewald parity verify all
six components. It may advertise `Capability::Stress` only after those checks
pass.

The public MD Plan then composes velocity integration, `NoseHooverNPT`,
NeighborList, `CoulLong`, and one KSpace component normally. NPT receives the
fully reduced total stress before its end-of-step cell update. No component
calls another component directly.

### Scope of the first case

The first regression is neutral and has no bonded exclusions, so it isolates
long-range/cell behavior. Special-bonds reciprocal correction is not silently
declared solved; it becomes active in link 04's molecular ForceField pipeline.
Only orthorhombic isotropic NPT is accepted here. Triclinic/anistropic behavior
remains an explicit CPU follow-up rather than an implicit fallback.

## Files to create or modify

- a narrow protocol under `include/atomiverse/cpu/component/`
- `include/atomiverse/cpu/potentials/coul_long.h`
- CPU Ewald/PME headers and implementations
- CPU MD Plan validation / long-range registration path
- KSpace stress, lifecycle, allocation, and pairing unit tests
- `regressions/md/npt_pme_nacl/` (new, with `REFERENCE.md`)
- `benchmarks/md/npt_pme/{validation,performance}/` (new)
- CPU full acceptance profile/CMake registrations
- `.claude/notes/notes.md`

## Tasks

- [x] Write failing tests for missing/duplicate/mismatched long-range partners
      and for PME's claimed-but-missing Stress contribution.
- [x] Add finite-strain stress tests and static NaCl Ewald/PME energy/force
      parity at the original and uniformly scaled cell.
- [x] Write the public charged NPT regression (Madelung + NPT conserved);
      LAMMPS validation remains open (ac-006).
- [x] Introduce the declarative pairing protocol and enforce it in CPU MD
      validation without welding KSpace into the Driver type hierarchy.
- [x] Implement complete PME reciprocal stress and reduce it consistently
      under MPI (rank-0 reciprocal energy/stress; replicated-grid Allreduce
      of the charge mesh matches Ewald).
- [x] Move PME atom scratch to persistent capacity and split cell refresh from
      grid/capacity rebuild; instrument refresh/replan/allocation counts.
- [x] Ensure NPT phase order consumes current total stress and mutates Box only
      after the force/stress pipeline for that step (public NPT regression
      exercises the existing MTK schedule with complete reciprocal stress).
- [x] Run two-rank unit parity (test_mpi_pme_nacl_parity n1/n2); controlled
      lifecycle/performance gate (ac-008) remains open.
- [x] Record the cell-cache and long-range pairing contracts in durable docs
      and `.claude/notes/notes.md`.

## Testing strategy

- **Unit**: protocol pairing, all six stress components by finite strain,
  cell-change invalidation, no-refresh unchanged-cell path, explicit regrid,
  and 1/2-rank reductions.
- **Regression**: one public NaCl charged variable-cell case, Madelung/static
  scaling checks plus short NPT conserved-quantity bound, ≤60 seconds.
- **Validation**: LAMMPS energy/force/pressure/stress at fixed selected frames;
  Ewald is an independent in-tree numerical cross-check for PME, not the sole
  scientific oracle.
- **Performance**: steady cell refresh versus rare grid replan reported
  separately; stable-grid NPT steps allocate zero after warmup.

## Deferred — GPU

- CUDA CoulLong/Ewald/PME pairing validation and complete stress parity.
- CUDA NPT cell-dependent cache refresh and device-side lifecycle counters.
- Single- versus multi-GPU charged NPT and performance parity.

No CUDA/Torch production source or test is changed by this link.

## Out of scope

- Triclinic/anistropic NPT and variable-cell minimization.
- Distributed FFT/PME decomposition; the accepted MPI fixture uses the current
  replicated-grid algorithm.
- Molecular 1-2/1-3/1-4 reciprocal exclusion correction (link 04).
- ForceField parsing/typing, Python NPT surface, or new electrostatic styles.
