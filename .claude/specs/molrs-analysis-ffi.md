# Trajectory Analysis via molrs FFI

status: draft
date: 2026-07-01
depends on: [manual-system](manual-system.md), [manual-discovery](manual-discovery.md)
policy: [[feedback-infrastructure-in-molrs]] — infrastructure lives in molrs, not C++.

## Problem

The manual's Analysis & Observables chapter (Ch 13) and the discovery report
flagged MSD / self-diffusion / VACF / Green-Kubo / structure-factor / ADF as
"missing." **They are not missing — they are unexposed.** molrs already ships a
full freud-class analysis library in `molrs/src/compute/` behind one stateless
`Compute` trait (`&[&Frame] + args → result`) plus a `Fit` trait
(`observable → derived quantity`):

- `compute::msd` → `MSD` / `MSDTimeSeries` / `MsdMode`
- `compute::rdf` → `RDF` / `RDFResult`
- `compute::van_hove`, `compute::distribution` (RDF/ADF), `compute::diffraction`
  (static structure factor, Debye + direct)
- `compute::fit` → `EinsteinDiffusion`, `GreenKuboDiffusion`,
  `GreenKuboConductivity`, `VACF`, `PowerSpectrum`, `LinearFit`, `Plateau`
- `compute::onsager`, `compute::hbond`, `compute::order` (Steinhardt, hexatic, …)

The C++ FFI bridge (`molrs-cxxapi/src/bridge.rs`) currently exposes only Frame /
MolRec / I-O. The gap is **the analysis kernels are not bridged to C++.**

## Design — the engine measures, molrs analyzes

The substrate already exists: `cpu::MolRec` (PIMPL, `cpu/io/molrec.h`) accumulates
a trajectory *inside molrs* via `commit_frame`. Analysis runs on that
molrs-resident trajectory. No new C++ analysis code, no C++ data structures.

```
C++ MD/AIMD loop → Recorder stages frames → MolRec (trajectory in molrs)
                                                     │
                                   molrs Compute (MSD/RDF/VACF) + Fit (diffusion)
                                                     │
                        FFI returns result arrays → thin C++ Observer exposes them
```

### FFI surface (add to `molrs-cxxapi/src/bridge.rs`)

Operate on the accumulated `MolRec` (or a frame handle sequence). Return plain
POD arrays (CXX `Vec<f64>` / `Vec<i64>`), never Rust analysis objects:

```
// Whole-trajectory observables (accumulate first, analyze at finalize):
fn analyze_rdf(rec: &AtvMolRec, r_max: f64, n_bins: i64, type_i: i32, type_j: i32) -> RdfOut; // r[], g[]
fn analyze_msd(rec: &AtvMolRec, mode: i32) -> MsdOut;         // tau[], msd[]  (mode = MsdMode)
fn fit_einstein_diffusion(tau: &[f64], msd: &[f64], dims: i32, t_lo: f64, t_hi: f64) -> DiffusionOut; // D, intercept, r2
fn analyze_vacf(rec: &AtvMolRec, max_lag: i64) -> VacfOut;    // tau[], vacf[]
fn analyze_structure_factor(rec: &AtvMolRec, q_max: f64, n_q: i64) -> SkOut; // q[], S[]
```

Each returns a small POD struct of `Vec<f64>`. Errors surface as a `Result`
mapped to a C++ exception (NOT `.expect()` — see the extxyz abort bug in
manual-discovery §12).

### C++ side — thin wrappers only

- A `role()==Observer` component `cpu::observer::TrajectoryAnalysis` (or discrete
  `MSD` / `VACF` observers) that: subscribes to snapshots, stages each into a
  `MolRec`, and at `finalize()` calls the molrs FFI and writes results to disk /
  exposes them via a small POD result struct.
- These are ≤~60-line wrappers with zero analysis math — all physics is in molrs.
- The existing hand-rolled C++ `cpu::observer::RDF` is **tech debt** under this
  policy: eventually route it through `analyze_rdf`; do not add more C++ analysis.

## Tasks (smallest-correct-first)

1. ✅ **MSD** (unlocks Ch 13, done 2026-07-01): `analyze_msd(&AtvMolRec) -> Vec<f64>`
   (molrs `compute::MSD` Direct = LAMMPS `compute msd`); C++ `cpu::MolRec::compute_msd()`
   + async `cpu::observer::MSD` staging snapshots via the existing snapshot ring,
   computing at `finalize()`; `test_cpu_observer_msd` passes vs analytic MSD.
   Recipe for adding bridge fns: [[project-molrs-ffi-recipe]] (edit build.rs schema,
   enable `compute` feature, rebuild twice). **Next:** `fit_einstein_diffusion`
   (molrs `fit::EinsteinDiffusion`) for D; capture image flags for PBC unwrapping
   (snapshot lacks them — MVP validated on a non-wrapping trajectory).
2. ✅ **RDF via molrs** (closed 2026-07-10): landed as the *streaming* route
   rather than batch-`analyze_rdf` — molrs `RdfAccumulator` behind
   `cpu::compute::RDFAccumulator` (trajectory_metrics.h). CPU observer
   migrated 2026-07-02; the CUDA tail (`cuda::RDF` + deletion of the last
   hand-rolled histogram `rdf_histogram.{h,cpp}`) is closed by
   **arch-refactor-07-molrs**. No g(r) analysis math remains in C++ on either
   backend; the batch `analyze_rdf` bridge fn is no longer needed for this line.
3. ✅ **VACF** (done 2026-07-01): `analyze_vacf(velocities, n_frames, n_dof, dt,
   resolution)` runs molrs `compute::VACF` (FFT-ACF, DOF-averaged) over a flat
   velocity buffer (VACF takes velocities in Args, ignores frames). C++
   `cpu::compute_vacf` + async `cpu::observer::VACF` (captures velocities from the
   snapshot ring). Test `test_cpu_observer_vacf`. **Next for this line:**
   Green-Kubo D (`fit::GreenKuboDiffusion` = running-integral of the VACF);
   cross-check against the Einstein D.
4. **Static structure factor** (Ch 13 stretch): bridge `analyze_structure_factor`.
5. **Extra trajectory formats** (DCD/TRR/PDB already in molrs): bridge the readers
   the discovery listed as "missing in C++."

### Cross-repo note — lammps_data.cpp replacement preconditions (out of scope here and in arch-refactor-07)

`src/cpu/io/lammps_data.cpp` (~259-line hand-written parser) is a **deliberate
divergence**: molrs's LAMMPS reader atom-style auto-detection mis-routes MSU
`atom_style=molecular` files onto the charge branch (mol_id read as type).
Replacing it requires, in order: (1) molrs-side detection-heuristic fix,
(2) a LAMMPS read FFI added to the bridge, (3) MSU PE benchmark regression in
Atomiverse. Until all three land, the C++ parser stays. Structure factor
(Task 4) and DCD/TRR/PDB (Task 5) also remain owned by this spec —
arch-refactor-07 closed only Task 2's CUDA tail.

## Validation

- Rust side: molrs already has unit tests per kernel — reuse as the source of
  truth; the C++ wrapper test asserts FFI round-trip equals the Rust result.
- Physical consistency: Einstein-D ≈ Green-Kubo-D within documented tolerance.
- Determinism: fixed trajectory in, identical arrays out (stateless `Compute`).

## Constraints

- Cross-repo: modifies molrs (Rust, currently `dev` 4 commits ahead of origin) +
  Atomiverse C++; requires a corrosion rebuild. Build locally is fine.
- Errors are catchable C++ exceptions, never process aborts.
- No analysis math in C++ — that would violate [[feedback-infrastructure-in-molrs]].
