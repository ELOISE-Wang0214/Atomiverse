# Atomiverse Manual — Discovery Report

status: reference (discovery phase — no rewrites)
date: 2026-07-01
purpose: Ground-truth map of the codebase before building the textbook-level
C++ manual. Answers the 13 discovery questions. Everything here is grep-verified
against the `dev` branch; file:line citations are evidence, not prescription.

> **Prime directive for the manual:** Atomiverse is *one engine with multiple
> methods and backends*, not a bundle of packages. The architecture below is
> genuinely unified (one `Driver`/`Component`/`Frame` model spanning classical
> MD, minimization, vibrations, Gaussian SCF/DFT, and plane-wave/GPW DFT). The
> manual's job is to make that unity visible.

---

## 1. Current architecture (Q1)

Two orthogonal axes; a concrete driver is one cell of the grid:

- **Backend axis** (`ATV::cpu` / `ATV::cuda` / `ATV::torch`) — selected at
  *compile time* via `DeviceTag` empty tags (`core/device.h`) + `backend_traits`
  specialization (`core/backend_traits.h`). Not runtime polymorphism; separate
  namespaces.
- **Category axis** (termination semantics) — `EvaluateDriver` (one pass),
  `DynamicsDriver` (step budget), `OptimizeDriver` (early-stop on converged),
  and an **undocumented 4th** `ResponseDriver` (finite-difference passes, used by
  vibrations). Bases in `cpu/driver/`.

Primitive: `template<PhaseEnum, DeviceTag> class Driver` (`facade/driver.h:92`).
Backend layer `cpu::CPUDriver<PhaseEnum>` adds `build()`→`install_defaults()`+
`validate_components()`, the step loop, and LAMMPS-shaped setters.

Loop skeleton: `setup_pass → for step: reset outer_status; step_once; consult
outer_status → teardown_pass`. Phases per category are fixed enums in
`facade/stage.h` (md=10 phases, minimizer, evaluate, vib). **No `scf::Phase`** —
SCF is an inner solver, not an outer driver.

## 2. Central abstractions (Q2)

| concept | real name (code) | where |
|---|---|---|
| calculation / loop owner | `Driver<PhaseEnum,DeviceTag>` + category bases | `facade/driver.h`, `cpu/driver/` |
| single-step operation | `cpu::ComponentBase` (`role()`, `capabilities()`, `reads()`, build/launch/destroy) | `cpu/component/component.h:50` |
| structural slot | `ComponentRole` enum (Evaluator, Integrator, Optimizer, Thermostat, Barostat, ConvergenceCheck, Recorder, Observer, None) | `facade/component_role.h:38` |
| written frame fields | `Capability` bitset (Energy, Force, Stress, Charge, Dipole, DensityMatrix, Orbitals, **Hessian**) | `facade/component_role.h:56` |
| device-resident state | backend `Frame` (cpu xtensor / cuda device ptrs / torch `at::Tensor`) | `cpu/frame/frame.h:47` |
| exchange state | `molrs::FrameRef` (Rust CXX-FFI), lowered once/run via `materialize`/`write_back` | `cpu/frame/materialize.h` |
| **result** | `SystemScalars` (energy, kinetic, conserved, stress, per-term energies) | `core/system_scalars.h:17` |
| particles | `cpu::Atom` SoA (x/y/z, LAMMPS per-atom layout), `frame.atoms` | `cpu/frame/atom.h` |
| cell / boundary | `cpu::Box` (min_image, triclinic, SIMD batch) | `cpu/frame/box.h:47` |
| observable | `role()==Observer` components | `cpu/observer/` |
| trajectory / I/O | `role()==Recorder` components | `cpu/recorder/` |
| iteration signal | `Context::outer_status {converged, stop_requested, reason}` | `facade/iteration_status.h` |

**Rules that hold in code:** Driver owns the loop; algorithms are Components
(LBFGS/CG are Optimizer components, VelocityVerlet/NoseHoover are Integrator/
Thermostat components), never driver subclasses; components tagged by
role+capability; no allocation in `launch()`.

## 3. How each method is represented (Q3)

All methods plug into the *same* Driver/Component/Frame model. A force method is
a `role()==Evaluator` component with a `Force` capability; you swap the physics
by swapping the evaluator, keeping the driver:

- **Classical MD** — `cpu::md::MD` (`DynamicsDriver`) + `PairLJCut`/`PairMorse`/
  `CoulLong`/bonded evaluators + `NeighborList` + integrator/thermostat/barostat.
- **Gaussian SCF/DFT** — `cpu::Evaluate` (`EvaluateDriver`) + `cpu::scf::
  SCFEvaluator` (wraps an inner `SCFSolver`). Same evaluator drops into
  `cpu::md::MD` for AIMD and `cpu::minimizer::Minimize` for geometry opt.
- **Plane-wave DFT** — `cpu::pw::PWEvaluator` + `PWSolver`.
- **GPW DFT** — `cpu::gpw::GPWEvaluator` + `GPWSolver`.
- **Minimization** — `cpu::minimizer::Minimize` (`OptimizeDriver`) + LBFGS/CG
  Optimizer component + any Energy/Force evaluator.
- **Vibrations** — `cpu::vib::Hessian` (`ResponseDriver`) drives 6N FD passes of
  any Force evaluator.

This "swap the evaluator, keep the driver" pattern is the manual's spine.

## 4. Backend abstraction (Q4) + parity

Compile-time; separate namespaces; user picks by constructing `cpu::md::MD` vs
`cuda::md::MD`. **CPU is the complete reference backend.** Coverage:

| domain | cpu | cuda | torch |
|---|---|---|---|
| driver taxonomy | ✅ 4 categories | ⚠️ present, doc says "planned" | ❌ none |
| component base | ✅ | ✅ | ❌ (Context only) |
| md | ✅ | ✅ | ❌ |
| minimizer | ✅ (LBFGS/CG) | ⚠️ partial | ❌ |
| scf (Gaussian) | ✅ (33 files) | ⚠️ thin (4) | ❌ |
| pw | ✅ | ⚠️ stub (1) | ❌ |
| gpw | ✅ | ❌ | ❌ |
| kspace (PME/Ewald) | ✅ | ✅ | ❌ |
| observer / recorder | ✅ | ⚠️ partial | ❌ |
| potentials | ✅ (13) | ✅ (12) | ⚠️ 2 (lj, harmonic_bond) |
| vib / xc / io / linalg | ✅ | ❌ | ❌ |

**Critical for manual honesty:** **Torch is NOT a third backend** — no
`device_tag`, no `backend_traits`, no `Driver`/`ComponentBase`. It is a
*differentiable-potential kernel library* (frame + neighbor + 2 potentials as
free functions). The manual must present cpu (reference) + cuda (partial mirror);
torch is an ML-potential adjunct, documented separately.

## 5. Method-status matrix (Q3 + Task #4)

| method / capability | status | evidence |
|---|---|---|
| Classical MD (NVE/NVT/NPT), LJ/Morse/Coulomb-long/Ewald/PME, bonded, cell-list NL, MPI DD/halo | **supported** | `cpu/md/`, `cpu/potentials/`, `cpu/comm/` (real MPI) |
| Minimization LBFGS / CG | **supported** | `cpu/minimizer/` |
| FIRE minimizer | absent (planned) | `minimize.h:6` comment only |
| EAM / tabulated potentials | absent | — |
| RHF / UHF | **supported** (PySCF-validated) | `rhf_solve.h`, `uhf_solve.h`, `test_cpu_scf_basis_vs_pyscf` |
| RKS / UKS (DFT), production default RKS-PBE | **supported** | `scf_solve.h:272`, `test_cpu_scf_ks_*` |
| Hybrid DFT (B3LYP/PBE0) | **supported** | `scf_solve.h:299,315`, `test_hybrid_gradient` |
| Native functionals LSDA/PBE (=libxc 1e-12) | **supported** | `xc_functional.h`, `test_xc_vs_libxc` |
| Analytic gradients (HF+KS+hybrid) | **supported** | `scf/gradient.h` |
| Vibrations (FD Hessian) + thermochem | **supported** (PySCF-validated) | `cpu/vib/`, `test_cpu_vib_*_vs_pyscf` |
| Basis: STO-3G…def2-SVP | **supported**; def2-TZVP (f) | **partial** | `scf/basis_data/` |
| MP2 / CI / CC / CCSD | **absent** | no code |
| Excited-state / TDDFT / response | **absent** | (`ResponseDriver`=FD Hessian, not electronic response) |
| Plane-wave DFT (LDA, Γ-validated), UPF pseudopotentials, Davidson SCF, k-points | **experimental** | `cpu/pw/`, `davidson.cpp`, `upf_reader.cpp` |
| PW no-UPF `V_loc` fallback | **toy/placeholder** | `pw_solver.h:55` |
| GPW DFT (closed-shell, Γ; PBE, GTH, MT-Poisson, OT, DIIS, D3, forces) | **experimental** | `cpu/gpw/` (29 tests) |
| PME / Ewald electrostatics (E/F/stress, CPU+CUDA) | **supported** | `cpu/kspace/`, `cuda/kspace/` |
| Classical analysis: RDF, kinetic VirialStress (C++) | **supported** | `cpu/observer/` |
| MSD / diffusion / VACF / Green-Kubo / S(k) / ADF / order params | **in molrs, unexposed to C++** | `molrs/src/compute/`; needs FFI bridge (molrs-analysis-ffi) |
| Trajectory I/O: XYZ / extxyz / Zarr (bridged) | **supported**; DCD/TRR/PDB in molrs (unexposed); HDF5/LAMMPS-dump absent | `cpu/recorder/` + molrs |

## 6. Examples inventory (Q5, Q6, Q7)

13 example sources; **all use current APIs — zero stale calls found**. Good
foundation. Gates in `examples/CMakeLists.txt`.

| # | file | purpose | gate | state |
|---|---|---|---|---|
| 01 | `01_scf.cpp` | bare SCF loop (pure xtensor) | unconditional | current |
| 02 | `02_rhf.cpp`/`.cu` | RHF single-point H₂ (CPU+CUDA twin) | RHF (+CUDA) | current |
| 03 | `03_uhf.cpp` | UHF + spin + Mulliken/XYZ | RHF | current |
| 04 | `04_uks.cpp` | UKS-LSDA | RHF | current |
| 05 | `05_geom_opt.cpp` | aspirin LBFGS + RHF forces | RHF | current |
| 06 | `06_aimd.cpp` | H₂ BOMD (NVT) | RHF | current (header says 200 steps, code 100 — cosmetic) |
| 07 | `07_md_lj.cpp`/`.cu` | NVE LJ argon (categorical setters) | unconditional (CPU) | current — canonical CPU lesson |
| 08 | `08_pw.cpp` | plane-wave DFT (Si) | PW | current |
| 09 | `09_gpw.cpp` | GPW DFT (H₂O) | GPW+PW | current, **invisible in all indexes**, hardcoded data path |
| 10 | `10_gpw_znoedl.cpp` | GPW 1214-atom ZnO-EDL | GPW+PW | current, **invisible in all indexes**, hardcoded paths |

**Stale/misleading (not the examples themselves, but their discoverability):**
examples 09/10 (the entire GPW method) are missing from `examples/README.md`,
`mainpage.dox`, and top-level `README.md` (all stop at 08). Examples 09/10 also
hardcode developer-specific data paths (same class of defect as the tests, §12).

## 7. Documentation inventory + consistency (Q10, Q11)

Docs already have a real tutorial arc (`docs/tutorial/00`–`08`), a user-guide,
a developer section, and a working Doxygen setup (`docs`/`docs-open` targets,
`ATV_DOCS_ONLY`). **Tutorials 00–08 are accurate** (every code snippet verified).
Rot is concentrated in three pockets:

- **`docs/api/index.md` + `tutorial/07`** disagree on I/O namespaces: code is
  `ATV::RDF` (not `cpu::RDF` nor `cpu::observer::RDF`) and `cpu::XYZ` (not
  `cpu::recorder::XYZ`); `cpu::Trajectory` (api/index.md:75) doesn't exist.
  `tutorial/07:290` uses `sim.observer(...)` — real method is `set_observer`.
- **`docs/developer/*`** carry pre-refactor drift: `architecture.md:234-242`
  still describes a `scf::SCF` driver + `scf::Phase` enum (both removed — SCF is
  `SCFSolver` now); `cpu-internals.md:50-56` omits `driver/` from header paths,
  puts `box.h`/LJ at old locations, cites 3 nonexistent `.claude/specs/*`.
- **`docs/user-guide/libtorch-backend.md`** has a non-compiling example
  (`atorch::lj(...)` missing its `NeighborList` arg; phantom `TorchConfig`) and
  dead test paths (`tests/test_torch_lj.cu` → real `tests/torch/test_lj.cu`).
- **`getting-started.md:119`** mislabels `07_md_lj` as GPU/`ATV_BUILD_CUDA` — it's
  the unconditional CPU target.

Doxygen/API-gen is present and current.

## 8. Benchmarks, tests, CI, validation data (Q8, Q9)

- **Benchmarks** (`benchmarks/`, gated `ATV_BUILD_BENCHMARKS`): `md/lj`, `md/pe`,
  `qc/rhf`, `pw/si`, each `{validation,performance}`. Reference libs: LAMMPS,
  PySCF, CP2K/ABINIT. **CI hygiene is correct** — timing/parity comparisons are
  out of the ctest gate.
- **Tests** (~209 ctest entries): SCF/GPW/MD/PW/MPI deep; parity refs embedded as
  headers (so `*_vs_pyscf` run in CI without the external tool). Local ctest on
  this machine: **117/119 pass** (2 explained — see §12).
- **CI**: `ci.yml` = CPU-only, matrix {float/RHF-off, double/RHF-on}; deliberately
  excludes CUDA/PW/GPW/KSPACE/TORCH. `benchmark.yml` = self-hosted GPU, perf only.

**Benchmark gaps (drive the benchmark-design spec):**
1. **Precision matrix was a no-op** — `ATV_REAL` unwired (fixed 2026-07-01 →
   tiered `core/precision.h`; CI now drives `ATV_REAL_STORAGE`).
2. **No cpu↔cuda backend-consistency job** — the "bit-identical" claim is
   manual-only.
3. **No automated perf-regression gate** (perf CSVs are committed history).
4. **Reference-data rot**: `md/pe` MSU data + all its `tools/` scripts are in
   SOURCES.md but **not committed**; qc generator writes to drifted paths.
5. **No `benchmarks/` case** for gpw / kspace / torch / minimizer / vib.
6. `benchmark.yml` runs a **non-existent `bench_scf`** (`|| true` hides it);
   `test_xc_vs_libxc.cpp` exists but **isn't registered**.

## 9. Terminology inconsistencies (feeds Task #4)

1. RDF: code `ATV::RDF`; docs say both `cpu::RDF` and `cpu::observer::RDF`.
2. XYZ: code `cpu::XYZ`; doc `cpu::recorder::XYZ`.
3. Driver headers live in `cpu/driver/`; several docs (and CLAUDE.md:198,202)
   use flat `cpu/…` paths; `box.h`/`box_simd.h`/`frame.h` are in `cpu/frame/`.
4. SCF is an inner solver (`SCFSolver`), not a `scf::SCF` driver — `architecture.md`
   is the laggard.
5. `Capability` (singular), not `Capabilities`.
6. CUDA minimizer taxonomy: `mainpage.dox:96` + CLAUDE.md say "SD/CG only" but
   `cuda::minimizer::Minimize` exists.
7. 4th driver category `ResponseDriver` + `Capability::Hessian` ship but are
   absent from `component_role.h`'s "3 categories" doc and CLAUDE.md.
8. `LJ`→`PairLJCut` migration is complete in docs/examples (only the cached
   CLAUDE.md snapshot shows the old name; on-disk CLAUDE.md is correct).

## 10. Missing tutorial-critical functionality (Q12)

- **Analysis/observables** — MSD, diffusion, VACF, Green-Kubo, S(k), ADF,
  structure factor, time-correlation are **already implemented in molrs**
  (`molrs/src/compute/`, a freud-class library behind the `Compute`/`Fit`
  traits) but **not bridged to C++** (the FFI exposes only Frame + I/O). Per
  [[feedback-infrastructure-in-molrs]] the fix is **FFI exposure + thin C++
  observer wrappers, not new C++ analysis** — see molrs-analysis-ffi.md. (The
  hand-rolled C++ `cpu::observer::RDF` duplicates molrs and is tech debt.)
- **GPW has no tutorial** (examples 09/10 exist, no narrative).
- **Vibrations has no tutorial or example** (code + tests exist).
- **NPT/barostats/advanced thermostats** only mentioned in passing.
- **Cross-method workflow chapter** (Ch. 8) — the "swap the evaluator" story is
  stated but never run end-to-end (e.g. classical-relax → AIMD → observable).
- Correlation methods (MP2/CC) absent in code — document as *planned/out-of-scope*,
  not as a gap.

## 11. What to fix BEFORE writing new tutorials (Q13)

Ordered, minimal, and mostly mechanical (do these in the cleanup phase, Task #5):
1. Settle & fix the **RDF/XYZ/observer namespace** (`ATV::RDF`, `cpu::XYZ`,
   `set_observer`) in `tutorial/07` + `api/index.md` — a new observables chapter
   will inherit whatever is written.
2. Remove `cpu::Trajectory` from `api/index.md`; point `cpu::XYZ` at
   `recorder/xyz.h`.
3. Repair `developer/cpu-internals.md` + `extending-cpu-subsystems.md` header
   paths (add `driver/`, `frame/`); fix `Minimize` header.
4. Rewrite the SCF section of `developer/architecture.md` to the `SCFSolver`
   inner-solver model; drop `scf::SCF`/`scf::Phase`.
5. Fix `getting-started.md:119` (`07_md_lj` is CPU/unconditional).
6. Surface **GPW (examples 09/10)** in every index + capability list.
7. Fix or quarantine `user-guide/libtorch-backend.md` (non-compiling example +
   dead paths); reframe torch as a differentiable-potential library, not a backend.
8. Correct small slips: `Capabilities`→`Capability`; `scale_diagonal` arity;
   nonexistent `.claude/specs/*` citations; `mainpage.dox:96` CUDA minimizer.
9. Document the 4th category (`ResponseDriver`) and `Capability::Hessian`.

## 12. Environment & portability findings (from the setup pass)

- Local env configured: OpenMP (libomp), MPI (Open MPI 5.0.9), BLAS/LAPACK
  (Accelerate), Rust 1.96, molrs at `../molrs`. Full CPU build (RHF+double) green.
- **Two hardcoded `jicli594` cluster paths** in `tests/CMakeLists.txt`
  (`ATV_CP2K_DATA_DIR`, `ATV_FINAL_FRAME_XYZ`) — **fixed** (opt-in + EXISTS-guard);
  the aborting `test_cpu_io_extxyz` now passes.
- **molrs robustness bug**: `xyz_read_first_frame` (`molrs-cxxapi/src/lib.rs:446`)
  `.expect()`s the read → a missing file **aborts the process** instead of
  surfacing the `std::runtime_error` that `load_extxyz`'s own tests expect. Fix
  belongs in molrs (separate repo, currently 4 commits ahead of origin/dev).
- **Precision knobs `ATV_REAL`/`ATV_INT` were dead** → replaced with LAMMPS-style
  tiered `core/precision.h` (`StorageReal`/`WorkReal`/`AccReal`,
  `LocalInt`/`TagInt`/`ImageInt`/`BigInt`). Migration of hot structures pending
  (see Task #9 / [[project-precision-tiered-lammps]]).
- `test_cpu_snapshot_publisher` fails only under `ctest -j` (parallel flake);
  passes standalone.

---

## Bottom line

The engine is **more capable and more unified than its own docs claim**: a real
4-category × multi-backend Driver/Component/Frame architecture; production
classical MD (incl. real MPI); PySCF-validated Gaussian HF/DFT/hybrid + gradients
+ vibrations; experimental-but-real plane-wave and GPW DFT; solid PME. The manual
can be built on a genuine foundation. The work is: (a) settle terminology + repair
developer-doc drift, (b) surface GPW + torch honestly, (c) **bridge molrs's
existing analysis library to C++** (FFI exposure, not new C++ code — see
[[feedback-infrastructure-in-molrs]]), and (d) build the missing benchmark
categories (backend-consistency, convergence, regression). No architectural
invention needed — align the manual to what already exists.
