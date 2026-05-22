---
title: CUDA driver taxonomy 04 — SCF as inner solver
status: approved
created: 2026-05-22
follows: cuda-driver-taxonomy-03-minimizer
---

# CUDA driver taxonomy 04 — SCF as inner solver

## Summary

This is the fourth link of the five-spec chain that migrates the Atomiverse CUDA backend onto the driver taxonomy already shipped on CPU. It demotes the CUDA SCF *driver* (`cuda::scf::SCF`, currently a `Driver<scf::Phase, cuda::device_tag>` subclass) into a public low-level inner solver `cuda::scf::SCFSolver` that owns the SCF self-consistency cycle internally — exactly mirroring `cpu::scf::SCFSolver`. The SCF step components (Hamiltonian, EigenSolver, DensityBuilder, RHFFock/RHFFockDevice, RHFDensity, ConvergenceCheck) stop being `Plan`-registered `ComponentBase` subclasses and become `cuda::scf::internal::` helpers owned by `SCFSolver`, keeping their CUDA kernels unchanged. A new `cuda::scf::SCFEvaluator` (a real `Evaluator`-role component) wraps an `SCFSolver`, and a new one-pass `cuda::Evaluate` driver plus a `cuda::single_point` helper make single-point CUDA SCF energy `cuda::Evaluate + cuda::scf::SCFEvaluator` — with no `cuda::scf::SinglePoint` and no surviving `cuda::scf::SCF` driver class. This is a pure structural refactor: SCF numerics are untouched and the existing CUDA SCF regression suite remains the proof of correctness.

## Domain basis

No new equations. The SCF self-consistency cycle (effective Hamiltonian build → diagonalize → density build → mixing → convergence check) and the RHF energy expression `E = ½ Tr[(Hcore + F) P] + V_nuc` are carried verbatim from the existing CUDA SCF components into `cuda::scf::SCFSolver` and its `internal::` helpers; no kernel is rewritten. The SCF energy and density results are preserved by construction — the same kernels run in the same order. Correctness is verified by the existing CUDA RHF regression tests (`test_cuda_scf_rhf`) against textbook Hartree-Fock reference energies (H2/STO-3G `-1.1167` Ha, He `-2.8078` Ha, HeH+ `-2.8609` Ha, LiH `-7.8634` Ha) and by the synthetic `test_cuda_scf` convergence/parity tests.

## Design

Entities touched and their post-migration shape:

- **`cuda::scf::SCFSolver`** (new, `include/atomiverse/cuda/scf/scf_solver.cuh`) — a public low-level class, **not** a `Driver<>` subclass. Construction → `build(reference)` (one-time device-workspace allocation: H0/H/P/P_old/C/eig, cuBLAS/cuSOLVER handles, X = S^{-1/2} for the RHF path) → `run(...)` which loops the SCF cycle internally to convergence or a max-iteration cap. It owns the `internal::` step helpers as members, runs them in fixed order, and zeroes the density matrix internally at the start of a run. It never references `scf::Phase` and never touches `cuda::Context::outer_status`. Mirrors `cpu::scf::SCFSolver`. Two run paths are preserved: the synthetic effective-Hamiltonian path (`Hamiltonian`) and the RHF integral path (`RHFFock`/`RHFFockDevice` + `RHFDensity`) — surfaced as constructor variants / configuration so both regression suites map onto `SCFSolver`.
- **`cuda::scf::internal::`** (new namespace, new directory `include/atomiverse/cuda/scf/internal/`) — `Hamiltonian`, `EigenSolver`, `DensityBuilder`, `RHFFock`, `RHFFockDevice`, `RHFDensity`, `ConvergenceCheck` move here. They stop inheriting `cuda::ComponentBase`, drop `default_bindings()` / `role()` / `name()`, keep `build`/`step`/`destroy`-style lifecycle methods and their CUDA kernels. They are owned by `SCFSolver`, never registered into any `Plan`, and must not be included from any facade-public production header (tests may include them directly).
- **`cuda::scf::SCFEvaluator`** (new, `include/atomiverse/cuda/scf/scf_evaluator.cuh`) — a real `cuda::ComponentBase` with `role() == ComponentRole::Evaluator` and `capabilities() == Energy | DensityMatrix | Orbitals`. It holds an `SCFSolver` member, builds it once in `build()`, and runs it in `launch()`, writing energy/density into the frame. **No `Force` capability and no `with_gradient`** — no CUDA SCF gradient kernel exists in `src/cuda/scf/` (gradient-on-CUDA is Out of scope).
- **`cuda::Evaluate`** (new, `include/atomiverse/cuda/evaluate.cuh`) — a concrete one-pass driver subclassing `cuda::EvaluateDriver` (the category base shipped by link 01). It is a thin concrete driver, not a new category; mirrors `cpu::Evaluate` minus CPU-only categorical setters that do not yet exist on CUDA.
- **`cuda::single_point`** (new, `include/atomiverse/cuda/single_point.cuh`) — a free function template (not a class) that constructs a `cuda::Evaluate`, binds the borrowed evaluator, builds, and runs. Mirrors `cpu::single_point`.
- **`ZeroDensity`** — removed as a `Plan`-auto-added infrastructure component. Density zeroing becomes an internal first step of `SCFSolver::run`.
- **`cuda::scf::SCF` driver and `cuda/scf/scf.cuh`** — deleted. After this spec, `scf::Phase` (in `facade/stage.h`) is referenced by nothing outside the legacy enum definition itself, so link 05 can delete it cleanly.

Lifecycle / ownership: `SCFEvaluator` owns one `SCFSolver`; `SCFSolver` owns the `internal::` helpers and all device buffers. Device allocation happens once in `SCFEvaluator::build()` → `SCFSolver::build()`; `launch()` performs no allocation. The chain decision that `cuda::Evaluate` validation accepts any `Evaluator`-role component with non-empty `capabilities()` is inherited from link 01's `cuda::EvaluateDriver`.

## Files to create or modify

- `include/atomiverse/cuda/scf/scf_solver.cuh` (new) — `cuda::scf::SCFSolver` declaration.
- `include/atomiverse/cuda/scf/scf_evaluator.cuh` (new) — `cuda::scf::SCFEvaluator` declaration.
- `include/atomiverse/cuda/evaluate.cuh` (new) — `cuda::Evaluate` concrete one-pass driver.
- `include/atomiverse/cuda/single_point.cuh` (new) — `cuda::single_point` free-function helper.
- `include/atomiverse/cuda/scf/internal/hamiltonian.cuh` (new — moved from `cuda/scf/hamiltonian.cuh`).
- `include/atomiverse/cuda/scf/internal/eigen_solver.cuh` (new — moved from `cuda/scf/eigen_solver.cuh`).
- `include/atomiverse/cuda/scf/internal/density_builder.cuh` (new — moved from `cuda/scf/density_builder.cuh`).
- `include/atomiverse/cuda/scf/internal/rhf_fock.cuh` (new — moved from `cuda/scf/rhf_fock.cuh`).
- `include/atomiverse/cuda/scf/internal/rhf_density.cuh` (new — moved from `cuda/scf/rhf_density.cuh`).
- `include/atomiverse/cuda/scf/internal/convergence_check.cuh` (new — moved from `cuda/scf/convergence_check.cuh`).
- `include/atomiverse/cuda/scf/scf.cuh` — deleted (`cuda::scf::SCF` driver + `ZeroDensity` removed).
- `include/atomiverse/cuda/scf/hamiltonian.cuh`, `eigen_solver.cuh`, `density_builder.cuh`, `rhf_fock.cuh`, `rhf_density.cuh`, `convergence_check.cuh` — deleted (relocated under `internal/`).
- `src/cuda/scf/scf_solver.cu` (new) — `SCFSolver` implementation (the internal SCF loop; absorbs `ZeroDensity`'s `cudaMemsetAsync` density-zeroing).
- `src/cuda/scf/scf_evaluator.cu` (new) — `SCFEvaluator` implementation.
- `src/cuda/scf/driver.cu` — deleted (`ZeroDensity` source removed; logic folded into `scf_solver.cu`).
- `src/cuda/scf/hamiltonian.cu`, `eigen_solver.cu`, `density_builder.cu`, `convergence_check.cu`, `rhf_fock.cu`, `rhf_density.cu` — updated to `cuda::scf::internal::` namespace; kernels unchanged.
- `src/cuda/scf/device_integrals.cu` — namespace updated if it currently lands in the relocated namespace; kernels unchanged.
- `tests/cuda/scf/test_scf.cu` — ported to `cuda::scf::SCFSolver` (and `cuda::Evaluate + SCFEvaluator` for the single-point assertion); add RED guard that `cuda::scf::SCF` no longer exists.
- `tests/cuda/scf/test_rhf.cu` — ported to `cuda::scf::SCFSolver` / `cuda::Evaluate + cuda::scf::SCFEvaluator`; add RED single-point RHF-energy assertion.
- `tests/cuda/scf/test_device_integrals.cu` — updated includes/namespaces for the relocated `internal::` headers.
- `CMakeLists.txt` — swap `src/cuda/scf/driver.cu` → `src/cuda/scf/scf_solver.cu` + `src/cuda/scf/scf_evaluator.cu` in the CUDA SCF source lists (lines ~225-229, ~277-279).

## Tasks

- [ ] Write failing tests in `tests/cuda/scf/test_rhf.cu` and `tests/cuda/scf/test_scf.cu`: assert single-point CUDA RHF energy via `cuda::Evaluate + cuda::scf::SCFEvaluator` matches the textbook HF reference, and add a compile-time / type-level guard that `cuda::scf::SCF` driver class no longer exists
- [ ] Implement `cuda::scf::SCFSolver` in `include/atomiverse/cuda/scf/scf_solver.cuh` + `src/cuda/scf/scf_solver.cu` — internal SCF loop, internal density zeroing, no `scf::Phase`, no `outer_status`
- [ ] Move SCF step components to `cuda::scf::internal::` (new `include/atomiverse/cuda/scf/internal/` headers + namespace-updated `src/cuda/scf/*.cu`), dropping `ComponentBase` inheritance and `default_bindings()`; delete `cuda/scf/scf.cuh` and `src/cuda/scf/driver.cu`
- [ ] Implement `cuda::scf::SCFEvaluator` in `include/atomiverse/cuda/scf/scf_evaluator.cuh` + `src/cuda/scf/scf_evaluator.cu` — `role()==Evaluator`, `capabilities()==Energy|DensityMatrix|Orbitals`, owns one `SCFSolver`
- [ ] Implement `cuda::Evaluate` one-pass driver in `include/atomiverse/cuda/evaluate.cuh` inheriting `cuda::EvaluateDriver`
- [ ] Implement `cuda::single_point` free-function helper in `include/atomiverse/cuda/single_point.cuh`
- [ ] Update `CMakeLists.txt` CUDA SCF source lists: replace `src/cuda/scf/driver.cu` with `src/cuda/scf/scf_solver.cu` and `src/cuda/scf/scf_evaluator.cu`
- [ ] Port `tests/cuda/scf/test_scf.cu`, `test_rhf.cu`, `test_device_integrals.cu` to `cuda::scf::SCFSolver` / `cuda::Evaluate + cuda::scf::SCFEvaluator` and the relocated `internal::` headers
- [ ] Add Doxygen per project doc style (`@brief`, `@param` with units) to every new public symbol in `scf_solver.cuh`, `scf_evaluator.cuh`, `evaluate.cuh`, `single_point.cuh`
- [ ] Run full check + test suite (`ctest --test-dir build -R cuda_scf --output-on-failure`)

## Testing strategy

- **Happy path** — single-point CUDA RHF energy for H2/STO-3G via `cuda::Evaluate + cuda::scf::SCFEvaluator` reproduces the textbook reference `-1.1167` Ha within the suite's existing tolerance; `cuda::single_point` produces the identical result as constructing `cuda::Evaluate` and calling `add()` directly.
- **Numeric preservation (regression)** — the full pre-existing CUDA SCF regression suite passes unchanged after the migration: `test_cuda_scf` (synthetic convergence, reproducibility, alpha-effect, CPU/CUDA parity at n=8/16), `test_cuda_scf_rhf` (H2/He/HeH+/LiH RHF energies vs textbook + CPU/CUDA parity), `test_cuda_scf_device_integrals` (GPU integral kernels). Energies must match the pre-migration values to the suite's existing tolerances — proof that no SCF numerics changed.
- **Edge cases** — `cuda::scf::SCFSolver::run` honours its max-iteration cap when convergence is not reached (synthetic high-alpha case); reproducibility test confirms two independent `SCFSolver` runs on identical input give bit-identical energy; `cuda::Evaluate` validation rejects a frame with no `Evaluator`-role component (inherited `cuda::EvaluateDriver` behaviour).
- **Structural guards** — type-level / compile-time assertions: `cuda::scf::SCF` driver class does not exist; no `cuda::scf::SinglePoint`; `cuda::scf::SCFSolver` is not a `Driver<>` subclass; `cuda::scf::SCFEvaluator::role()` returns `Evaluator`; `cuda::Evaluate` derives from `cuda::EvaluateDriver`; `cuda::scf::internal::` helpers are not `cuda::ComponentBase` subclasses.

## Out of scope

- **CUDA SCF analytic gradient / `Force` capability** — no CUDA SCF gradient kernel exists in `src/cuda/scf/`; `cuda::scf::SCFEvaluator` is energy/density only, with no `with_gradient` toggle. Adding a CUDA SCF gradient is a future feature spec, not part of this structural migration.
- **Deletion of the legacy `scf::Phase` enum** from `facade/stage.h` — this spec leaves `scf::Phase` referenced by nothing outside its own definition; the actual deletion is link 05 (`cuda-driver-taxonomy-05-scrub-and-docs`).
- **CUDA categorical setter API** (`set_pair` / `set_neighbor` etc. on `cuda::Evaluate`) — CPU `Evaluate` inherits these from `CPUDriver`; the equivalent CUDA categorical surface is not built here. `cuda::Evaluate` exposes only `add()` / `on(phase).add()`.
- **ASPC predictor-corrector warm-start** — `cpu::scf::SCFEvaluator` has `predictor()`; the CUDA `SCFEvaluator` does not gain it here (no behavioural feature additions in a structural migration).
- **Changes to SCF numerics, kernels, or convergence algorithm** — kernels are relocated verbatim.
