> **Moved to `.claude/notes/notes.md`.** This file kept for backward reference.
> New content is written to `.claude/notes/notes.md` via `/mol:note` or `/av-note`.

# Architecture Notes (legacy copy)

Accumulated architectural decisions and corrections from development iterations.
Agents read this alongside CLAUDE.md. When a note stabilizes, promote it into CLAUDE.md.

---

## [2026-05-06] CPU driver taxonomy redesign — three categories + SCFSolver inner solver

Landed `.claude/specs/cpu-driver-taxonomy.md` (phases 1-6). Replaces the
single `Driver<PhaseEnum, DeviceTag>` template-with-algorithm-subclasses
scheme with a category-driven layer.

**Three CPU outer driver categories** (decided by termination semantic,
not by algorithm):

  - `cpu::EvaluateDriver` — one pass; concrete: `cpu::Evaluate`
  - `cpu::DynamicsDriver` — full step budget; concrete: `cpu::md::MD`
  - `cpu::OptimizeDriver` — early stop on convergence; concrete:
    `cpu::minimizer::Minimize` (renamed from `cpu::minimizer::LBFGS`)

The underlying `Driver<PhaseEnum, DeviceTag>` template is now a private
plan/dispatch primitive — user code subclasses one of the three
category bases via `cpu::CPUDriver<PhaseEnum>` (the CPU backend
layer, in `include/atomiverse/cpu/cpu_driver.h`). Subclassing
`Driver<>` directly is forbidden in CPU code. CUDA migration will
add a parallel `cuda::CUDADriver<PhaseEnum>` sibling.

**SCF is not an outer category.** `cpu::scf::SCFSolver` is a public
low-level inner-solver class (NOT a Driver subclass), used internally
by `cpu::scf::SCFEvaluator`. Its 5 step helpers
(`HamiltonianBuilder`, `EigenSolver`, `DensityBuilder`, `DIISMixer`,
`DualConvergenceCheck`) live in `cpu::scf::internal::` and never
inherit `cpu::ComponentBase`. Tests may include the internal headers
directly; production facade headers must not. There is no
`SinglePoint` class — single-point energy is `cpu::Evaluate +
cpu::scf::SCFEvaluator(.with_gradient(false))`, plus the optional
`cpu::single_point(frame, evaluator)` helper.

**Algorithm-as-component.** LBFGS / CG / FIRE are Optimizer
components plugged into `Minimize`. NoseHoover / Langevin /
Berendsen are Thermostat / Barostat components plugged into `MD`.
The new `cpu::minimizer::CG` (Polak-Ribière+) is hot-swappable with
`cpu::minimizer::LBFGS` inside the same `Minimize` driver.

**Component tagging — two orthogonal axes**:
  - `role()` — structural slot (`Integrator`, `Optimizer`,
    `Evaluator`, `Thermostat`, `Barostat`, `ConvergenceCheck`,
    `Recorder`, `Observer`, `None`)
  - `capabilities()` — bitset of frame fields written
    (`Energy | Force | Stress | Charge | Dipole | DensityMatrix |
    Orbitals`)

`DynamicsDriver` validates "≥ 1 component with `Capability::Force`"
regardless of role — `FixLangevin` (Thermostat role, Force-capable)
satisfies a pure-Langevin run.

**Outer iteration status channel**:
`cpu::Context::outer_status` is a `{converged, stop_requested,
reason}` struct, owned by the category driver, reset before each
iteration. ConvergenceCheck-tagged components write it; SCF
internals never touch it.

**Refactor was full no-backwards-compat**: `cpu::scf::DIIS` →
`cpu::scf::internal::DIISMixer` (renamed in place; one test caller
migrated). `cpu::minimizer::LBFGSStep` → `cpu::minimizer::LBFGS`.
`cpu::minimizer::LBFGS` driver → `cpu::minimizer::Minimize`. All 6
in-tree call sites updated in one commit per phase. No
`[[deprecated]]` shims.

**CUDA migration deferred**: `cuda::md::MD`, `cuda::minimizer::CG`/
`SD`, `cuda::scf::SCF` keep the per-algorithm shape until
`.claude/specs/cuda-driver-taxonomy.md` lands.

**`scf_solve()`** survives as a ≤30-line shim that constructs a
transient `SCFSolver` and runs it. Decision (vs deletion): retained
because `rhf_solve` / `uhf_solve` / `uks_solve` and the gradient
finite-difference path still call it; revisit during CUDA SCF
migration.

**Phase 5b deferred**: IntegralPack/ShellList split (drop `shells` +
`shell_data` from `IntegralPack` into a separate `ShellList` view)
and the 5 large file splits (`native_integrals.cpp` 1710 lines etc.)
remain TODO. Mechanical work; better as dedicated commits.

**Rule promoted to CLAUDE.md**: the three-category taxonomy, the
`role()`/`capabilities()` axes, and the "SCF is an inner solver, not
an outer category" decision.

## [2026-03-25] DeviceSystem removed — metadata lives in cuda::Frame

DeviceSystem class was deleted. Host-cached metadata (n_atoms, n_bonds, n_angles, n_dihedrals) now lives directly in `cuda::Frame` as plain Int fields with accessor methods. The device-side atom count pointer is `d_n_atoms` (not `n_atoms`). (Note: `DeviceFrame` was later renamed to `cuda::Frame`.)

**Rule**: Contexts (BuildContext, LaunchContext, DestroyContext) hold `cuda::Frame& frame` directly. Components access `ctx.frame.n_atoms()` for host-cached count and `ctx.frame.fx` etc. for device pointers. No wrapper class around `cuda::Frame`.

## [2026-04-06] CPU SCF force evaluators → Component refactoring

RHFForce, UHFForce, KSForce bypass the Driver/Component model — each owns a velocity Verlet loop (170 lines duplicated 3x), ASPC predictor logic (140 lines duplicated), and gradient skeleton (90 lines duplicated). Total ~590 lines of unnecessary duplication.

**Completed refactoring:**
1. ✓ Deduplicated gradient.cpp: extracted `accumulate_W`, `accumulate_common_gradient`, `make_result` helpers; unified FD gradient into single callback-based function
2. ✓ Unified `scf_gradient()` replacing separate rhf_gradient/uhf_gradient/ks_gradient — dispatches HF vs KS based on SCFConfig. Key subtlety: restricted HF stores P_alpha = P_total (factor 2 included), unrestricted stores per-spin densities
3. ✓ Created `SCFEvaluator` (cpu::ComponentBase at md::Phase::force) — single component for all methods, with optional ASPC predictor
4. ✓ Deleted legacy `RHFForce`, `UHFForce`, `KSForce` standalone classes
5. ✓ ASPC predictor lives inside SCFEvaluator (no separate wrapper needed)

**Rule promoted to CLAUDE.md:** "Driver-owns-the-loop" and "No duplication across methods" rules added to Component model section.

**Systematic bugs found during L-BFGS implementation:**
1. CUDA Driver classes use `Driver<Phase>` (1 arg) but template requires 2 `<PhaseEnum, DeviceTag>`. Latent bug — compiles only because CUDA sources not compiled in CPU-only builds. Fix deferred to CUDA session.
2. Infrastructure components (`md::ZeroForces`) hardcode `md::Phase::pre_force` in default_bindings — can't be reused by minimizer. Fix: create `cpu::minimizer::ZeroForces`/`ZeroEnergy` with bindings to `minimizer::Phase::force`.
3. CPU Driver subclass headers MUST include `cpu/backend_traits.h` for the `backend_traits<cpu::device_tag>` specialization. Without it, `this->add()` / `this->on()` fail at template instantiation.

## [2026-04-04] molrs integration complete

molrs CXX FFI integration is done. Key decisions that diverged from original spec:
- No `MolFrame` type in C++ -- all frame data lives in `molrs::MolRec` on the Rust side
- `cpu::MolRec` uses PIMPL (not direct `rust::Box`), keeping `rust/cxx.h` out of public headers
- Bridge is precision-templated (`bridge.rs` generated from template via `build.rs`)
- `ZarrWriter` delegates to `molrs::MolRec::write_zarr` (no separate `zarrs_ffi` dependency)
- Mulliken analysis stays in C++ (xtensor-based), not migrated to molrs
- `ATV_HAS_ZARR` / `ATV_BUILD_ZARR` flags removed -- Zarr is always available via molrs

## [2026-04-06] Test directory mirrors src/ layout

Tests are structured as `tests/{cpu,cuda,torch,facade}/...` mirroring the `src/` directory layout. Each test file uses GoogleTest fixtures (`TEST_F`) to organize tests by component/driver.

**Directory mapping:**
- `tests/cpu/md/` → CPU MD components (VelocityVerlet, LJ, CellList, bonded)
- `tests/cpu/scf/` → CPU SCF (integrals, solvers, gradient, grid, XC, evaluator)
- `tests/cpu/minimizer/` → CPU minimizer (L-BFGS)
- `tests/cuda/{md,potentials,minimizer,scf}/` → CUDA backend tests
- `tests/facade/` → Plan builder tests
- `tests/torch/` → LibTorch backend tests

**CMake naming convention:** `test_{backend}_{domain}_{component}`, e.g. `test_cpu_scf_gradient`, `test_cuda_md_thermostat`. Helper functions `atv_add_test(name, source_path)` and `atv_add_cpp_test(name, source_path)` take explicit source paths.

**Rule**: New test files go in the directory matching their source module. Use `TEST_F` fixtures to group tests by component. Each fixture class = one component or driver.
