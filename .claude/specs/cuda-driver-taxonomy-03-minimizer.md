---
title: CUDA Driver taxonomy 03 — minimizer migration
status: approved
created: 2026-05-22
follows: cuda-driver-taxonomy-02-md
---

# CUDA Driver taxonomy 03 — minimizer migration

## Summary

This is the third link of the five-spec chain that migrates the Atomiverse
CUDA backend onto the driver taxonomy already shipped on CPU. It removes the
algorithm-as-driver anti-pattern from the CUDA minimizer domain: the
`cuda::minimizer::SD` and `cuda::minimizer::CG` *driver* classes are deleted,
the steepest-descent step is promoted to a proper `Optimizer`-role component
renamed `cuda::minimizer::SteepestDescent` (matching the CPU `LBFGS` / `CG`
component naming), and a new concrete `cuda::minimizer::Minimize` driver
subclassing `cuda::OptimizeDriver` becomes the single user-facing CUDA
minimizer. Algorithms are now hot-swappable components plugged into one
driver, exactly mirroring `cpu::minimizer::Minimize`. After this lands, CUDA
energy minimization is `cuda::minimizer::Minimize` + a `SteepestDescent`
component, and the existing CUDA steepest-descent regression continues to
relax a perturbed LJ lattice to a minimum unchanged.

## Domain basis

Pure structural refactor — no equations, no new physics. The steepest-descent
update kernel (`x += alpha * f / m * kFtm2v`) is moved verbatim, not modified.
The minimization trajectory is preserved by construction: the same kernel,
the same `Zero*` infrastructure, and the same phase sequence
(`search_direction → line_search → force → convergence_check → end_of_step`)
are retained. Behavior preservation is verified by the existing CUDA
steepest-descent regression test (`test_cuda_minimizer_sd`), which asserts the
energy of a displaced LJ lattice decreases monotonically toward a minimum.

## Design

Affected domain: CUDA minimizer (`cuda/minimizer/`).

Current state (post-01/02): `cuda/minimizer/sd.cuh` holds both
`SteepestDescentStep` (a `ComponentBase`, `role()` defaulting to `None`,
search_direction phase) and an `SD` *driver* class that subclasses raw
`Driver<Phase>` and auto-adds `ZeroForces` / `ZeroEnergy`.
`cuda/minimizer/cg.cuh` holds only a `CG` *driver* class — an empty shell;
there is no CG step component and no CG kernel anywhere under
`src/cuda/minimizer/` (confirmed: `src/cuda/minimizer/` contains only
`sd.cu`).

Changes:

- **Delete the driver classes.** Remove `cuda::minimizer::SD` from `sd.cuh`.
  Delete `cuda/minimizer/cg.cuh` entirely — it carried only the empty `CG`
  driver shell with no kernel behind it; nothing else references CG.
- **Promote and rename the step component.** Rename `SteepestDescentStep` →
  `cuda::minimizer::SteepestDescent` for parity with the CPU `Optimizer`
  component names (`cpu::minimizer::LBFGS`, `cpu::minimizer::CG` are bare
  algorithm names without a `Step` suffix). Add
  `ComponentRole role() const override { return ComponentRole::Optimizer; }`.
  Keep its `default_bindings()` on `search_direction @ 0`, its `alpha`
  ctor parameter, and the `sd_step` kernel unchanged.
- **New `cuda::minimizer::Minimize` driver** in a new header
  `cuda/minimizer/minimize.cuh`. It subclasses `cuda::OptimizeDriver` (the
  category base introduced in spec 01). Its constructor binds the
  phase-agnostic `Zero*` infrastructure components (made phase-agnostic in
  spec 02) at `minimizer::Phase::force` via explicit `on(phase).add()`,
  mirroring `cpu::minimizer::Minimize`'s constructor. It adds no
  categorical-setter API — `cuda::OptimizeDriver` provides only the
  loop/validation core, and users register the `SteepestDescent` component
  and force evaluators via the inherited public `add()` / `on(phase).add()`.
- **Convergence.** CUDA has no `ConvergenceCheck` component, so
  `cuda::minimizer::Minimize` runs the full `max_iter` budget every time
  (`outer_status.converged` never set) — identical to the behavior of the
  deleted `SD` driver. A `ConvergenceCheck` component is a future addition.

Lifecycle / ownership: unchanged. `SteepestDescent` is a stack-constructed
component handed to the driver via `add()`; `Minimize` owns the registered
components through the underlying `Driver<>` plan, exactly as `SD` did.

## Files to create or modify

- `include/atomiverse/cuda/minimizer/sd.cuh` — delete the `SD` driver class;
  rename `SteepestDescentStep` to `SteepestDescent`; add `role()==Optimizer`;
  update the file-level and class-level Doxygen.
- `include/atomiverse/cuda/minimizer/cg.cuh` — delete the file entirely
  (empty `CG` driver shell, no kernel).
- `include/atomiverse/cuda/minimizer/minimize.cuh` (new) — new
  `cuda::minimizer::Minimize` driver subclassing `cuda::OptimizeDriver`,
  binding `Zero*` at `minimizer::Phase::force`.
- `src/cuda/minimizer/sd.cu` — update the `SteepestDescentStep::launch`
  definition to `SteepestDescent::launch`; update the file-level comment.
- `tests/cuda/minimizer/test_sd.cu` — port both tests to the
  `Minimize + SteepestDescent` shape; add a RED test asserting the new shape
  and the absence of the old driver classes.

## Tasks

- [ ] Write failing tests for the new CUDA minimizer shape in `tests/cuda/minimizer/test_sd.cu` — assert `cuda::minimizer::Minimize` runs a `SteepestDescent` Optimizer component to a relaxed geometry, and that the test compiles only against the new API (no `cuda::minimizer::SD` / `CG` driver classes, no `SteepestDescentStep`)
- [ ] Delete the `cuda::minimizer::SD` driver class from `include/atomiverse/cuda/minimizer/sd.cuh`
- [ ] Delete `include/atomiverse/cuda/minimizer/cg.cuh` entirely and remove any include of it from headers, sources, and `tests/`
- [ ] Rename `SteepestDescentStep` to `SteepestDescent` and add `role()==Optimizer` in `include/atomiverse/cuda/minimizer/sd.cuh`; update its Doxygen
- [ ] Update the kernel-launch definition from `SteepestDescentStep::launch` to `SteepestDescent::launch` in `src/cuda/minimizer/sd.cu`
- [ ] Implement `cuda::minimizer::Minimize` in new `include/atomiverse/cuda/minimizer/minimize.cuh`, subclassing `cuda::OptimizeDriver` and binding `Zero*` at `minimizer::Phase::force` in its constructor
- [ ] Add Doxygen per the project doc style (`@brief`, `@param` with units) to `SteepestDescent` and `cuda::minimizer::Minimize`
- [ ] Port both `tests/cuda/minimizer/test_sd.cu` tests to the `Minimize + SteepestDescent` shape and verify against the `test_cuda_minimizer_sd` regression (energy decreases to a minimum)
- [ ] Run full check + test suite (`ctest --test-dir build -R cuda_minimizer --output-on-failure`, plus the full build)

## Testing strategy

- **Happy path:** `cuda::minimizer::Minimize` with a `SteepestDescent`
  component and an LJ force evaluator (`NeighborList` + `LJ` bound at
  `minimizer::Phase::force`) relaxes a perturbed cubic LJ lattice — final
  energy is strictly less than the energy after a single force evaluation.
  This is the existing `Minimizer.SDEnergyDecreases` test ported to the new
  shape.
- **Behavior preservation:** the existing `test_cuda_minimizer_sd` regression
  must still pass under the new `Minimize + SteepestDescent` shape; the
  relaxed-geometry energy trajectory matches the pre-migration result (same
  kernel, same `Zero*`, same phase sequence).
- **Edge case — build-time validation:** constructing a `cuda::minimizer::Minimize`,
  adding only a force evaluator (no `Optimizer` component), and calling
  `build()` throws — `cuda::OptimizeDriver::validate_components()` requires
  ≥ 1 `Optimizer`-role component.
- **Edge case — thermo integration:** the existing `Minimizer.SDWithThermo`
  test ported to the new shape — `Thermo` bound at `end_of_step` prints
  `Step` / `PotEng` and the final energy is finite.
- **Compile-time shape guard:** the ported test references only
  `cuda::minimizer::Minimize` and `cuda::minimizer::SteepestDescent`; it
  fails to compile against the old `SD` / `CG` driver classes or the old
  `SteepestDescentStep` name.

Test command: `ctest --test-dir build -R cuda_minimizer --output-on-failure`
(requires a CUDA toolkit; the test target is `test_cuda_minimizer_sd`).

## Out of scope

- **CG as a component.** CUDA `CG` was only ever an empty driver shell with
  no kernel under `src/cuda/minimizer/`. Per the chain's no-new-CUDA-features
  rule, this spec does not invent a steepest-descent-to-CG algorithm — it
  ships `SteepestDescent` only. A `cuda::minimizer::ConjugateGradient`
  Optimizer component is deferred to a future spec.
- **A CUDA `ConvergenceCheck` component.** CUDA minimization currently has no
  convergence check; `cuda::minimizer::Minimize` runs the full `max_iter`
  budget exactly as the deleted `SD` driver did. Adding a `ConvergenceCheck`
  component that writes `outer_status` is a future addition.
- **Categorical-setter API for `cuda::minimizer::Minimize`.** `cuda::OptimizeDriver`
  provides only the loop/validation core; users register components via
  `add()` / `on(phase).add()`. A LAMMPS-shaped `set_optimizer` / `set_pair`
  API for CUDA is not introduced here.
- **CUDA MD and SCF migration** — links 02 (landed) and 04 of this chain.
