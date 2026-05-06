---
title: Phase 2 — SinglePoint driver as first-class citizen
status: superseded
superseded_by: cpu-driver-taxonomy
created: 2026-04-14
superseded_at: 2026-05-06
---

> **Superseded by [`cpu-driver-taxonomy.md`](cpu-driver-taxonomy.md).** The
> SinglePoint driver design here is preserved verbatim, but lifted into a
> broader four-class taxonomy (Evaluate / Dynamics / Optimize / SCF). See
> Phase 3 of the successor spec for the SinglePoint scope.

# Phase 2 — SinglePoint Driver

## Summary

Promote "single-point energy" from a call to `scf_solve()` inside evaluators to a **first-class CPU driver type** (`cpu::scf::SinglePoint`) that composes with SCF-phase components. Mixed-level workflows ("optimize at A, single-point at B") are then expressed by the user **sequentially invoking two drivers at the call-site** — no composite wrapper class.

## Binding rules (apply here)

- **No `OptimizeThenSinglePoint` class. No `WorkflowRunner`. No `PipelineDriver`.** The word "workflow" appears nowhere in the public API.
- The only compositional primitive is sequential driver invocation in user code, each driver constructed with its own config.
- `SinglePoint` is a `Driver<scf::Phase, cpu::device_tag>` just like `cpu::minimizer::LBFGS` is a `Driver<minimizer::Phase, cpu::device_tag>`. Same contract: `add` / `on(phase).add` / `build` / `run`.
- PySCF `pyscf.scf.RHF(mol).kernel()` is the ground truth for single-point energies.

## Driver design

```cpp
namespace ATV::cpu::scf {

/// @brief Single-point electronic-structure driver.
///
/// Runs the SCF cycle to self-consistency over a single geometry.
/// Composes with SCF-phase components (orbital observers, cube writers,
/// convergence checkers). Does NOT own a geometry-optimization loop.
class SinglePoint : public Driver<ATV::scf::Phase, cpu::device_tag> {
public:
    explicit SinglePoint(SCFConfig config);

    /// Auto-added defaults: ZeroDensity (pre_hamiltonian), SCFStep, DIISMixer,
    /// StandardConvergenceCheck. User can override via add().
private:
    void default_bindings() override;

    SCFConfig config_;
};

} // namespace
```

Internal components that make up the default binding set (all new, one responsibility each):

| Component | Phase | Responsibility |
|---|---|---|
| `cpu::scf::ZeroDensity` | `build_hamiltonian` (order -1000) | Clear density matrix on entry |
| `cpu::scf::HamiltonianBuilder` | `build_hamiltonian` | Build `H = h_core + J - αK + V_xc` |
| `cpu::scf::EigenSolver` | `diagonalize` | `dsygvd` generalized eigen on `(F, S)` |
| `cpu::scf::DensityBuilder` | `build_density` | `P = 2 C_occ C_occ^T` (or αβ for UHF) |
| `cpu::scf::DIISMixer` | `build_density` (order +1) | Pulay subspace extrapolation |
| `cpu::scf::DualConvergenceCheck` | `convergence_check` | `|ΔE|<tol ∧ ‖err‖_F<√tol` |

The existing `scf_solve` monolithic free function survives **as an internal implementation helper only**; its call sites in `SCFEvaluator` migrate to build + run a `SinglePoint` instance instead. This ensures a single source of truth for the SCF cycle.

## Mixed-level workflow (no wrapper)

```cpp
// 1) Cheap optimization at PBE/def2-SVP
cpu::minimizer::LBFGS opt;
opt.add(cpu::minimizer::LBFGSStep{/*...*/});
opt.on(minimizer::Phase::force).add(
    cpu::scf::SCFEvaluator{cpu::scf::default_config(integrals_A)});
opt.on(minimizer::Phase::end_of_step).add(
    cpu::minimizer::TrajectoryWriter{"traj.xyz", 1});
opt.build();
opt.run(frame, 200);

// 2) Accurate single-point at B3LYP/6-31G* (Phase 3 functional; def2-TZVP interim)
auto integrals_B = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Def2TZVP);
cpu::scf::SinglePoint sp{cpu::scf::b3lyp_config(integrals_B)};  // Phase 3
sp.on(scf::Phase::end_of_step).add(cpu::scf::OrbitalInfoObserver{});
sp.build();
sp.run(frame, /*max_iter=*/100);
```

**That is the workflow.** Two driver constructions, two `run` calls. No glue class.

## SCFEvaluator refactor

Currently `SCFEvaluator::launch` calls `scf_solve(integrals, config)` directly. Refactor so it constructs a transient `SinglePoint` driver on first `launch()` (cached in `build()`), reuses density warm-start (ASPC) across `launch()` calls. This keeps the composition pattern consistent and eliminates the parallel code path.

## Acceptance

- **Correctness:** `SinglePoint{default_config}.run(frame)` produces identical (bit-for-bit) energy + density as current `scf_solve(integrals, default_config)` on all 12 benchmark molecules.
- **Correctness:** single-point energies at PBE/def2-SVP match PySCF `mf.kernel()` to ≤ 1e-7 Ha on all 12 benchmarks.
- **Correctness:** mixed-level workflow (PBE/def2-SVP opt → PBE/def2-TZVP SP) reproduces PySCF `pyscf.geomopt.geometric_solver.optimize` + second `RHF.kernel()` to ≤ 1e-5 Ha final SP energy on ≥ 5 molecules.
- **Performance:** `SinglePoint::run` adds ≤ 2% overhead vs the current `scf_solve` free-function path (verified by ctest benchmark).
- **Performance:** single-point PBE/def2-SVP on aspirin ≤ 1.2× PySCF wall time.
- **Architecture:** grep check — zero occurrences of `Workflow`, `Runner`, `Pipeline`, `Composite` class names in new code. Zero callers of `scf_solve` outside `SinglePoint::run` and transitional backward-compatibility tests.

## Files to create / modify

| File | Action |
|---|---|
| `include/atomiverse/cpu/scf/single_point.h` | new (driver class) |
| `src/cpu/scf/single_point.cpp` | new |
| `include/atomiverse/cpu/scf/components/*.h` | new (HamiltonianBuilder, EigenSolver, DensityBuilder, DIISMixer, DualConvergenceCheck, ZeroDensity) |
| `src/cpu/scf/components/*.cpp` | new |
| `src/cpu/scf/scf_solve.cpp` | demoted to helper; public wrapper deprecated |
| `src/cpu/scf/scf_evaluator.cpp` | refactor to hold `SinglePoint` member |
| `tests/cpu/scf/test_single_point_driver.cpp` | new (parity with scf_solve + PySCF) |
| `tests/cpu/scf/test_mixed_level_composition.cpp` | new (opt→SP sequencing) |
| `examples/tutorial/06_opt_then_sp.cpp` | new (demonstrates composition pattern) |

## Dependencies

- **Phase 1 required:** `BasisSet`, `build_integrals`, `default_config`, recorder components.

## Risks

- **Component granularity drift:** splitting `scf_solve` into six components risks over-decomposition. Mitigation: each component has an explicit `launch()` responsibility documented in the header; code reviewer (`av-architect`) validates no component straddles two.
- **ASPC predictor integration:** `SCFEvaluator` currently stashes predictor state in its own fields. Moving to `SinglePoint` requires the predictor to be either a stateful component or plumbed through the driver — design decision tracked in the spec implementation phase.

## Exit gate

All acceptance criteria green + `SCFEvaluator` successfully migrated to internal `SinglePoint` usage without test regression + architecture grep check clean + user sign-off.
