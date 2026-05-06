---
title: CPU SCF / Driver-Component model alignment
status: superseded
superseded_by: cpu-driver-taxonomy
created: 2026-05-06
approved: 2026-05-06
superseded_at: 2026-05-06
supersedes: phase2-singlepoint-driver (absorbs the SinglePoint driver design and extends it)
---

> **Superseded by [`cpu-driver-taxonomy.md`](cpu-driver-taxonomy.md).** The
> reframing introduces four driver category bases (Evaluate / Dynamics /
> Optimize / SCF), which subsumes the SCF-specific alignment originally
> targeted here. All work items from this spec survive there, distributed
> across phases 3-5 of the new spec.

# CPU SCF Architecture Alignment

## Summary

Bring the CPU SCF subsystem into structural parity with CPU MD and CPU Minimizer
under the `Driver<PhaseEnum, DeviceTag>` / Component / Frame contract documented in
`CLAUDE.md`. Today, CPU SCF is the only domain that **bypasses** the canonical
model: the SCF iteration is a free function (`scf_solve`), `scf::Phase` is dead
code on CPU, no `cpu::scf::ZeroDensity` exists, and `SCFEvaluator::launch`
allocates the entire integral pack + result on every call.

This spec drives three converging fixes, sequenced by ROI:

1. **`SCFEvaluator` no-alloc-in-launch** (highest ROI, no API change).
2. **`cpu::scf::SinglePoint` Driver** with the canonical phase components
   (absorbs `phase2-singlepoint-driver.md`; that spec is folded in here).
3. **`IntegralPack` discipline** — split fundamental quantities from derived
   pointer views; cleave oversized monolith files (`scf_solve.cpp`,
   `gradient.cpp`, `native_integrals.cpp`, `xc_functional.cpp`,
   `integral_derivatives.cpp`) toward the 200-400 line CLAUDE.md target.

When this lands, the four CPU domains (`md`, `minimizer`, `scf`, future `vib`)
share the same Driver / Component / Phase / Zero* / build-launch-destroy idiom.

## Domain basis

No new physics. Numerical ground truth stays at the existing oracles:
PySCF (RHF/UKS energies) and libcint (gradients). All refactors must preserve
bit-for-bit (1e-12 Hartree) parity on the existing benchmark molecules
(aspirin and the 12-molecule benchmark set introduced in phase 1).

## Architectural targets (acceptance contract — short form)

The detailed pass criteria live in `cpu-scf-driver-alignment.acceptance.md`.
Headline contract:

- **Driver symmetry.** `cpu::scf::SinglePoint : Driver<scf::Phase, cpu::device_tag>`
  exists, drives the SCF cycle through phase-bound components, and is the only
  call-path inside `SCFEvaluator`. `scf_solve()` survives only as an internal
  helper or is removed; no production caller invokes it directly.
- **Lifecycle.** `SCFEvaluator::build()` allocates the `IntegralPack` buffers,
  the SCF working set (Fock, density, DIIS history), and the gradient
  workspace. `SCFEvaluator::launch()` performs **zero** heap allocation
  on the steady-state path (verified by a counter or sanitizer hook test).
- **Infrastructure components.** `cpu::scf::ZeroDensity` exists, is auto-added
  by `SinglePoint` at `build_hamiltonian @ -1000`, and is documented in the
  CLAUDE.md "Implemented components" table.
- **Hot-path data discipline.** `IntegralPack` carries only fundamental
  physical quantities (S, T, V, Hcore, V_nuc, n_basis, n_electrons, atoms).
  Shell pointer tables / `shell_data` aliasing move into a derived view
  (`ShellList` / `ShellPairList`) constructed on demand. No raw pointer that
  aliases another field of the same struct.
- **File size.** No `src/cpu/scf/*.cpp` exceeds 800 lines. Targeted splits:
  `native_integrals.cpp` (1710), `rys_roots.cpp` (1455),
  `integral_derivatives.cpp` (1146), `rys_eri.cpp` (1043),
  `xc_functional.cpp` (639). `scf_solve.cpp` (413) becomes the implementation
  body of `SinglePoint`.
- **No regression.** Aspirin SCF wall time within ±5% of current on the
  reference machine; energy + gradient bit-for-bit identical to current
  baseline on all 12 benchmark molecules.

## Design

### Phase 1 — `SCFEvaluator` build/launch split (no API change)

Move all per-step heap allocation up to `build()`:

| Currently in `launch()` | Move to `build()` |
|---|---|
| `build_integrals(frame, basis, charge)` | Allocate `IntegralPack` once. Add `rebuild_integrals_inplace(frame, pack)` that overwrites buffers without reallocating xtensor storage (use `resize` only on first build, then `assign` / `noalias` going forward). |
| `SCFResult result;` (per call, owns several xtensors) | Pre-allocate one `SCFResult` member; `launch()` only updates fields. |
| `auto grad = scf_gradient(...)` | Add `scf_gradient_inplace(frame, ints, result, config, GradientWorkspace&)` that fills a pre-allocated `GradientWorkspace`. |
| Per-call `density_guess` xtensors | Reuse `last_P_alpha_` / `last_P_beta_` storage; `SCFInitialGuess` becomes a view, not an owner. |

Acceptance gate after Phase 1: a unit test asserts allocation count is constant
across `launch()` calls 2..N (e.g. via a simple `operator new` counter wrapper
in test mode, mirroring `tests/cpu/test_md_no_alloc.cpp` pattern).

### Phase 2 — `cpu::scf::SinglePoint` Driver (absorbs `phase2-singlepoint-driver.md`)

Promote the SCF cycle from a free function to a `Driver<scf::Phase, cpu::device_tag>`:

```cpp
namespace ATV::cpu::scf {

class SinglePoint : public Driver<ATV::scf::Phase, cpu::device_tag> {
public:
    explicit SinglePoint(SCFConfig config);
    // Auto-adds: ZeroDensity, HamiltonianBuilder, EigenSolver,
    //            DensityBuilder, DIISMixer, DualConvergenceCheck.
};

}
```

Component fan-out (one launch responsibility each, each
implements `build` / `launch` / `destroy`):

| Component | Phase | Order | Responsibility |
|---|---|---|---|
| `ZeroDensity` | `build_hamiltonian` | -1000 | Zero P on iteration 0 |
| `HamiltonianBuilder` | `build_hamiltonian` | 0 | F = Hcore + J - αK + V_xc |
| `EigenSolver` | `diagonalize` | 0 | dsygvd generalized eigen on (F, S) |
| `DensityBuilder` | `build_density` | 0 | P = 2 C_occ C_occ^T (or αβ) |
| `DIISMixer` | `build_density` | +1 | Pulay extrapolation |
| `DualConvergenceCheck` | `convergence_check` | 0 | \|ΔE\| < tol ∧ ‖err‖_F < √tol |

`SCFEvaluator` migrates from holding state to holding a built `SinglePoint`
member; `SCFEvaluator::launch` only reseats geometry-dependent inputs and
calls `single_point_.run(frame, max_iter)`. `scf_solve()` is reduced to a
shim that constructs a temporary `SinglePoint` and invokes `run()`, OR is
deleted from the public surface entirely (decision tracked in Task 2.7).

### Phase 3 — `IntegralPack` split + file granularity

Split the IntegralPack into a fundamental struct and a derived view:

```cpp
struct IntegralPack {                  // fundamentals only
    xt::xtensor<double, 2> S, T, V, Hcore;
    xt::xtensor<double, 4> ERI;        // optional, only when requested
    double V_nuc; int n_basis, n_electrons;
    std::vector<native::Atom> atoms;
};

struct ShellList {                     // derived view, computed from atoms+basis
    std::vector<native::Shell> shells; // pointers stable for the view's lifetime
    std::vector<std::vector<double>> shell_data;
    // Constructed by a free function `build_shell_list(atoms, basis)`;
    // never aliased into IntegralPack.
};
```

Consumers that today read `IntegralPack::shells` (Fock build, gradient,
basis_eval) take a `ShellList` parameter explicitly. `SCFEvaluator` owns one
`IntegralPack` and one `ShellList` member, both built once in `build()` and
refreshed in place each step.

File splits (target ≤ 800 lines each):

- `native_integrals.cpp` (1710) → `native_one_electron.cpp`,
  `native_two_electron_ssss.cpp`, `native_two_electron_general.cpp`,
  `native_dispatch.cpp`.
- `rys_roots.cpp` (1455) → `rys_roots_low_n.cpp`,
  `rys_roots_high_n.cpp` (or by Ltot bucket; `roots_for_x0.inl` stays as
  inline data).
- `integral_derivatives.cpp` (1146) → `derivatives_one_electron.cpp`,
  `derivatives_two_electron.cpp`, `derivatives_dispatch.cpp`.
- `rys_eri.cpp` (1043) → split by Ltot dispatch family.
- `xc_functional.cpp` (639) → `xc_lda.cpp`, `xc_gga.cpp`, `xc_dispatch.cpp`.

Splits are mechanical (move-and-rename), not redesigns; existing tests cover
the regression risk.

## Files

### Created

```
include/atomiverse/cpu/scf/single_point.h
include/atomiverse/cpu/scf/components/zero_density.h
include/atomiverse/cpu/scf/components/hamiltonian_builder.h
include/atomiverse/cpu/scf/components/eigen_solver.h
include/atomiverse/cpu/scf/components/density_builder.h
include/atomiverse/cpu/scf/components/diis_mixer.h
include/atomiverse/cpu/scf/components/dual_convergence_check.h
include/atomiverse/cpu/scf/shell_list.h
src/cpu/scf/single_point.cpp
src/cpu/scf/components/zero_density.cpp
src/cpu/scf/components/hamiltonian_builder.cpp
src/cpu/scf/components/eigen_solver.cpp
src/cpu/scf/components/density_builder.cpp
src/cpu/scf/components/diis_mixer.cpp
src/cpu/scf/components/dual_convergence_check.cpp
src/cpu/scf/shell_list.cpp
tests/cpu/scf/test_scf_evaluator_no_alloc.cpp
tests/cpu/scf/test_single_point_driver.cpp
tests/cpu/scf/test_zero_density_lifecycle.cpp
tests/cpu/scf/test_integral_pack_no_alias.cpp
```

### Modified

```
include/atomiverse/cpu/scf/scf_solve.h        (demote to internal helper or delete public surface)
include/atomiverse/cpu/scf/scf_evaluator.h    (hold SinglePoint + workspaces)
include/atomiverse/cpu/scf/integrals.h        (drop shells/shell_data from IntegralPack)
include/atomiverse/cpu/scf/gradient.h         (add scf_gradient_inplace + GradientWorkspace)
src/cpu/scf/scf_solve.cpp                     (becomes SinglePoint body or deleted)
src/cpu/scf/scf_evaluator.cpp                 (no-alloc launch)
src/cpu/scf/gradient.cpp                      (split + in-place variant)
src/cpu/scf/native_integrals.cpp              (split)
src/cpu/scf/rys_roots.cpp                     (split)
src/cpu/scf/integral_derivatives.cpp          (split)
src/cpu/scf/rys_eri.cpp                       (split)
src/cpu/scf/xc_functional.cpp                 (split)
CLAUDE.md                                     (add SinglePoint to driver table; ZeroDensity to components)
.claude/NOTES.md                              (record the alignment decision)
.claude/specs/INDEX.md                        (this entry; mark phase2-singlepoint-driver superseded)
.claude/specs/phase2-singlepoint-driver.md    (status → superseded, pointer here)
```

## Tasks

### Phase 1 — `SCFEvaluator` no-alloc (week 1)

- [ ] **1.1** Add `tests/cpu/scf/test_scf_evaluator_no_alloc.cpp`: build an
      evaluator, run 3 launches on a small molecule (H2O / STO-3G), assert
      that the per-launch new/delete count is constant on launches 2 and 3
      (mirror the pattern in `tests/cpu/test_md_no_alloc.cpp`). Test runs
      RED first.
- [ ] **1.2** Add `rebuild_integrals_inplace(frame, basis, charge, pack)` in
      `integrals.h`/`integrals.cpp` that overwrites an existing
      `IntegralPack` without reallocating its xtensors when shapes match.
- [ ] **1.3** Add `GradientWorkspace` + `scf_gradient_inplace(frame, ints,
      result, config, ws, out)` in `gradient.h`/`gradient.cpp`. Existing
      `scf_gradient` becomes a thin wrapper that owns a workspace.
- [ ] **1.4** Refactor `SCFEvaluator`: hoist `IntegralPack`, `SCFResult`,
      `GradientWorkspace`, density-guess buffers, and ASPC scratch into
      members allocated in `build()`. `launch()` becomes a sequence of
      in-place updates.
- [ ] **1.5** Test 1.1 turns GREEN. All existing AIMD / SCF tests still pass.
- [ ] **1.6** Benchmark: aspirin SCF (existing benchmark) within ±5% of
      current wall time on the reference machine.

### Phase 2 — `cpu::scf::SinglePoint` Driver (weeks 2-3)

- [ ] **2.1** Create `cpu::scf::ZeroDensity` infrastructure component
      (`include/atomiverse/cpu/scf/components/zero_density.h` +
      implementation). Default binding: `build_hamiltonian @ -1000`. Unit
      test (`test_zero_density_lifecycle.cpp`): only zeros on iteration 0,
      build/launch/destroy contract honoured.
- [ ] **2.2** Carve `HamiltonianBuilder`, `EigenSolver`, `DensityBuilder`,
      `DIISMixer`, `DualConvergenceCheck` out of the current `scf_solve.cpp`
      monolith. Each component lives in `cpu/scf/components/`. The phase-1
      no-alloc workspaces become member buffers of these components,
      allocated in their `build()`.
- [ ] **2.3** Implement `cpu::scf::SinglePoint`. Constructor auto-adds the
      six default components above; `add(...)` / `on(phase).add(...)` honour
      the canonical Driver fluent API.
- [ ] **2.4** Add `tests/cpu/scf/test_single_point_driver.cpp`:
      `SinglePoint{default_config}.run(frame, 100)` reproduces current
      `scf_solve(integrals, default_config)` energy + density bit-for-bit
      (1e-12 Ha) on aspirin and ≥ 5 of the 12-benchmark molecules.
- [ ] **2.5** Migrate `SCFEvaluator` to hold a `SinglePoint` member built
      once in `build()`; `launch()` reseats geometry inputs and calls
      `single_point_.run(...)`. ASPC predictor stays in `SCFEvaluator` (it
      is an MD-cadence concern, not an SCF-cycle concern).
- [ ] **2.6** Verify all AIMD tests, geometry-optimization tests, and the
      mixed-level workflow path still pass.
- [ ] **2.7** Decide: delete `scf_solve()` from the public header, or keep
      it as a one-line shim. Update `rhf_solve`/`uhf_solve`/`uks_solve`
      accordingly. Record the decision in `.claude/NOTES.md`.
- [ ] **2.8** Update `CLAUDE.md`: add `cpu::scf::SinglePoint` to the
      Drivers table; list `ZeroDensity` + the six SCF components.

### Phase 3 — `IntegralPack` split + file granularity (week 4)

- [ ] **3.1** Add `tests/cpu/scf/test_integral_pack_no_alias.cpp`: a copy of
      `IntegralPack` does not require pointer rebinding; no field of
      `IntegralPack` aliases another field of the same struct.
- [ ] **3.2** Introduce `cpu::scf::ShellList` (header + impl). Free function
      `build_shell_list(atoms, basis)` returns a value with stable internal
      pointers for its lifetime.
- [ ] **3.3** Move `shells` + `shell_data` out of `IntegralPack` into the
      new `ShellList`. Update every consumer
      (`operators.cpp`, `gradient.cpp`, `basis_eval.h`,
      `xc_grid_engine.h`, `SCFEvaluator`) to take `const ShellList&`
      explicitly. Test 3.1 turns GREEN.
- [ ] **3.4** Split `native_integrals.cpp` (1710 → ≤ 4 files, each ≤ 600
      lines). Build + tests still green.
- [ ] **3.5** Split `rys_roots.cpp` (1455 → ≤ 2 files; keep
      `roots_for_x0.inl` as data-only inline).
- [ ] **3.6** Split `integral_derivatives.cpp` (1146 → ≤ 3 files).
- [ ] **3.7** Split `rys_eri.cpp` (1043 → ≤ 2 files by Ltot family).
- [ ] **3.8** Split `xc_functional.cpp` (639 → 3 files: `xc_lda.cpp`,
      `xc_gga.cpp`, `xc_dispatch.cpp`).
- [ ] **3.9** Re-run the whole CPU SCF + AIMD test suite. Record final
      wall-time deltas (must remain within ±5% of pre-Phase-3 baseline).

### Phase 4 — Documentation + spec hygiene

- [ ] **4.1** Update `CLAUDE.md` "Implemented components" and "Conventions"
      sections to reflect SinglePoint, ZeroDensity, and the SCF component
      list.
- [ ] **4.2** Append a NOTE to `.claude/NOTES.md` recording: "CPU SCF is
      now a Driver-model citizen; `scf_solve()` is no longer the public
      entry point."
- [ ] **4.3** Mark `.claude/specs/phase2-singlepoint-driver.md` as
      `status: superseded`, pointer to this file.
- [ ] **4.4** Mark this spec `status: done` and remove from
      `INDEX.md` once all acceptance criteria pass.

## Testing strategy

- **Lifecycle.** `test_zero_density_lifecycle.cpp` —
  build/launch/destroy contract, only-zeros-on-iter-0 behaviour.
- **No-alloc.** `test_scf_evaluator_no_alloc.cpp` — operator-new counter
  asserts steady-state allocation count is constant across launches 2..N.
  Same pattern as `tests/cpu/test_md_no_alloc.cpp`.
- **Driver parity.** `test_single_point_driver.cpp` — bit-for-bit
  (1e-12 Ha) match with current `scf_solve` on aspirin + 5 of the 12
  benchmarks. Energy decomposition + final density both checked.
- **Pointer discipline.** `test_integral_pack_no_alias.cpp` — copy
  semantics, no internal aliasing.
- **Regression.** Existing `test_aimd_*`, `test_scf_*`, `test_uks_*`,
  `test_gradient_vs_libcint.cpp` all stay green.
- **Domain.** PySCF parity for SinglePoint at PBE/def2-SVP and HF/STO-3G.
  libcint parity for analytic gradient.

## Out of scope

- New physics (functionals, basis sets, post-HF) — none in this spec.
- CUDA SCF changes — `cuda::scf::SCF` already follows the contract; no
  edits planned.
- Mixed-level workflow / AIMD physics — these inherit the alignment
  automatically once SinglePoint lands.
- Performance optimization beyond "match current within ±5%" — perf gains
  are tracked in `eri-hpc-opt.md` and the SCF perf-gap memory.

## Risks

- **Component over-decomposition.** Splitting `scf_solve` into six
  components adds dispatch overhead. Mitigation: bench gate (Phase 2
  acceptance) ≤ 5% regression.
- **Pointer-rebinding bugs in IntegralPack split.** Catching the bug requires
  a copy-and-use test. Mitigation: Task 3.1 adds it RED-first.
- **File splits hide review intent.** Mechanical move-and-rename only — the
  reviewer needs to confirm no semantic edits hide in the diff. Mitigation:
  separate "moves" and "edits" commits.

## Exit gate

All Phase 1-4 tasks ticked, `cpu-scf-driver-alignment.acceptance.md`
criteria all green, CLAUDE.md updated, `phase2-singlepoint-driver.md`
marked superseded, INDEX entry removed, user sign-off recorded.
