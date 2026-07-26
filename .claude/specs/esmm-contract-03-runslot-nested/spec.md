---
slug: esmm-contract-03-runslot-nested
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: driver / pipeline-contract
queue: esmm-campaign
deferred: false
depends_on: []
unblocks: [rhf-minimize-linesearch, rhf-in-neb]
grilled: true
---

# Spec: `esmm-contract-03-runslot-nested` — run_stage replays nested driver nodes

**Slug:** `esmm-contract-03-runslot-nested`
**Domain:** CPU Runtime stage replay (nested-node dispatch)
**Queue:** esmm-campaign (Wave 1 — composition contract)
**Status:** approved. Grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

`Runtime::run_slot_` — the replay path behind `run_stage` and therefore behind
Minimize's line-search callback and PathDriver's per-replica force pass —
**skips nested driver nodes**: `runtime.cpp:245-247` `if (op.is_driver())
continue;`. The consequence is worse than stale forces, it is **zeroed**
forces: `reeval_forces` (`minimize.h:260-263`) runs the Zero* components,
skips the nested SCF, runs Reduce — every line-search trial sees `E = 0`,
`f = 0`. StrongWolfe is the **default** line-search mode for both LBFGS and CG
(`lbfgs.h:169`, `cg.h:94`), so `examples/05_geom_opt.cpp` and
`tests/cpu/minimize/test_cg.cpp:28` have been limping on degenerate line
searches (the outer force slot rescues each iteration because
`Runtime::step` *does* dispatch nested drivers, `runtime.cpp:137-149`).
PathDriver's per-replica `run_stage(path::Phase::force, s)` has the same hole:
a nested-RHF band silently computes no forces.

Fix: dispatch driver ops in `run_slot_` exactly the way `step_plan` does.
Ships regression **R2** — the first geometry optimization in the tree whose
line search provably works over a nested QM source.

## Domain basis

Re-entrancy is verified safe: `Driver::run_nested` (`driver.h:715-742`) saves
and restores `ctx.loop` on the stack, so calling `run_nested_erased` from
inside an optimizer component's `launch()` (the `reeval_forces` path) nests
correctly and an inner loop still cannot stop its enclosing loop. Frame
column-resolution windows are saved/restored per component on every path
(`runtime.cpp:211-215, 249-253`), and `MeanFieldMethod` parks its own window
(`mean_field_method.cpp:868, 891`).

Timing: `run_slot_` deliberately never touches the timing registry (its own
comment, `runtime.cpp:238-241`) and `step_plan` records no nested timings
either (`runtime.cpp:179-181`). Nested-op timing is **explicitly deferred** —
extend those comments so the deferral is a recorded decision, not an
oversight.

## Design

Replace the skip with the `step_plan` dispatch shape:

```cpp
if (op.is_driver()) {
    static_cast<cpu::DriverNode*>(op.nested_driver)
        ->run_nested_erased(context_, *frame_, op.nested_budget);
    continue;
}
```

What then works:

- **Minimize line search** — φ(α) trials re-run Zero* → nested SCF → Reduce at
  the trial geometry; StrongWolfe/Armijo see real energies and slopes. Each
  trial costs one SCF; users who cannot afford it select
  `LineSearchMode::Fixed` (release-note this).
- **PathDriver per-replica pass** — nested RHF evaluates per replica (with
  `esmm-contract-01`'s per-replica Zero* first). The shared warm-start density
  now hops between adjacent images — a good guess (neighboring geometries) and
  correctness-neutral (SCF converges to the geometry's fixed point).
- **vib replay** — unchanged (Hessian/Phonon use `step_once` → `Runtime::step`,
  which already dispatched nested nodes).

Rejected alternative: teaching Minimize to special-case nested force sources —
the Runtime replay path is the single owner of slot execution; forking the
semantics per driver would leave PathDriver and every future `run_stage`
caller broken.

## Files to create or modify

- `src/cpu/driver/runtime.cpp` — `run_slot_` (`:237-256`): the dispatch; keep
  + extend the no-timing comments (recorded deferral) in both `run_slot_` and
  `step_plan`
- `tests/cpu/driver/test_driver_nesting.cpp` — new case
  "RunStageReplaysNestedDriver" (U2)
- `tests/cpu/minimize/` — new case: StrongWolfe line-search trials re-enter a
  nested stub node (U5)
- NEW `regressions/minimizer/rhf_h2_linesearch/{reg_rhf_h2_linesearch.cpp,REFERENCE.md}` (R2)
- `regressions/CMakeLists.txt` — register R2 under `ATV_BUILD_HF`
  (**gated step** — see Notes)

`runtime.cpp` and the test files are outside the operator's in-flight dirty
set (verified 2026-07-26); `regressions/CMakeLists.txt` is **inside** it.

## Tasks

- [ ] RED: `test_driver_nesting.cpp` "RunStageReplaysNestedDriver" — nested
  inner driver carrying a call-counting component, bound on the outer force
  phase; invoke the outer backend's `run_stage(force)`; assert the counter
  incremented (fails today: stays 0)
- [ ] RED (U5): Minimize + nested stub node + StrongWolfe — after `run`, the
  stub's call counter exceeds the outer iteration count (line-search trials
  re-entered it; fails today: equal)
- [ ] Implement the `run_slot_` dispatch; extend the timing-deferral comments
- [ ] GREEN both units + full ctest (existing minimize/path regressions must
  pass; `05_geom_opt` behavior may improve, not degrade)
- [ ] R2: write `regressions/minimizer/rhf_h2_linesearch` — H₂, RHF/STO-3G,
  LBFGS with default StrongWolfe; assert `r_e = 0.7122 Å` (CCCBDB RHF/STO-3G)
  within 5e-3 Å AND iteration count ≤ a bound only reachable with a working
  line search (record the measured bound in REFERENCE.md)
- [ ] REFERENCE.md: cite NIST CCCBDB (RHF/STO-3G H₂ equilibrium bond length
  0.712 Å), document the iteration-bound rationale
- [ ] Register R2 in regressions CMake (gate: file must be committed-clean
  first)

## Testing strategy

**Unit:** U2 (replay dispatch) + U5 (line-search re-entry), both RED-first on
stubs — no SCF in `tests/`.

**Regression (R2):** `minimizer/rhf_h2_linesearch` — public API only, ≤ 60 s,
literature value from CCCBDB; the iteration bound is the witness that trials
carry real energies (with E=0 trials the search degenerates and the count
blows past the bound).

## Out of scope

- Nested-op timing (recorded deferral)
- forces_requested plumbing (esmm-contract-04) — R2 runs with forces on
- Zeroing additions (esmm-contract-01) — independent
- CUDA runtime

## Notes for implementers

- HARD GATE for the CMake step only: `regressions/CMakeLists.txt` is part of
  the operator's uncommitted libcint/C++20 sweep (~760 files, 2026-07-26). If
  `git status` still shows it dirty when you reach that step, STOP and report
  — do not rebase over uncommitted operator work. All other edits in this
  spec are on clean files and may proceed.
- `op.nested_budget` already stores the bound budget (`core/plan.h`); `< 0`
  was resolved to `default_budget()` at bind time — do not re-resolve.
- Perf note for the release notes: StrongWolfe + nested SCF is now genuinely
  expensive per trial (it used to be free because it was wrong).
