---
title: CUDA Driver taxonomy — migrate cuda::md::MD / cuda::minimizer::CG / SD / cuda::scf::SCF to category bases
status: planned
created: 2026-05-06
follows: cpu-driver-taxonomy
---

# CUDA Driver Taxonomy

## Summary

`cpu-driver-taxonomy.md` (now `done`) introduced three CPU outer-driver
categories — `cpu::EvaluateDriver`, `cpu::DynamicsDriver`,
`cpu::OptimizeDriver` — plus the `cpu::scf::SCFSolver` inner solver,
component role/capability tagging, and the
`cpu::Context::outer_status` iteration channel. CUDA was deliberately
left out of scope.

This spec migrates the CUDA backend to the same three-category
taxonomy. After it lands, the user-facing facade is symmetric across
backends:

| Category | CPU concrete | CUDA concrete (after this spec) |
|---|---|---|
| `EvaluateDriver` | `cpu::Evaluate` | `cuda::Evaluate` |
| `DynamicsDriver` | `cpu::md::MD` | `cuda::md::MD` |
| `OptimizeDriver` | `cpu::minimizer::Minimize` | `cuda::minimizer::Minimize` |

CUDA SCF moves to the same inner-solver pattern: a `cuda::scf::SCFSolver`
class, used internally by a `cuda::scf::SCFEvaluator` ForceProvider
component. The current per-algorithm CUDA driver shape is dismantled:

  - `cuda::md::MD` (currently `Driver<Phase>` directly) → `DynamicsDriver`
  - `cuda::minimizer::SD` / `CG` (algorithm-as-driver) → `Minimize` +
    `SD` / `CG` Optimizer components
  - `cuda::scf::SCF` (currently `Driver<Phase>`, the only outer SCF
    driver) → demoted to `cuda::scf::SCFSolver` (inner) + use
    `cuda::Evaluate + cuda::scf::SCFEvaluator` for outer single-point.

## Binding rules (apply here, inherited from cpu-driver-taxonomy)

- Three current outer driver categories. No fourth without distinct
  outer-loop semantic.
- No `cuda::scf::SinglePoint`. Single-point energy is `cuda::Evaluate
  + cuda::scf::SCFEvaluator`.
- Algorithms (SD, CG, FIRE) become Optimizer components — never
  driver subclasses.
- SCF internals live in `cuda::scf::internal::` and are not
  registered into any outer Plan.
- Capability-over-role for evaluator validation. CUDA components get
  the same `role()` / `capabilities()` virtuals.

## Out of scope (deliberate)

- Performance optimization beyond "match current within ±5%" wall
  time on the CUDA reference machine.
- New CUDA features. This is purely a structural migration.
- Python bindings. Tracked separately.

## Pre-existing latent bugs to fix during migration

`.claude/NOTES.md` records two CUDA-side issues observed during the
CPU work:

  1. CUDA Driver classes use `Driver<Phase>` (1 template arg) but the
     `Driver<PhaseEnum, DeviceTag>` template requires 2. This compiles
     today only because CUDA sources are not exercised in the CPU-only
     CI. Fix: provide explicit `<Phase, cuda::device_tag>` instantiation
     OR migrate cleanly to `cuda::DynamicsDriver` / `OptimizeDriver` /
     `EvaluateDriver` bases that fix the second template arg.

  2. `cuda::md::ZeroForces` etc. hardcode `md::Phase::pre_force` in
     `default_bindings()` — fine for `cuda::md::MD` but not reusable by
     a CUDA minimizer/Evaluate driver. Mirror the CPU resolution: each
     concrete driver auto-adds Zero* via explicit `on(phase).add()` so
     the components are phase-agnostic.

Both should be addressed as part of this migration, not separately.

## High-level phases (sketch — full task list deferred to in-progress)

1. **Foundations.** Add `cuda::CategoryDriverBase<PhaseEnum>` and the
   three CUDA category bases (`cuda::EvaluateDriver` /
   `DynamicsDriver` / `OptimizeDriver`). Tag CUDA components with
   `role()` / `capabilities()` (mirrors CPU Phase 2).
2. **MD migration.** `cuda::md::MD` parent → `DynamicsDriver`. Fix
   the latent `Driver<Phase>` template-arg issue above.
3. **Minimizer migration.** Rename `cuda::minimizer::SD` and
   `cuda::minimizer::CG` from drivers to Optimizer components. New
   `cuda::minimizer::Minimize` driver.
4. **SCF migration.** Demote `cuda::scf::SCF` to
   `cuda::scf::SCFSolver` (inner solver, public). Add
   `cuda::scf::SCFEvaluator` ForceProvider component. Add
   `cuda::Evaluate` for single-point.
5. **Scrub `scf::Phase`.** With CPU and CUDA both off it for outer
   work, the legacy enum becomes unused — delete it from
   `facade/stage.h` (was kept for one release window).
6. **Docs / cleanup.** Update CLAUDE.md drivers table, NOTES entry,
   INDEX status → done.

## Acceptance (sketch)

  - Three CUDA category bases inherit `Driver<PhaseEnum,
    cuda::device_tag>`.
  - No `cuda::scf::SinglePoint`. No `cuda::scf::SCF` driver class.
  - All existing CUDA tests pass (if CUDA toolkit present in CI).
  - `scf::Phase` removed from `facade/stage.h` once CUDA migration
    is in.
  - PyTorch backend (if `ATV_BUILD_TORCH=ON`) — if it depends on the
    Driver template directly, a lightweight migration parallel to
    CUDA.

## Status

**Planned.** No driver work yet. CPU side is done and committed
(commits `588ac8a`, `<phase 2 sha>`, `<phase 3 sha>`, `<phase 4 sha>`,
`<phase 5a sha>`, `<phase 6 sha>`).

This spec is a *placeholder* until someone with CUDA toolkit access
(and time) picks it up. Promote `status: planned` → `draft` when
work begins.
