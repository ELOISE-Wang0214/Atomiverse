---
slug: esmm-contract-01-zeroing
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: driver / pipeline-contract
queue: esmm-campaign
deferred: false
depends_on: []
unblocks: [hessian-classical, path-classical, esmm-contract-05-accumulate-switch]
grilled: true
---

# Spec: `esmm-contract-01-zeroing` — zeroing is the driver's job, in every driver

**Slug:** `esmm-contract-01-zeroing`
**Domain:** CPU driver pipeline (zeroing contract)
**Queue:** esmm-campaign (Wave 1 — composition contract)
**Status:** approved. Grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

`ZeroForces`/`ZeroEnergy`/`ZeroStress` are bound by MD, Evaluate, Minimize and
Phonon — but **not uniformly**: `vib::Hessian` binds neither Zero, `PathDriver`
binds neither, MC binds `ZeroEnergy` only, and Minimize/Evaluate bind no
`ZeroStress`. Three of these are **live bugs today** and all four block the
Wave-1 accumulate switch (`esmm-contract-05`):

1. **Hessian**: an accumulating classical `Potential` integrates garbage across
   the 6N displacements (`HessianAccumulator` snapshots raw frame spans,
   `hessian_accumulator.h:63-72`). The artifact is **fossilized in an oracle**:
   `regressions/architecture/harmonic_driver_matrix/reg_harmonic_driver_matrix.cpp:177-183`
   hedges "Scale may be k or k/2 depending on FD assembly convention" — the k/2
   is not a convention, it is exactly the central difference of running sums of
   a linear force.
2. **Minimize/Evaluate stress**: `PairLJCut` does `+=` on `kStress`, so the LJ
   cell relax (`tests/cpu/minimize/test_cell_relax.cpp:48-73`) accumulates
   stress across all 80 minimizer iterations, passing only because the oracle
   is `L1 > L0 && isfinite`.
3. **PathDriver**: no Zero anywhere on `path::Phase::force`; `NebCoupling`
   zeroes only the frozen endpoints (`neb.cpp:61-71`). NEB works today by
   accident of the regression-local `MullerBrown` stub being an overwriter
   (`reg_neb_muller_brown.cpp:70-73`). Any `+=` Potential on a band accumulates
   across band iterations.
4. **MC**: forces are unused but grow without bound once sources accumulate
   (NaN/overflow trap; debugging trap).

This spec lands the zeroing everywhere, tightens the fossilized oracle, and
updates the plan-shape anchor tests — **before** any evaluator switches from
`=` to `+=`.

## Domain basis

Target contract (binding for the whole esmm campaign): **the hosting driver
owns zeroing of the accumulator channels** (`kEnergy`, every `kEnergy*`
per-term key, `kAtom_fx/fy/fz`, `kStress` virial part). A force/energy source
never zeroes — it would erase a co-bound source. The Zero* components
deliberately declare no `default_bindings()` (`zero.h:5-9`): the driver picks
the phase.

Rejected alternative: auto-injecting Zero* from the base `Driver::build()` for
every driver — SCF has no force phase at all, and a blanket injection would
surprise nested method-drivers.

## Design

| Driver | Today | Change |
| ------ | ----- | ------ |
| `md::MD` | ZeroForces/Energy/Stress @pre_force −1000/−999/−998 (`md.h:155-158`) | none |
| `Evaluate` | ZeroForces/Energy @prepare (`evaluate.h:73-75`) | **add `ZeroStress` @prepare −998** |
| `minimizer::Minimize` | ZeroForces/Energy @force (`minimize.h:137-139`) | **add `ZeroStress` @force −998** |
| `vib::Phonon` | ZeroForces −1000 / ZeroEnergy −999 @vib force (`phonon.h:78-79`) | none (the model to mirror) |
| `vib::Hessian` | **nothing** (`hessian.cpp:49-75`) | **add ZeroForces −1000 / ZeroEnergy −999 @`vib::Phase::force`** in `wire_full_chain` AND the `WithoutAccumulator` ctor |
| `mc::McBase` | ZeroEnergy @energy −999 (`mc_base.h:146`) | **add ZeroForces @energy −1000** |
| `PathDriver` | **nothing** (`path_driver.h:88`) | **add ZeroForces −1000 / ZeroEnergy −999 @`path::Phase::force`** in the ctor (replays per replica inside the per-replica `run_stage` pass — exactly the right scope) |

Behavioral audit (must hold): the only existing case whose numbers change is
`harmonic_driver_matrix` (H00 goes k/2 → k, which the loosened oracle already
accepts) — tighten it in the same change. ES-based Hessian regressions
(`regressions/scf/*hessian*`, `h2o_rhf_freq`, `h2o_rhf_ir`) are sole
overwriting sources → unchanged. NEB regression: overwrite-after-zero →
unchanged. `test_cell_relax`: pressure now correct per iteration; assertions
still pass.

## Files to create or modify

- `src/cpu/vibration/hessian.cpp` — the two `on(vib::Phase::force)` Zero adds
  in `wire_full_chain` + `WithoutAccumulator` ctor; include `component/zero.h`
- `include/atomiverse/cpu/path/path_driver.h` — ctor Zero pair on
  `path::Phase::force`; include `component/zero.h`
- `include/atomiverse/cpu/montecarlo/mc_base.h` — `.add(cpu::ZeroForces{}, -1000)`
  before the existing ZeroEnergy on `sampling::Phase::energy`
- `include/atomiverse/cpu/minimize/minimize.h` — `.add(cpu::ZeroStress{}, -998)`
  on `minimizer::Phase::force`
- `include/atomiverse/cpu/evaluate/evaluate.h` — `.add(cpu::ZeroStress{}, -998)`
  on `prepare`
- `tests/cpu/test_driver_binding_kind.cpp` — update the three plan-shape
  anchors: `HessianPlanShapeUnchanged` (:237), `McPlanShapeUnchanged` (:265),
  `PathPlanShapeUnchanged` (:314)
- `regressions/architecture/harmonic_driver_matrix/reg_harmonic_driver_matrix.cpp`
  — delete the k/2 hedge (:154-156, :177-183); assert `H00 == k` (small
  `err_cart`); keep the isotropy check
- NEW `tests/cpu/vibration/test_hessian_zeroing.cpp` (unit U1a)
- NEW `tests/cpu/path/test_path_zeroing.cpp` (unit U1b)

All of these are **outside** the operator's in-flight dirty set (verified
2026-07-26: `hessian.cpp`, `path_driver.h`, `mc_base.h`, `minimize.h`,
`evaluate.h`, `runtime.cpp` clean).

## Tasks

- [ ] RED: `tests/cpu/vibration/test_hessian_zeroing.cpp` — an accumulate-style
  (`+=`) harmonic stub potential under `vib::Hessian` must yield `H00 == k`
  (fails today with k/2-flavored garbage)
- [ ] RED: `tests/cpu/path/test_path_zeroing.cpp` — two PathDriver iterations
  with a `+=` stub: iteration-2 replica forces equal iteration-1's at the same
  geometry (fails today: doubled)
- [ ] Add the five driver Zero bindings per the table
- [ ] Update the three plan-shape anchors in `test_driver_binding_kind.cpp`
- [ ] Tighten `harmonic_driver_matrix`: delete the k/2 hedge, assert the exact
  scale; run the regression and record H00 == k
- [ ] Audit run: full ctest + the vib/path/mc/minimize regression subset —
  confirm only `harmonic_driver_matrix` changed numbers

## Testing strategy

**Unit:** U1a (Hessian `+=` stub → H00 == k), U1b (PathDriver two-iteration
force idempotence). Both are RED first.

**Regression:** no new case — this spec's e2e witness is the **tightened**
`harmonic_driver_matrix` oracle (deleting a hedge that fossilized a bug is the
regression-side deliverable).

**Behavioral-equivalence gate:** whole existing suite green; only
`harmonic_driver_matrix` numbers move.

## Out of scope

- Any `=` → `+=` evaluator change (esmm-contract-05)
- run_slot_ nested dispatch (esmm-contract-03)
- MC validation walk (esmm-contract-02)
- CUDA drivers

## Notes for implementers

- Zero* zero **owned rows only**; ghosts are zeroed by
  `HaloSend::zero_ghost_forces` — do not touch that split.
- Keep the Zero orders exactly −1000/−999/−998 (mirrors MD/Phonon) so the
  plan-shape anchors stay canonical.
- Regression CMake is NOT touched here (no new case).
