---
slug: esmm-cpmd-06-single-eval
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: electronic-propagation / cpmd
queue: esmm-campaign
deferred: false
depends_on: [esmm-cpmd-03-observability]
unblocks: []
grilled: true
---

# Spec: `esmm-cpmd-06-single-eval` — one KS evaluation per GPW CPMD step

**Slug:** `esmm-cpmd-06-single-eval`
**Domain:** GPW CPMD performance (CPU)
**Status:** approved (grilled via plan-mode design review + operator Q&A, 2026-07-26)

## Summary

Every GPW CPMD production step currently evaluates the fixed-density KS
matrix **twice**: once only to obtain S for `align_to_overlap` (Γ RKS
`gpw_evaluator.cpp:~568-570`, UKS `:~516-524`, multi-k `:~389-397`), then
again at the aligned coefficients. Verified: nothing between the two evals
reads anything from the first result except `.S` — the first eval's H and E
are fully discarded. S is geometry-only. Build it standalone, align, then run
**one** KS evaluation: mathematically identical, and **bit-identical** if the
extraction reuses the same basis-assembly loop (same image counts, same
summation order). Expected win ≈ 2× on the dominant collocate→Poisson→XC→
integrate kernel of every production step.

## Design

### New standalone overlap on `GpwDft` (`gpw/gpw_solver.h`)

```cpp
/// @brief Standalone overlap S at the given geometry — basis assembly only,
///        no collocation/Poisson/XC. Bit-identical to the S a fixed-density
///        run() produces at the same geometry (same image counts and
///        summation order) — CPMD aligns the propagated C against it BEFORE
///        the step's single KS evaluation. @return S [nb × nb].
linalg::RealMatrix overlap(const std::vector<GpwAtom>& atoms) const;
/// @brief Multi-k overlaps S_k (Hermitian nb × nb, one per sampled k).
std::vector<linalg::CplxMatrix> overlap_k(const std::vector<GpwAtom>& atoms) const;
```

Extracted from `run()`'s own basis-assembly + S loop (Γ path over
`internal::overlap(a, b, L, m, reach)` in `one_electron.h`; k path via
`internal/kpoints.h::build_overlap_kinetic_k`) with identical image counts
and summation order — extraction, not re-implementation.

### CPMD step restructure (all three branches: Γ RKS, Γ UKS, multi-k RKS)

Old: eval#1 (discarded H/E, keeps S) → align → eval#2 → energy/forces →
advance. New: `S = solver_->overlap(atoms)` → `align_to_overlap(S)` →
**one** fixed-density eval → energy/forces → `advance(H, S, dt, occ)`.
The meta-key writes from esmm-cpmd-03 ride the single eval unchanged.

## Files to create or modify

- `include/atomiverse/cpu/gpw/gpw_solver.h`, `src/cpu/gpw/gpw_solver.cpp`
  (`overlap` / `overlap_k` extraction)
- `src/cpu/gpw/gpw_evaluator.cpp` (three CPMD branch restructures)
- `tests/cpu/gpw/test_gpw_cpmd.cpp` (overlap-parity test)

## Tasks

- [ ] RED: overlap-parity test — `solver.overlap(atoms)` equals the `last_.S`
      of a fixed-density `run()` at the same geometry element-wise
      (0 / ≤1e-15); multi-k twin against `last_.S_k`
- [ ] Extract `overlap` / `overlap_k` from the `run()` assembly loop
      (shared code path, not a copy — `run()` calls the same routine)
- [ ] Restructure the Γ RKS, Γ UKS, and multi-k CPMD branches to
      S→align→single-eval
- [ ] Re-run the full cpmd regression set — drift pins unchanged is the
      physics gate; note the step-time change in the close-out summary
- [ ] Scoped build + ctest green

## Testing strategy

**Unit:** overlap parity (Γ + multi-k); a KS-eval **call-count** assertion —
instrument via the existing `GpwResult` fixed-density payload path in the
test (one `run()` per production step observed through the solver's public
counters if present, else assert indirectly: the parity test + branch review
in the PR). **Regression:** no new case (performance is not a feature with a
literature value); the gates are the existing cpmd-water / cpmd-uks /
cpmd-kpoints / cpmd-rattle-water / cpmd-electronic-thermostat drift pins
passing unchanged.

## Out of scope

- Any change to BOMD / SGCP / XL-BOMD branches (they converge SCF — the
  double-eval pathology is CPMD-only)
- Algorithmic changes to align_to_overlap or the propagators
- New public counters (only if the call-count test genuinely needs one —
  then a minimal public member per the regression-witness rule, named in
  the PR)

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/gpw/gpw_evaluator.{h,cpp}`,
  `src/cpu/scf/mean_field_method.{h,cpp}`, `scf/internal/gradient.*`, or any
  CMakeLists, verify `git status` shows them clean. As of 2026-07-26 the
  operator's libcint-removal/C++20 work (~760 files) is uncommitted and
  includes them. If dirty: STOP and report.
- Bit-identical S is the design requirement: if the extracted routine and
  `run()` ever diverge in image counts or summation order, the parity test
  must catch it — do not loosen that tolerance.
- Land LAST in Wave 3: its branch surgery should sit on top of the
  esmm-cpmd-02/03 edits to the same functions, not under them.
