---
slug: esmm-cpmd-05-scf-xlbomd
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: electronic-propagation / cpmd
queue: esmm-campaign
deferred: false
depends_on: [esmm-cpmd-04-ao-cpmd]
unblocks: []
grilled: true
---

# Spec: `esmm-cpmd-05-scf-xlbomd` — XL-BOMD on the Gaussian SCF methods

**Slug:** `esmm-cpmd-05-scf-xlbomd`
**Domain:** Gaussian-basis XL-BOMD (CPU)
**Status:** approved (grilled via plan-mode design review + operator Q&A, 2026-07-26)
**Ordering note:** `depends_on: esmm-cpmd-04` is churn-ordering only (both
rewrite `mean_field_method.*`); the logical dependency is esmm-cpmd-01
(`electronic::XlBomd`).

## Summary

GPW has extended-Lagrangian BOMD (`XlBomd` auxiliary-density propagation);
the Gaussian methods have only ASPC. With the propagators now shared in
`cpu/electronic/` (esmm-cpmd-01), porting XL-BOMD to `MeanFieldMethod` is
~100 lines: two members, a guess-ladder hook, one advance call, and the
exclusivity gates. `XlBomd` is unit-agnostic dense algebra on AO densities —
no new physics. Directly serves the operator's method×propagator
free-combination directive at the same seam CPMD opened.
Niklasson, *J. Chem. Phys.* **130**, 214109 (2009).

## Design

### Public members on `MeanFieldMethod<Derived>`

```cpp
/// @brief Extended-Lagrangian BOMD (Niklasson 2009): propagate an auxiliary
///        density as the SCF guess so warm steps converge in a few
///        iterations without ASPC's corrector bias. Nested-MD use; exclusive
///        with predictor() and cpmd (typed throws). Default false.
bool xl_bomd = false;
/// @brief XL-BOMD dissipation order K (history depth). Default 5.
int xl_bomd_k = 5;
```

Privates: `std::optional<electronic::XlBomd> xl_alpha_, xl_beta_;`.

### Guess ladder in `drive_scf_` (single ordered resolution)

`explicit guess` > `xl aux()` > `ASPC predict` > `last-P reuse`. After each
**converged** nested run: `xl_alpha_->advance(P_alpha_)` (+β for UHF/UKS;
restricted advances the factor-2 total density, matching the loop
convention). Lazy construction on first nested run (needs `nbasis`).

### Gates (typed throws + guard tests, same change)

`xl_bomd && use_predictor_` → throw (two guess owners); `xl_bomd && cpmd` →
throw (already half-present from esmm-cpmd-04's lattice — this spec makes it
symmetric); standalone `run()`/`forces()` with `xl_bomd` allowed **only** as
plain BOMD-style calls? No — mirror GPW: XL state advances on nested runs
only; standalone runs leave the aux history untouched (document). ROHF/GHF:
allowed (XL is a guess strategy, not a propagator of the SCF state) — but β
channel only exists for UHF/UKS; ROKS/ROHF advance the total density.

## Files to create or modify

- `include/atomiverse/cpu/scf/mean_field_method.h`, `src/cpu/scf/mean_field_method.cpp`
- NEW `tests/cpu/scf/test_scf_xl_bomd.cpp`
- NEW `regressions/dynamics/xl-bomd-rhf-water/` (`reg_*.cpp` + `REFERENCE.md`)
- `regressions/CMakeLists.txt`, `tests/CMakeLists.txt`
- NEW `docs/tutorial/methods/xl-bomd-rhf-water.md` (+ index row)

## Tasks

- [ ] RED: `test_scf_xl_bomd.cpp` — guess-ladder order (explicit beats XL
      beats ASPC beats last-P, probed via a stub sequence); exclusivity
      throws (xl×predictor, xl×cpmd); aux advances only after converged
      nested runs
- [ ] Implement members + lazy `electronic::XlBomd` + ladder hook + advance
- [ ] Regression `xl-bomd-rhf-water`: RHF/STO-3G water NVE MD; pins
      (a) warm-step SCF-iteration reduction vs a cold-start twin run,
      (b) energy-drift envelope, (c) Niklasson JCP 130, 214109 citation in
      REFERENCE.md
- [ ] Tutorial page + roadmap row; scoped build + ctest + regression green

## Testing strategy

**Unit:** ladder order, exclusivity, advance-on-converged-only, restricted
vs unrestricted advance conventions. **Regression:** `xl-bomd-rhf-water`
(≤60 s, `ATV_BUILD_HF AND ATV_BUILD_MD`, public API) — iteration-count
reduction is the XL signature; drift envelope guards stability.

## Out of scope

- Linearized / zero-SCF shadow XL-BOMD (needs the shadow-energy force
  assembly; stays GPW-only)
- Multi-k anything; CUDA
- Changes to GPW's XL wiring

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/gpw/gpw_evaluator.{h,cpp}`,
  `src/cpu/scf/mean_field_method.{h,cpp}`, `scf/internal/gradient.*`, or any
  CMakeLists, verify `git status` shows them clean. As of 2026-07-26 the
  operator's libcint-removal/C++20 work (~760 files) is uncommitted and
  includes them. If dirty: STOP and report.
- The ladder must be one readable resolution block — no scattered ifs.
- Guard tests for the two new exclusivity gates land in the same change
  (house rule).
