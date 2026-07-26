---
slug: esmm-contract-02-validation-walks
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: driver / pipeline-contract
queue: esmm-campaign
deferred: false
depends_on: []
unblocks: [rhf-in-mc, rhf-in-neb, esmm-contract-04-forces-requested]
grilled: true
---

# Spec: `esmm-contract-02-validation-walks` — MC and PathDriver see nested method drivers

**Slug:** `esmm-contract-02-validation-walks`
**Domain:** CPU driver validation (energy-source recognition)
**Queue:** esmm-campaign (Wave 1 — composition contract)
**Status:** approved. Grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

Energy-source acceptance is a `dynamic_cast` walk against the kind set
{Potential, Bias, PathCoupling, HF, Semiempirical, PlaneWave, GPW}.
MD/Evaluate/Minimize/Hessian use the **node-aware** walk
`for_each_pending_kind(comp_fn, driver_fn)` (`driver.h:317-324`), whose
driver-node lambda recognizes a nested `MeanFieldMethod` through its `HF` tag
(`md.h:374-382`). Two drivers still use the **component-only**
`for_each_pending`:

- `McBase::validate_components` (`mc_base.h:160-181`) — a nested
  `cpu::scf::RHF` bound at `sampling::Phase::energy` **fails MC validation**
  even though it runs fine. RHF-in-MC is a row of the free-combination matrix.
- `PathDriver::validate_components` (`path_driver.h:158-178`) — worse: with
  `NebCoupling` present (a PathCoupling kind), a nested-RHF-only band *passes*
  validation for the wrong reason, and (until `esmm-contract-03`) `run_slot_`
  then silently skips the RHF — zero forces, instant garbage "convergence".
  The walk fix makes the nested HF node *count* as the force source it is.

This spec switches both drivers to the node-aware walk. It is header-only plus
tests, and it is a prerequisite for `esmm-contract-04` (whose R3 regression
binds a nested RHF inside MC).

## Domain basis

Kind recognition is by `dynamic_cast` to kind classes — never name matching,
never a free-function classifier (component-model hard cut,
`.claude/notes/notes.md:438-481`). `HF` is a pure polymorphic tag
(`component/hf.h:13-16`), not a Component; concrete Gaussian methods are
`Driver<scf::Phase> + HF` and bind as nested driver nodes. A validation walk
that skips driver nodes therefore cannot see any Gaussian method.

Rejected alternative: making `HF` a Component so the component walk sees it —
would force a launch() on a Driver-shaped method and contradict the shipped
node-erasure design (`driver_node.h:4-7`).

## Design

Both call sites switch to the exact MD pattern:

```cpp
for_each_pending_kind(
    [&](const cpu::Component& c) { /* existing kind casts, unchanged */ },
    [&](const cpu::DriverNode* d) {
        if (dynamic_cast<const HF*>(d)) has_energy = true;   // has_force in PathDriver
    });
```

`for_each_pending_kind` lives on `cpu::Driver` (`driver.h:317-324`), and both
`McBase` (`Driver<sampling::Phase>`) and `PathDriver` (`Driver<path::Phase>`)
derive it — verified, no plumbing needed.

## Files to create or modify

- `include/atomiverse/cpu/montecarlo/mc_base.h` — `validate_components`
  (`:160-181`): `for_each_pending` → `for_each_pending_kind` + HF-node lambda
- `include/atomiverse/cpu/path/path_driver.h` — `validate_components`
  (`:158-178`): same transformation for `has_force`
- `tests/cpu/test_mc_scaffold.cpp` — new case: nested HF node satisfies the MC
  energy requirement (stub `Driver<scf::Phase> + HF`; check
  `tests/cpu/driver_test_helpers.h` for an existing stub shape first)
- `tests/cpu/path/` — new/extended case: nested HF node satisfies the
  PathDriver force requirement

All files are outside the operator's in-flight dirty set (verified
2026-07-26).

## Tasks

- [ ] RED: MC scaffold test — `MC` with only a nested stub `Driver + HF` at
  `sampling::Phase::energy` must `build()` (fails today with the
  missing-energy-source throw)
- [ ] RED: PathDriver test — band with `NebCoupling` removed and only a nested
  stub HF node must `build()` (fails today)
- [ ] Switch both validate_components to `for_each_pending_kind`
- [ ] GREEN + full ctest

## Testing strategy

**Unit (U3):** the two RED-first build()-passes above, plus a negative case
each (no source at all still throws — the walk change must not weaken the
guard).

**Regression:** none — validation-only change; the e2e witness arrives with
`esmm-contract-04`'s R3 (nested RHF actually sampled under MC).

## Out of scope

- run_slot_ dispatch of nested drivers (esmm-contract-03 — without it a
  nested-HF-only PathDriver band still computes zero forces; this spec only
  fixes *recognition*)
- forces_requested plumbing (esmm-contract-04)
- Any evaluator accumulation change

## Notes for implementers

- Keep the existing component-side kind list byte-identical in both drivers —
  only the walk shape and the added node lambda change.
- MD's `langevin` name-escape (`md.h`) is MD-specific; do not copy it into MC
  or Path.
