---
spec: esmm-cpmd-02-timestep-protocol
created: 2026-07-26
criteria:
  - id: ac-001
    summary: cpmd_dt is gone from the tree
    type: code
    pass_when: |
      grep -rn "cpmd_dt" over include/ src/ tests/ regressions/ docs/ returns
      zero hits; GPWEvaluator has no public dt member.
    status: pending
  - id: ac-002
    summary: Conflicting or missing providers throw at MD build
    type: runtime
    pass_when: |
      tests/cpu/dynamics/test_timestep_protocol.cpp shows (a) two providers
      with dt differing >1e-12 relative throw naming both components,
      (b) a timestep_required consumer with zero providers throws,
      (c) a stub consumer receives the pushed dt before its build().
    status: pending
  - id: ac-003
    summary: CPMD outside MD throws at component build
    type: runtime
    pass_when: |
      A cpmd-enabled GPWEvaluator hosted by Evaluate (or any driver without a
      TimestepProvider) throws at build with the "CPMD outside MD is
      meaningless" message; guard test in test_gpw_cpmd.cpp.
    status: pending
  - id: ac-004
    summary: cpmd && compute_stress throws instead of silently ignoring
    type: runtime
    pass_when: |
      GPWEvaluator::build throws a typed exception when cpmd and
      compute_stress are both set; EXPECT_THROW guard added in
      test_gpw_cpmd.cpp (replacing today's silent drop).
    status: pending
  - id: ac-005
    summary: The 5 cpmd regressions pass on a single dt truth
    type: runtime
    pass_when: |
      cpmd-water, cpmd-electronic-thermostat, cpmd-kpoints, cpmd-uks,
      cpmd-rattle-water pass with dt set ONLY on VelocityVerlet; their
      physics pins (drift, iters==0) are unchanged.
    status: pending
  - id: ac-006
    summary: Guard tests updated in the same change
    type: code
    pass_when: |
      The stale negative-cpmd_dt EXPECT_THROW is deleted; set_timestep
      validation (non-finite/<=0 throws) is added; no orphaned guards
      reference the deleted member.
    status: pending
out_of_scope:
  - MeanFieldMethod consumer wiring (esmm-cpmd-04)
  - Thermostat dt cross-validation (separate operator decision)
  - Variable-cell CP physics
---

# Acceptance — `esmm-cpmd-02-timestep-protocol`

**Done** means the nuclear timestep has exactly one owner (the integrator),
CPMD consumers receive it through `TimestepProvider`/`TimestepConsumer` at MD
build, `cpmd_dt` no longer exists anywhere, and the two silent failure modes
(mismatched dt, ignored stress request) are build-time throws.

## AC-001 — duplicated truth deleted

The protocol exists to kill `cpmd_dt`; keeping it validated would preserve
the failure mode.

## AC-002 — push semantics

Conflict and absence are loud, named, and pre-build.

## AC-003 — CPMD is an MD composition

Evaluate/Minimize hosting a CPMD evaluator is meaningless and must say so.

## AC-004 — no silent capability drops

Gate-lattice honesty: a user asking for stress must not silently lose it.

## AC-005 — physics unchanged

The 5 shipping cpmd regressions are the no-change witnesses.

## AC-006 — guard hygiene

House rule: gates and their EXPECT_THROW twins move together.
