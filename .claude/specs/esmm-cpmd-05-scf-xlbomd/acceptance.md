---
spec: esmm-cpmd-05-scf-xlbomd
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Guess ladder resolves in the specified order
    type: runtime
    pass_when: |
      tests/cpu/scf/test_scf_xl_bomd.cpp proves resolution order
      explicit > xl aux() > ASPC > last-P (stubbed sequence probes which
      source seeded the SCF).
    status: pending
  - id: ac-002
    summary: Exclusivity gates throw and are guard-tested
    type: runtime
    pass_when: |
      xl_bomd with predictor() throws; xl_bomd with cpmd throws (both
      directions); each has an EXPECT_THROW in test_scf_xl_bomd.cpp or
      test_scf_cpmd.cpp.
    status: pending
  - id: ac-003
    summary: Aux density advances only on converged nested runs
    type: runtime
    pass_when: |
      Unit test shows xl aux() history is untouched by a standalone run()
      and by an unconverged nested run, and advances once per converged
      nested run (restricted: factor-2 total; UHF/UKS: per-channel).
    status: pending
  - id: ac-004
    summary: xl-bomd-rhf-water regression pins the XL signature
    type: runtime
    pass_when: |
      regressions/dynamics/xl-bomd-rhf-water passes: warm-step mean SCF
      iterations strictly below the cold-start twin by the margin frozen in
      REFERENCE.md; energy-drift envelope holds; REFERENCE.md cites
      Niklasson JCP 130, 214109 (2009); <= 60 s, public API only.
    status: pending
out_of_scope:
  - Linearized / zero-SCF XL (GPW-only)
  - Multi-k, CUDA, GPW wiring changes
---

# Acceptance — `esmm-cpmd-05-scf-xlbomd`

**Done** means any Gaussian method warm-starts MD through a propagated
auxiliary density with the documented guess-ladder precedence, exclusive
with ASPC and CPMD, witnessed by an iteration-reduction regression.

## AC-001 — one ladder, one order

Guess precedence is the entire contract; it must be provable, not implied.

## AC-002 — one guess owner at a time

Two propagation schemes feeding one SCF is the bug class the gates prevent.

## AC-003 — history integrity

Rejected/standalone evaluations must not corrupt the time-continuous aux
history.

## AC-004 — the XL signature

Fewer warm-step iterations at stable drift is what XL-BOMD *is*.
