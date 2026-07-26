---
spec: esmm-cpmd-06-single-eval
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Standalone overlap is bit-identical to the run() overlap
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_gpw_cpmd.cpp shows GpwDft::overlap(atoms) equals the
      last_.S of a fixed-density run() at the same geometry element-wise
      (exact or <= 1e-15), and overlap_k matches last_.S_k per sampled k.
    status: pending
  - id: ac-002
    summary: CPMD production steps perform exactly one KS evaluation
    type: code
    pass_when: |
      All three CPMD branches (Gamma RKS, Gamma UKS, multi-k RKS) in
      gpw_evaluator.cpp follow S -> align_to_overlap -> single fixed-density
      eval -> forces -> advance; no branch calls the solver twice per step
      (verified by review + the parity test; a public counter only if added
      as a named minimal member).
    status: pending
  - id: ac-003
    summary: Physics pins unchanged
    type: runtime
    pass_when: |
      cpmd-water, cpmd-uks, cpmd-kpoints, cpmd-rattle-water,
      cpmd-electronic-thermostat pass with their existing drift/iters==0
      pins untouched (no tolerance edits).
    status: pending
  - id: ac-004
    summary: Extraction is shared, not copied
    type: code
    pass_when: |
      run()'s S assembly and overlap()/overlap_k() call the same underlying
      routine (single definition of the image-summed basis loop); no
      duplicated overlap loop exists in gpw_solver.cpp.
    status: pending
out_of_scope:
  - BOMD / SGCP / XL-BOMD branches
  - align_to_overlap or propagator changes
---

# Acceptance — `esmm-cpmd-06-single-eval`

**Done** means every GPW CPMD production step builds S standalone, aligns,
and runs exactly one fixed-density KS evaluation, with the standalone S
provably identical to the in-solve S and all shipping CPMD pins untouched.

## AC-001 — the identity that makes it safe

Bit-identical S is what turns "optimization" into "deletion of discarded
work".

## AC-002 — the actual win

One eval per step on all three branches — the ~2x on the dominant kernel.

## AC-003 — physics is the gate

Drift pins are the regression discipline for a performance change.

## AC-004 — no second truth

A copied loop would drift from run() eventually; sharing is the requirement.
