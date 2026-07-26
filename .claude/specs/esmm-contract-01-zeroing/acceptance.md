---
spec: esmm-contract-01-zeroing
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Hessian zeroes per displacement — accumulate stub gives H00 == k
    type: runtime
    pass_when: |
      tests/cpu/vibration/test_hessian_zeroing.cpp passes: a `+=`-style
      harmonic stub under vib::Hessian yields H00 == k (relative error
      ≤ 1e-8), and both Hessian ctors (full chain + WithoutAccumulator)
      carry ZeroForces@-1000 / ZeroEnergy@-999 on vib::Phase::force.
    status: pending
  - id: ac-002
    summary: PathDriver zeroes per replica pass
    type: runtime
    pass_when: |
      tests/cpu/path/test_path_zeroing.cpp passes: two PathDriver iterations
      with a `+=` stub produce identical replica forces at identical
      geometry (no doubling).
    status: pending
  - id: ac-003
    summary: MC zeroes forces; Minimize/Evaluate zero stress
    type: code
    pass_when: |
      mc_base.h binds ZeroForces@sampling::Phase::energy order -1000;
      minimize.h and evaluate.h bind ZeroStress order -998 on their eval
      phases; updated plan-shape anchors in test_driver_binding_kind.cpp
      (Hessian/Mc/Path) pass.
    status: pending
  - id: ac-004
    summary: Fossilized k/2 oracle deleted and tightened
    type: runtime
    pass_when: |
      reg_harmonic_driver_matrix no longer contains the "k or k/2" hedge and
      asserts the exact analytic scale (H00 == k within its stated tol);
      the regression passes.
    status: pending
  - id: ac-005
    summary: Behavioral audit — nothing else moves
    type: runtime
    pass_when: |
      Full ctest suite green; ES-based Hessian regressions
      (regressions/scf/*hessian*, h2o_rhf_freq, h2o_rhf_ir), NEB regression,
      and test_cell_relax pass unchanged (only harmonic_driver_matrix
      numbers changed).
    status: pending
out_of_scope:
  - Evaluator `=`→`+=` switches (esmm-contract-05)
  - run_slot_ nested dispatch (esmm-contract-03)
  - CUDA
---

# Acceptance — `esmm-contract-01-zeroing`

**Done** means every CPU driver that hosts force/energy/stress sources zeroes
the accumulator channels itself (MD, Evaluate, Minimize, MC, Hessian, Phonon,
PathDriver), the k/2 fossil oracle is tightened to the true analytic value,
and nothing else in the suite moves.

## AC-001 / AC-002 — the two live accumulation bugs die

Hessian and PathDriver are the drivers where `+=` sources are wrong today.

## AC-003 — uniformity

MC force growth and Minimize/Evaluate stress accumulation closed; anchors
updated so the plan shape is pinned, not accidental.

## AC-004 — no fossilized hedges

A regression oracle that accepts a bug's output is worse than none.

## AC-005 — the switch precondition

esmm-contract-05 may only land on a tree where zero-then-accumulate is
equivalent to overwrite for every sole-source configuration.
