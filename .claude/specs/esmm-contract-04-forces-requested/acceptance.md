---
spec: esmm-contract-04-forces-requested
created: 2026-07-26
criteria:
  - id: ac-001
    summary: MC energy phase launches with forces_requested == false
    type: runtime
    pass_when: |
      U4 unit test: a stub component at sampling::Phase::energy observes
      ctx.forces_requested == false; the same stub at md::Phase::force,
      minimizer::Phase::force, and evaluate::Phase::evaluate observes true.
    status: pending
  - id: ac-002
    summary: Energy-only veto propagates through nesting
    type: runtime
    pass_when: |
      U4 nested case: a stub nested driver bound inside the MC energy phase
      launches ITS OWN slots with forces_requested == false (ExecServices
      AND channel), and restores true afterwards for outer phases.
    status: pending
  - id: ac-003
    summary: ASPC predictor under MC throws typed
    type: runtime
    pass_when: |
      U6 unit test: MeanFieldMethod with predictor() reached via
      run_nested_erased with ctx.mc.proposal non-null throws
      std::runtime_error naming predictor() and sampling; plain last-density
      warm start (no predictor) does not throw.
    status: pending
  - id: ac-004
    summary: R3 pins QM-in-MC against the analytic variance
    type: runtime
    pass_when: |
      regressions/mc/mc_bias_es builds under ATV_BUILD_HF, runs ≤ 60 s,
      nested RHF H2 + HarmonicBias umbrella; measured ⟨(r−r0)²⟩ matches the
      REFERENCE.md analytic k_BT/k_eff within its stated tolerance; no SCF
      gradient work is triggered (energy-only path).
    status: pending
  - id: ac-005
    summary: No behavior change outside MC
    type: runtime
    pass_when: |
      Full ctest green; MD/Minimize/Evaluate/Hessian regressions unchanged
      (their slots all stamp true).
    status: pending
out_of_scope:
  - Evaluator accumulation switch (esmm-contract-05)
  - Rollback-notification hooks
---

# Acceptance — `esmm-contract-04-forces-requested`

**Done** means the energy-only policy is a compiled Plan property that reaches
every source (including arbitrarily nested method drivers), MC trials with QM
engines skip gradients, and the one genuinely unsound MC warm-start path (the
ASPC predictor) is a loud typed error instead of a silent detailed-balance
violation.

## AC-001 / AC-002 — plumbing

Slot-stamped policy, Context-carried veto — components keep reading the one
existing flag.

## AC-003 — sampling stays unbiased

Last-density reuse is guess-only (correctness-neutral); extrapolated history
is not. The line between them is enforced, not documented.

## AC-004 — first QM-in-MC witness

Analytic pin (stiff-spring variance), fast because trials are energy-only.

## AC-005 — everyone else unchanged

The default is true everywhere a force consumer exists.
