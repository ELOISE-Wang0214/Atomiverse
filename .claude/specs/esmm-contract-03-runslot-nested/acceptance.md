---
spec: esmm-contract-03-runslot-nested
created: 2026-07-26
criteria:
  - id: ac-001
    summary: run_stage replays nested driver nodes
    type: runtime
    pass_when: |
      tests/cpu/driver/test_driver_nesting.cpp "RunStageReplaysNestedDriver"
      passes: after backend run_stage(force), the nested node's call counter
      incremented (was 0 before the fix).
    status: pending
  - id: ac-002
    summary: Line-search trials re-enter the nested force source
    type: runtime
    pass_when: |
      Minimize + StrongWolfe unit test passes: nested stub call count >
      outer iteration count after run() (trial evaluations re-entered it).
    status: pending
  - id: ac-003
    summary: Re-entrancy and loop isolation preserved
    type: runtime
    pass_when: |
      Existing test_driver_nesting cases still pass (inner loop cannot stop
      the outer); no resolution-window throw appears in the new paths; full
      ctest green.
    status: pending
  - id: ac-004
    summary: R2 regression pins a working QM line search
    type: runtime
    pass_when: |
      regressions/minimizer/rhf_h2_linesearch builds under ATV_BUILD_HF,
      runs ≤ 60 s, asserts |r_e − 0.7122 Å| ≤ 5e-3 Å (CCCBDB RHF/STO-3G)
      and iteration count ≤ the REFERENCE.md bound; public headers only.
    status: pending
  - id: ac-005
    summary: Timing deferral is a recorded decision
    type: docs
    pass_when: |
      run_slot_ and step_plan carry comments stating nested-op timing is
      deliberately deferred (not an oversight), citing this spec's slug.
    status: pending
out_of_scope:
  - Nested-op timing implementation
  - forces_requested plumbing (esmm-contract-04)
---

# Acceptance — `esmm-contract-03-runslot-nested`

**Done** means every `run_stage`-based replay (Minimize line search,
PathDriver per-replica force pass, any future caller) executes nested method
drivers exactly like the main step loop does, pinned by a unit that counts
re-entries and a regression whose convergence is only reachable with real
trial energies.

## AC-001 / AC-002 — the E=0/f=0 bug dies

Trial evaluations must carry the nested SCF's real energy and forces.

## AC-003 — nothing loosens

`ctx.loop` isolation and resolution-window discipline survive the new
dispatch path.

## AC-004 — literature-pinned witness

CCCBDB r_e plus an iteration bound that a degenerate line search cannot meet.

## AC-005 — honest deferral

What is not done is written down where the next reader will look.
