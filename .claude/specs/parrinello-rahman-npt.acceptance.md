---
spec: parrinello-rahman-npt
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Public PR NPT Integrator with full cell matrix
    type: code
    pass_when: |
      cpu::md::ParrinelloRahmanNPT exists as Integrator; couples full h (not
      diagonal-only); no public factory or Config bag.
    status: pending
  - id: ac-002
    summary: Unit tests for PR integrator
    type: runtime
    pass_when: |
      tests/cpu/dynamics/test_parrinello_rahman_npt passes under ctest.
    status: pending
  - id: ac-003
    summary: LJ PR-NPT regression gate
    type: runtime
    pass_when: |
      regressions/md/pr_npt_lj exits 0 with hard-coded physical window.
    status: pending
out_of_scope:
  - NPH
  - CUDA
---

# Acceptance — parrinello-rahman-npt

Done means CPU MD can run fully flexible-cell NPT via a public Integrator, with unit + regression gates.
