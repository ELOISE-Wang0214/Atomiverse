---
spec: ri-mp2
created: 2026-07-23
criteria:
  - id: ac-001
    summary: RI-MP2 energy on RHF public mp2 path
    type: code
    pass_when: |
      RHF::mp2 (or documented member) yields MP2Result using density fitting
      without materializing full AO ERI when RI enabled; no public IntegralPack.
    status: pending
  - id: ac-002
    summary: Unit RI vs canonical MP2
    type: runtime
    pass_when: |
      tests/cpu/scf/test_ri_mp2 passes with documented tolerance.
    status: pending
  - id: ac-003
    summary: H2O RI-MP2 regression
    type: runtime
    pass_when: |
      regressions/scf/h2o_ri_mp2 exits 0 with hard-coded e_corr gate.
    status: pending
out_of_scope:
  - MP2 gradients
  - UMP2
  - CUDA
---

# Acceptance — ri-mp2

Done means density-fitted closed-shell MP2 energy is available on the RHF API and regression-gated.
