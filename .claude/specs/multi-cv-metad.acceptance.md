---
spec: multi-cv-metad
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Multi-CV well-tempered metadynamics Bias
    type: code
    pass_when: |
      Public Bias accepts ≥2 CollectiveVariables and deposits multidimensional
      hills with WT heights; no factory.
    status: pending
  - id: ac-002
    summary: Unit tests for ND metadynamics
    type: runtime
    pass_when: |
      tests/cpu/bias/test_metadynamics_nd passes.
    status: pending
  - id: ac-003
    summary: 2D metadynamics regression
    type: runtime
    pass_when: |
      regressions/sampling/metad_2d exits 0 with hard-coded gates.
    status: pending
out_of_scope:
  - OPES
  - WHAM
  - CUDA
---

# Acceptance — multi-cv-metad

Done means well-tempered metadynamics works for more than one CV on CPU with a regression gate.
