---
spec: full-grid-alda
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Full MO-product ALDA fills RKS TDA A-matrix
    type: code
    pass_when: |
      RKS LSDA TDA path builds (ia|f_xc|jb) via grid MO products, not only
      diagonal Dirac shift; ExcitedStatesResult.alda is true for LSDA.
    status: pending
  - id: ac-002
    summary: Unit tests pin symmetry and diagonal limit
    type: runtime
    pass_when: |
      tests/cpu/scf/test_full_grid_alda and updated test_alda_diff pass under ctest.
    status: pending
  - id: ac-003
    summary: H2 full-grid ALDA regression gates hard-coded omega
    type: runtime
    pass_when: |
      regressions/scf/h2_full_grid_alda embeds hard-coded S1 [eV] + tolerance,
      exits 0, no live third-party oracle at run time.
    status: pending
out_of_scope:
  - Hybrid TDDFT
  - Analytic ES forces / NAC vectors
  - CUDA
---

# Acceptance — full-grid-alda

Done means LSDA TDA on CPU uses a real-space full-grid ALDA kernel on the A matrix, unit-tested and regression-gated on H2, without claiming hybrid or force completeness.
