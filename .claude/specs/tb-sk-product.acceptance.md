---
spec: tb-sk-product
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Multi-element static SK tables for DFTB2
    type: code
    pass_when: |
      DFTB2 supports at least H/C/N/O via in-tree static tables; unsupported
      pairs throw; no runtime .skf parser; no public factory.
    status: pending
  - id: ac-002
    summary: Unit tests for multi-element DFTB2
    type: runtime
    pass_when: |
      tests/cpu/semiempirical/test_dftb2_sk passes.
    status: pending
  - id: ac-003
    summary: H2O DFTB2 mio regression
    type: runtime
    pass_when: |
      regressions/semiempirical/h2o_dftb2_mio exits 0 with hard-coded gates.
    status: pending
out_of_scope:
  - DFTB3
  - GFN2 full tables
  - CUDA
---

# Acceptance — tb-sk-product

Done means DFTB2 is a multi-element static-table method with forces and a water regression, not H-only toy.
