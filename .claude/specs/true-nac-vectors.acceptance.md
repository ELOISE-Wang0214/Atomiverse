---
spec: true-nac-vectors
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Public TDA NAC vectors 3N per atom pair of roots
    type: code
    pass_when: |
      RHF provides tda_nac_vectors (or equivalent) returning 3N components [1/Å]
      between two TDA roots; no public Config bag.
    status: pending
  - id: ac-002
    summary: Unit tests for antisymmetry and size
    type: runtime
    pass_when: |
      tests/cpu/scf/test_nac_vectors passes.
    status: pending
  - id: ac-003
    summary: H2 NAC vectors regression
    type: runtime
    pass_when: |
      regressions/scf/h2_nac_vectors exits 0 with hard-coded gates.
    status: pending
out_of_scope:
  - Surface hopping MD
  - TDDFT NAC
  - CUDA
---

# Acceptance — true-nac-vectors

Done means true first-order NAC vectors exist on the RHF TDA public surface and are regression-gated on a textbook diatomic.
