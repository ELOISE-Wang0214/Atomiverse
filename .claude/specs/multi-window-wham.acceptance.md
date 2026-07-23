---
spec: multi-window-wham
created: 2026-07-23
criteria:
  - id: ac-001
    summary: WHAM math lives in molrs with thin C++ API
    type: code
    pass_when: |
      WHAM self-consistency is implemented in molrs (not duplicated as full
      numeric core in cpu/bias); cpu exposes a thin public entry for rebuild.
    status: pending
  - id: ac-002
    summary: Unit tests for WHAM on synthetic data
    type: runtime
    pass_when: |
      tests/cpu/observe/test_wham passes.
    status: pending
  - id: ac-003
    summary: Multi-window WHAM regression
    type: runtime
    pass_when: |
      regressions/sampling/wham_double_well exits 0 with hard-coded F(s) gates.
    status: pending
out_of_scope:
  - MBAR
  - Multi-dimensional WHAM
  - CUDA
---

# Acceptance — multi-window-wham

Done means multi-window WHAM PMF reconstruction works end-to-end with math owned by molrs and a regressions/ gate.
