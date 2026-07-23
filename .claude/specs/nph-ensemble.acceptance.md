---
spec: nph-ensemble
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Thermostat-free NPH Integrator on CPU MD
    type: code
    pass_when: |
      Public Integrator runs NPH without fused atom thermostat; no factory API.
    status: pending
  - id: ac-002
    summary: Unit tests for NPH integrator
    type: runtime
    pass_when: |
      tests/cpu/dynamics/test_nph_integrator passes.
    status: pending
  - id: ac-003
    summary: LJ NPH regression enthalpy gate
    type: runtime
    pass_when: |
      regressions/md/nph_lj exits 0 with hard-coded H drift bound.
    status: pending
out_of_scope:
  - CUDA
  - Full anisotropic NPH if isotropic v1 ships
---

# Acceptance — nph-ensemble

Done means NPH is reachable on CPU MD with a conserved-enthalpy regression gate.
