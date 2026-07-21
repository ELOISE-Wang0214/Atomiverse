---
spec: gaff-md-path-06-protocol-smoke
created: 2026-07-21
criteria:
  - id: ac-001
    summary: Smoke regression completes without NaN
    type: runtime
    pass_when: |
      reg_gaff_protocol_smoke exits 0 within ≤120 s on a single rank Release
      build; no NaN/Inf in final printed thermo fields (temp, pe, volume/density).
    status: pending
  - id: ac-002
    summary: Protocol includes minimize and at least one NVT and one NPT segment
    type: code
    pass_when: |
      Source uses cpu::minimizer::Minimize (or equivalent) and cpu::md::MD with
      both a thermostat-only stage and NoseHooverNPT (or documented barostat)
      stage; REFERENCE lists step counts.
    status: pending
  - id: ac-003
    summary: Wall-time budget documented and ctest-registered
    type: runtime
    pass_when: |
      Target is registered under regressions and completes under the budget in
      REFERENCE.md on the project default CI-like CPU config.
    status: pending
out_of_scope:
  - Paper-length 30 ns runs
  - Rg / density plateau matching flexible/ figures
  - Multi-replica statistics
---

# Acceptance — gaff-md-path-06-protocol-smoke

Done means the GAFF stack survives a short staged min/NVT/NPT smoke on a
periodic charged polymer fragment.
