---
spec: dfpt-gga-si
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Pure GGA accepted; hybrids still rejected
    type: runtime
    pass_when: |
      test_dfpt: PBE0 (hybrid) throws invalid_argument; empty UPF with PBE
      still throws; LSDA path still accepted.
    status: met
  - id: ac-002
    summary: GGA kernel used in Sternheimer response (not LDA-only gate)
    type: code
    pass_when: |
      dfpt_gamma_force_constants / dfpt_q_force_constants no longer hard-require
      XCType::LSDA; pure GGA uses XCBuilder FD kernel; pure LDA keeps local fxc.
    status: met
  - id: ac-003
    summary: reg_dfpt-gga-si matches offline energy-FD pin
    type: runtime
    pass_when: |
      regressions/response/dfpt-gga-si builds with ATV_BUILD_PW, never launches
      CP2K, and |C00_dfpt − C00_energy_FD| ≤ tol in REFERENCE.md; optical /
      acoustic / ASR gates pass; wall ≤ 60 s.
    status: met
  - id: ac-004
    summary: No public Config / factories
    type: code
    pass_when: |
      API remains DfptPhonon members + compute(PW, Frame); no *Config bag.
    status: met
  - id: ac-005
    summary: perf gate
    type: runtime
    pass_when: |
      T_atv/T_oracle ≤ 10 when CP2K wall recorded; else perf_exempt with reason
      in REFERENCE.md (no CP2K module on host).
    status: met
---
