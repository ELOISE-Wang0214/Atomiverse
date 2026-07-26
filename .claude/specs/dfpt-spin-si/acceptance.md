---
spec: dfpt-spin-si
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Collinear spin accepted; smear + hybrid still rejected
    type: runtime
    pass_when: |
      test_dfpt: spin_polarized=true no longer throws at the spin gate for
      pure LSDA/GGA + fixed occupations; smearing + spin still throws;
      hybrids still throw.
    status: met
  - id: ac-002
    summary: nspin=2 Sternheimer response in dfpt_gamma_force_constants
    type: code
    pass_when: |
      dfpt_gamma_force_constants runs collinear fixed-occupation SCF and
      spin-resolved Sternheimer (Δn↑/Δn↓, spin XC kernel); no
      "collinear spin unsupported" throw on the Γ path.
    status: met
  - id: ac-003
    summary: reg_dfpt-spin-si matches offline energy-FD pin
    type: runtime
    pass_when: |
      regressions/response/dfpt-spin-si builds with ATV_BUILD_PW, never
      launches CP2K, and |C00_dfpt − C00_energy_FD| ≤ tol in REFERENCE.md;
      optical / acoustic / ASR gates pass; wall ≤ 60 s.
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
      T_atv/T_oracle ≤ 10 when CP2K wall recorded; else perf_exempt with
      reason in REFERENCE.md (no CP2K module on host).
    status: met
---
