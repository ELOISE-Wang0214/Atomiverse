spec: gpw-d3bj-stress
created: 2026-07-31
criteria:
  - id: ac-001
    summary: Offline CP2K D3(BJ)+plain NUMERICAL stress pin recorded
    type: docs
    pass_when: |
      tools/cp2k_oracle/inputs/gpw-d3bj-stress/ holds both .inp twins;
      regressions/gpw/gpw-d3bj-stress/REFERENCE.md records CP2K 2026.1
      version string, host, wall, the printed dispersion energy, both
      numerical stress tensors, and their difference with units.
    status: verified (agent-auto 2026-07-31)
  - id: ac-002
    summary: E_disp matches CP2K printed dispersion energy
    type: runtime
    pass_when: |
      reg_gpw-d3bj-stress gates |E_disp(ATV) - E_disp(CP2K)| <= 1e-8 Ha on
      the twin geometry.
    status: verified (agent-auto 2026-07-31)
  - id: ac-003
    summary: D3(BJ) dispersion stress matches CP2K stress difference
    type: runtime
    pass_when: |
      Diagonal dispersion stress from the geometry-only D3(BJ) strain-FD in
      the regression matches sigma(D3BJ)-sigma(plain) from the two CP2K
      NUMERICAL stress runs within the REFERENCE.md tolerance.
    status: verified (agent-auto 2026-07-31)
  - id: ac-004
    summary: Dual gate analytic vs multi-SCF strain-FD with d3bj=true
    type: runtime
    pass_when: |
      Product analytic stress (GPWEvaluator d3bj=true, stress_use_fd=false)
      matches gpw_stress_fd on the same method within the REFERENCE.md
      tolerance (<= 2e-5 Ha/Bohr^3), off-diagonals ~0, diag above floor,
      kStress eV/A^3 conversion exact.
    status: verified (agent-auto 2026-07-31)
  - id: ac-005
    summary: Unit-tier D3BJ stress fold gate on cheap H2
    type: runtime
    pass_when: |
      GpwStress.AnalyticMatchesFdWithD3Bj passes in test_cpu_gpw_stress
      (cheap MT H2, no oracle, seconds tier).
    status: verified (agent-auto 2026-07-31)
  - id: ac-006
    summary: Perf within 10x of oracle or exempt with stated reason
    type: performance
    pass_when: |
      REFERENCE.md states T_atv and T_cp2k for comparable work and the ratio
      <= 10, or perf_exempt with the reason recorded.
    status: verified (agent-auto 2026-07-31)
