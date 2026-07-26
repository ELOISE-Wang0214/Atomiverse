---
spec: ot-dftu
created: 2026-07-26
criteria:
  - id: ac-001
    summary: UKS OT + plus_u no longer throws
    type: code
    pass_when: |
      src/cpu/gpw/gpw_solver.cpp UKS OT path does not throw on do_dftu;
      gpw_solver.h / gpw_evaluator docs no longer claim DFT+U+OT throws at Γ.
    status: met
  - id: ac-002
    summary: Unit OT+PLUS_U converges with non-trivial E_U
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_ot.cpp GpwOt.UksOtPlusUConverges: small UKS
      + plus_u, both OT and DIIS converge, e_dftu > 0, energy class match.
    status: met
  - id: ac-003
    summary: Regression vs offline CP2K OT+PLUS_U pin
    type: runtime
    pass_when: |
      regressions/gpw/ot-dftu builds as reg_ot-dftu / run_reg_ot-dftu;
      |E_atv - E_cp2k| ≤ tol in REFERENCE.md; never launches CP2K; ≤ 60 s.
    status: met
  - id: ac-004
    summary: REFERENCE documents pin + perf
    type: docs
    pass_when: |
      REFERENCE.md cites CP2K version, input path, E_tot, optional wall ratio
      or perf_exempt reason.
    status: met
  - id: ac-005
    summary: API iron laws
    type: code
    pass_when: |
      No public *Config / make_* / build_* for this feature; use_ot + plus_u
      members only; public Frame energy still eV.
    status: met
out_of_scope:
  - OT + k-points
  - Live oracle in ctest
  - CUDA
---

# Acceptance — `ot-dftu`

**Done** means Γ-UKS OT with Dudarev PLUS_U is a supported product path,
pinned offline to CP2K, with unit OT/DIIS parity and no architecture debt.
