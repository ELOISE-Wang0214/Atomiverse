---
spec: gpw-pbe-d3bj-ice
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Offline CP2K DFTD3(BJ) pin recorded
    type: docs
    pass_when: |
      regressions/gpw/gpw-pbe-d3bj-ice/REFERENCE.md cites CP2K 2026.1,
      geometry, TYPE DFTD3(BJ) / PBE params, E_disp = -0.00137627769481 Ha,
      and the twin input under tools/cp2k_oracle/inputs/gpw-pbe-d3bj-ice/.
    status: passed
  - id: ac-002
    summary: Public GPW d3bj API without Config bags
    type: code
    pass_when: |
      GpwDft and GPWEvaluator expose bool d3bj (and GpwResult::e_disp);
      no new public *Config; no make_*/build_* factories; d3 and d3bj
      mutually exclusive (throw if both true).
    status: passed
  - id: ac-003
    summary: Periodic water-dimer E_disp matches CP2K pin
    type: runtime
    pass_when: |
      Unit and/or regression assert |E_disp - (-0.00137627769481)| < 1e-8 Ha
      for the pin geometry, periodic 12 Å cell, PBE-D3(BJ), rcut = 2×
      CP2K R_CUTOFF 20 Å (~75.59 Bohr).
    status: passed
  - id: ac-004
    summary: e_total includes e_disp; forces FD-consistent
    type: runtime
    pass_when: |
      With d3bj on, GpwResult.e_total includes e_disp; analytic D3(BJ)
      forces match central FD of the same energy to 1e-7 Ha/Bohr (unit test).
    status: passed
  - id: ac-005
    summary: E2e regression registered under regressions/gpw/
    type: runtime
    pass_when: |
      regressions/gpw/gpw-pbe-d3bj-ice/reg_gpw-pbe-d3bj-ice exits 0;
      ATV_BUILD_GPW enables atv_add_regression(gpw gpw-pbe-d3bj-ice);
      public headers only in the regression main.
    status: passed
  - id: ac-006
    summary: Dual perf gate
    type: runtime
    pass_when: |
      Fair wall vs CP2K when binary available: ratio T_atv/T_cp2k ≤ 10 for a
      comparable job, OR perf_exempt with reason in REFERENCE.md.
    status: passed
out_of_scope:
  - Absolute e_tot micro-Ha parity with CP2K QS (OT/grid differences)
  - Full ice Ih supercell
  - D4 / ATM / CUDA
---

# Acceptance — gpw-pbe-d3bj-ice

Done means GPW can turn on PBE-D3(BJ) for a periodic water cell, match the
offline CP2K dispersion energy pin, and expose the result without Config bags.
