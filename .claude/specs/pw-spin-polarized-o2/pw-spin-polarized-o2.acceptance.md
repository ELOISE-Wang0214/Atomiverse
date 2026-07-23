---
spec: pw-spin-polarized-o2
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Public collinear spin API on cpu::pw::PW members
    type: code
    pass_when: |
      PW exposes spin_polarized (bool) and multiplicity (int) as public members
      (defaults false / 1). No public *Config bag, no make_*/build_* factory.
      PWResult exposes magnetization (electrons). Restricted default path unchanged.
    status: met
  - id: ac-002
    summary: Unit test for spin API and O2 magnetization identity
    type: runtime
    pass_when: |
      tests/cpu/planewave/test_pw_spin (or equivalent) passes: multiplicity-3
      O2 yields |magnetization - 2| below unit tol with fixed occupations.
    status: met
  - id: ac-003
    summary: Offline pin REFERENCE for O2 PW LSD e_tot and mag
    type: docs
    pass_when: |
      regressions/pw/pw-spin-polarized-o2/REFERENCE.md cites method (functional,
      ecut, cell, UPF, multiplicity), CP2K version if used else literature /
      fixed-occupation identity with explicit method string, and hard-coded
      e_tot + mag gates. No live third-party at ctest.
    status: met
  - id: ac-004
    summary: E2E regression vs offline pin
    type: runtime
    pass_when: |
      regressions/pw/pw-spin-polarized-o2 exits 0: energy within stated tol of
      pin and |mag - 2| within stated tol; public headers only.
    status: met
  - id: ac-005
    summary: Perf gate or exempt
    type: docs
    pass_when: |
      REFERENCE.md records T_atv/T_cp2k ≤ 10 when CP2K twin wall available,
      or perf_exempt=true with reason.
    status: met
out_of_scope:
  - k-point UKS / smearing
  - noncollinear spin
  - spin-polarized forces/stress
  - GPW UKS
  - CUDA
---

# Acceptance — pw-spin-polarized-o2

Done means Γ-point collinear spin-polarized plane-wave KS on O₂ triplet reports
magnetization 2 and total energy matching the offline pin under the dual gate.
