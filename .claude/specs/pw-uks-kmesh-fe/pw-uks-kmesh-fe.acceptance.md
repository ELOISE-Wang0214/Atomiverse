---
spec: pw-uks-kmesh-fe
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Multi-k collinear UKS uses existing public PW members
    type: code
    pass_when: |
      Regression and unit use cpu::pw::PW members spin_polarized, multiplicity,
      and k_mesh (≥ 2). No new public *Config bag, no make_*/build_* factory.
      PWResult.magnetization remains the mag observable.
    status: met
  - id: ac-002
    summary: Unit test multi-k magnetization identity
    type: runtime
    pass_when: |
      tests/cpu/planewave/test_pw_spin (or equivalent) includes a multi-k case
      (k_mesh ≥ 2) with multiplicity 3 O2-like fill yielding
      |magnetization - 2| below unit tol under fixed occupations.
    status: met
  - id: ac-003
    summary: Offline pin REFERENCE for multi-k UKS e_tot and mag
    type: docs
    pass_when: |
      regressions/pw/pw-uks-kmesh-fe/REFERENCE.md cites method (functional,
      ecut, cell, UPF, k_mesh≥2, multiplicity), offline pin source (CP2K
      version if used else frozen ATV / literature mag identity), and hard-coded
      e_tot + mag gates. No live third-party at ctest.
    status: met
  - id: ac-004
    summary: E2E regression vs offline pin
    type: runtime
    pass_when: |
      regressions/pw/pw-uks-kmesh-fe exits 0: SCF converged; energy within
      stated tol of pin; |mag - 2| within stated tol; public headers only;
      wall ≤ 60 s on the pin host class.
    status: met
  - id: ac-005
    summary: Perf gate or exempt
    type: docs
    pass_when: |
      REFERENCE.md records T_atv/T_cp2k ≤ 10 when CP2K twin wall available,
      or perf_exempt=true with reason.
    status: met
out_of_scope:
  - spin + metallic smearing / bulk Fe metal pin
  - noncollinear spin
  - spin-polarized forces/stress
  - GPW k-points
  - CUDA
---

# Acceptance — pw-uks-kmesh-fe

Done means collinear spin-polarized plane-wave KS with **k_mesh ≥ 2** on the
O₂ cubic cell reports magnetization 2 and total energy matching the offline
pin under the dual gate (correctness + perf/exempt).
