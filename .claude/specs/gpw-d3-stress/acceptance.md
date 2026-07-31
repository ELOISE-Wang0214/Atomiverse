---
spec: gpw-d3-stress
created: 2026-07-31
criteria:
  - id: ac-001
    summary: Offline CP2K D3-zero + plain NUMERICAL stress pins recorded
    type: docs
    pass_when: |
      regressions/gpw/gpw-d3-stress/REFERENCE.md records CP2K 2026.1 version,
      inputs, total energies, Dispersion energy, and both NUMERICAL stress
      tensors (plus the derived diagonal Δσ pin), re-runnable from
      tools/cp2k_oracle/inputs/gpw-d3-stress/. Pins are numbers CP2K printed,
      never seeded from ATV.
    status: pending
  - id: ac-002
    summary: E_disp matches the CP2K pin
    type: runtime
    pass_when: |
      reg_gpw-d3-stress passes |e_disp − pin| ≤ 1e-8 Ha with d3=true on the
      periodic water dimer (kD3Cp2kPairCutoffBohr convention).
    status: pending
  - id: ac-003
    summary: D3 dispersion stress matches the CP2K difference pin
    type: runtime
    pass_when: |
      max diagonal |Δσ_atv − Δσ_cp2k| ≤ the REFERENCE-stated tolerance, where
      Δσ_atv = σ_an(d3) − σ_an(plain) and Δσ_cp2k = σ_num(D3) − σ_num(plain).
    status: pending
  - id: ac-004
    summary: Internal dual gate analytic vs multi-SCF strain-FD with d3=true
    type: runtime
    pass_when: |
      reg_gpw-d3-stress passes max diag |σ_an − σ_fd| ≤ 2e-5 Ha/Bohr³ on the
      same d3=true method; off-diagonal product ≈ 0; diag floor; kStress unit
      conversion check.
    status: pending
  - id: ac-005
    summary: Unit gates for the D3-zero stress fold
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_stress.cpp gains AnalyticMatchesFdWithD3 and
      D3ChangesGammaDiagonal (cheap fixture, seconds) and test_cpu_gpw_stress
      passes.
    status: pending
  - id: ac-006
    summary: Perf within 10x or exempt with stated reason
    type: performance
    pass_when: |
      REFERENCE.md records T_atv / T_cp2k for comparable stress work with
      ratio ≤ 10, or perf_exempt with a concrete reason.
    status: pending
  - id: ac-007
    summary: Queue, roadmap, tutorial updated
    type: docs
    pass_when: |
      docs/developer/cp2k-method-queue.md row gpw-d3-stress → done;
      docs/roadmap/dft.md dispersion-stress residual prose refreshed;
      docs/tutorial/methods/gpw-d3-stress.md exists.
    status: pending
out_of_scope:
  - D3(BJ) stress (sibling gpw-d3bj-stress)
  - ATM/C9 stress, shear components, dedicated multi-k D3 stress pin
---
