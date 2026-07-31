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
    status: pass
  - id: ac-002
    summary: E_disp matches the CP2K pin
    type: runtime
    pass_when: |
      reg_gpw-d3-stress passes |e_disp − pin| ≤ 1e-8 Ha with d3=true on the
      periodic water dimer (kD3Cp2kPairCutoffBohr convention).
    status: pass
  - id: ac-003
    summary: D3 dispersion stress matches the CP2K difference pin
    type: runtime
    pass_when: |
      max diagonal |Δσ_atv − Δσ_cp2k| ≤ the REFERENCE-stated tolerance, where
      Δσ_atv = σ_an(d3) − σ_an(plain) and Δσ_cp2k = σ_num(D3) − σ_num(plain).
    status: pass
  - id: ac-004
    summary: Internal dual gate analytic vs multi-SCF strain-FD with d3=true
    type: runtime
    pass_when: |
      reg_gpw-d3-stress passes max diag |σ_an − σ_fd| ≤ 2e-5 Ha/Bohr³ on the
      same d3=true method; off-diagonal product ≈ 0; diag floor; kStress unit
      conversion check.
    status: pass
  - id: ac-005
    summary: Unit gates for the D3-zero stress fold
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_stress.cpp gains AnalyticMatchesFdWithD3 and
      D3ChangesGammaDiagonal (cheap fixture, seconds) and test_cpu_gpw_stress
      passes.
    status: pass
  - id: ac-006
    summary: Perf within 10x or exempt with stated reason
    type: performance
    pass_when: |
      REFERENCE.md records T_atv / T_cp2k for comparable stress work with
      ratio ≤ 10, or perf_exempt with a concrete reason.
    status: pass
  - id: ac-007
    summary: Queue, roadmap, tutorial updated
    type: docs
    pass_when: |
      docs/developer/cp2k-method-queue.md row gpw-d3-stress → done;
      docs/roadmap/dft.md dispersion-stress residual prose refreshed;
      docs/tutorial/methods/gpw-d3-stress.md exists.
    status: pass
out_of_scope:
  - D3(BJ) stress (sibling gpw-d3bj-stress)
  - ATM/C9 stress, shear components, dedicated multi-k D3 stress pin
---

# Verification (first green, 2026-07-31, wave 2)

- ac-001: REFERENCE.md carries CP2K 2026.1 pins (E_disp, both NUMERICAL
  tensors, derived Δσ), re-runnable from `tools/cp2k_oracle/inputs/gpw-d3-stress/`.
- ac-002: `reg_gpw-d3-stress` |E_disp − pin| = **5.53e−11 Ha** (tol 1e−8).
- ac-003: max diag |Δσ_atv − Δσ_cp2k| = **2.27e−10 Ha/Bohr³** (tol 2e−9;
  xx signal 1.71e−8 = 75× the residual).
- ac-004: dual gate max|σ_an − σ_fd| = **7.39e−07 Ha/Bohr³** (tol 2e−5);
  off-diag ~0; diag floor 5.83e−5; kStress conversion check green.
  Wall 47.7 s ≤ 60 s.
- ac-005: `test_cpu_gpw_stress` green (19 tests incl. AnalyticMatchesFdWithD3
  + D3ChangesGammaDiagonal, 11.6 s, ATV_TEST_HEAVY=ON).
- ac-006: T_atv/T_cp2k = 18.5/33.6 = **0.55** ≤ 10 (REFERENCE.md Speed).
- ac-007: queue row done (wave 2), roadmap §8 + header prose refreshed,
  tutorial page exists + indexed.
