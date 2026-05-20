---
slug: phase1-productization-02-ks-dshell-gradient
criteria:
  - id: ac-001
    summary: d-shell KS-gradient FD regression test exists and is registered
    type: code
    pass_when: |
      tests/cpu/scf/test_ks_dshell_gradient.cpp exists, contains a PBE
      RKS d-shell gradient-vs-FD test on water, and tests/CMakeLists.txt
      registers test_cpu_scf_ks_dshell_gradient under the ATV_BUILD_RHF
      block.
    status: pending
  - id: ac-002
    summary: analytic KS gradient matches FD to <= 1e-2 eV/A at PBE/6-31G* on water
    type: runtime
    pass_when: |
      In test_cpu_scf_ks_dshell_gradient, every force component of the
      analytic KS gradient for PBE/6-31G* water agrees with
      scf_finite_difference_gradient within 1e-2 eV/A; the test passes
      under ctest.
    status: pending
  - id: ac-003
    summary: analytic KS gradient matches FD to <= 1e-2 eV/A at PBE/def2-SVP on water
    type: runtime
    pass_when: |
      In test_cpu_scf_ks_dshell_gradient, every force component of the
      analytic KS gradient for PBE/def2-SVP water agrees with
      scf_finite_difference_gradient within 1e-2 eV/A; the test passes
      under ctest.
    status: pending
  - id: ac-004
    summary: max-force at PySCF-optimized water geometry drops to <= 1e-3 Ha/Bohr
    type: runtime
    evaluator_hint: scientific
    pass_when: |
      In test_cpu_scf_ks_dshell_gradient, the analytic KS gradient at
      the PySCF-optimized water geometry (PBE/def2-SVP) has max-force
      component <= 1e-3 Ha/Bohr (down from ~1 eV/Bohr); the test passes
      under ctest.
    status: pending
  - id: ac-005
    summary: tight analytic-vs-FD agreement on a fixed dense grid
    type: runtime
    evaluator_hint: scientific
    pass_when: |
      In test_cpu_scf_ks_dshell_gradient, with SCF converged to 1e-9 Ha,
      FD step h = 1e-3 Bohr, and the SAME fixed grid for + and -
      displacements, analytic vs FD agree to <= 1e-5 Ha/Bohr per
      component on PBE/def2-SVP water; the test passes under ctest.
    status: pending
  - id: ac-006
    summary: translational invariance and C2v symmetry of the KS gradient hold
    type: runtime
    pass_when: |
      In test_cpu_scf_ks_dshell_gradient, sum_A dE/dR_A = 0 per axis to
      <= 1e-7 and the H1/H2 water force components are mirror-equal; the
      test passes under ctest.
    status: pending
  - id: ac-007
    summary: full ctest suite stays green with no s/p or HF gradient regression
    type: runtime
    pass_when: |
      cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON && cmake
      --build build && ctest --test-dir build completes with all tests
      passing, including test_cpu_scf_gradient, test_cpu_scf_def2svp_gradient,
      and test_cpu_scf_gradient_vs_libcint.
    status: pending
---

# Acceptance criteria

- **ac-001** — The RED-before-GREEN entry point: the new test file exists and is
  wired into CTest. Verifiable by inspecting the two files; no run required.
- **ac-002 / ac-003** — The core fix. Analytic XC gradient matches a central-
  difference oracle to the project's existing d-shell tolerance (1e-2 eV/Å, the
  bar used by `test_def2svp_gradient.cpp`) at both d-bearing bases. Two criteria
  because 6-31G* and def2-SVP exercise different d-shell contraction patterns.
- **ac-004** — The user-facing symptom. The PySCF-optimized geometry must be
  recognized as near-stationary: max-force drops from ~1 eV/Bohr to
  <= 1e-3 Ha/Bohr, removing the reliance on Armijo backtracking.
- **ac-005** — The tighter scientific bar from the FD protocol. A fully correct
  GGA gradient agrees with FD to <= 1e-5 Ha/Bohr on a fixed dense grid. If grid
  noise prevents reaching 1e-5 the test records the achieved bound; ac-002/ac-003
  remain the binding correctness bar.
- **ac-006** — Physical sanity invariants independent of any reference code.
- **ac-007** — No-regression guard for the s/p analytic gradient path and the
  libcint-backed HF gradient cross-check.
