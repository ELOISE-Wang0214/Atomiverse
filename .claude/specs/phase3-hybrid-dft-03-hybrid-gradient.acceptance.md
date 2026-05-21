---
slug: phase3-hybrid-dft-03-hybrid-gradient
criteria:
  - id: ac-001
    summary: B3LYP analytic force matches PySCF grad.RKS to 1e-5 Ha/Bohr
    type: runtime
    pass_when: |
      A test in tests/cpu/scf/test_hybrid_gradient.cpp computes the B3LYP/6-31G*
      analytic force for H2O, CH4, and NH3 and asserts max|F_analytic - F_pyscf|
      <= 1e-5 Ha/Bohr against checked-in PySCF grad.RKS(xc='b3lyp') reference
      arrays; the test passes under ctest.
    status: pending
  - id: ac-002
    summary: PBE0 analytic force matches PySCF grad.RKS to 1e-5 Ha/Bohr
    type: runtime
    pass_when: |
      A test in tests/cpu/scf/test_hybrid_gradient.cpp computes the PBE0/6-31G*
      analytic force for H2O, CH4, and NH3 and asserts max|F_analytic - F_pyscf|
      <= 1e-5 Ha/Bohr against checked-in PySCF grad.RKS(xc='pbe0') reference
      arrays; the test passes under ctest.
    status: pending
  - id: ac-003
    summary: B3LYP/PBE0 analytic force matches moving-grid finite difference
    type: runtime
    pass_when: |
      A test in tests/cpu/scf/test_hybrid_gradient.cpp asserts the B3LYP and
      PBE0 analytic forces match scf_finite_difference_gradient (moving grid,
      h=1e-3 Bohr, SCF tol 1e-9) to <= 1e-5 Ha/Bohr; the test passes under ctest.
    status: pending
  - id: ac-004
    summary: Hybrid gradient is translationally invariant to 1e-8 Ha/Bohr
    type: runtime
    pass_when: |
      A test in tests/cpu/scf/test_hybrid_gradient.cpp asserts |sum_A F_A| per
      Cartesian axis is <= 1e-8 Ha/Bohr for both B3LYP and PBE0; the test
      passes under ctest.
    status: pending
  - id: ac-005
    summary: B3LYP/6-31G* geometry optimization converges and matches PySCF
    type: runtime
    pass_when: |
      A test in tests/cpu/scf/test_hybrid_gradient.cpp runs cpu::minimizer::Minimize
      with SCFEvaluator{b3lyp_config} + LBFGS on a small molecule, the optimizer
      reports converged, the optimized geometry RMS is <= 1e-3 Angstrom and the
      final energy <= 1e-6 Ha vs the checked-in PySCF optimized reference; the
      test passes under ctest.
    status: pending
  - id: ac-006
    summary: K-channel scales by alpha_hf in both 2e-gradient engines
    type: code
    pass_when: |
      accumulate_two_electron_gradient in src/cpu/scf/integral_derivatives.cpp
      multiplies factor_K by hybrid_exchange_fraction(config.xc_type) in both the
      libcint branch and the native Rys block-deriv branch, and the K loop is
      skipped when that value is 0.0 (LSDA/PBE).
    status: pending
  - id: ac-007
    summary: XC gradient routes B3LYP/PBE0 through the GGA vsigma path
    type: code
    pass_when: |
      accumulate_xc_gradient in src/cpu/scf/gradient.cpp determines is_gga via
      is_gga(config.xc_type) (true for PBE, B3LYP, PBE0), with no remaining
      XCType::PBE-only branch in the XC-gradient assembly.
    status: pending
  - id: ac-008
    summary: Full ctest suite stays green with no gradient regression
    type: runtime
    pass_when: |
      A clean build with -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON followed by
      ctest --output-on-failure reports zero failures, including all pre-existing
      HF, LSDA, and PBE gradient tests.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** — Reference parity. The binding bar is checked-in genuine PySCF reference data; the test file must record the PySCF version and `xc`/basis used, mirroring `test_def2svp_gradient.cpp`.
- **ac-003** — Independent oracle. Finite difference uses the moving-grid path (default `xc_grad_project_translation == false`); this catches any PySCF-convention mismatch without trusting PySCF.
- **ac-004** — Translational invariance is exact (not just `1e-5`) because the merged Becke grid-response term is already present; a `1e-8` residual would signal a regression in that term.
- **ac-005** — End-to-end composition proof. No new production code is exercised here beyond `b3lyp_config` flowing into the existing `Minimize`/`SCFEvaluator`/`LBFGS` stack.
- **ac-006 / ac-007** — Code-shape criteria pinning the two surgical changes; verified by reading the diff, independent of test outcomes.
- **ac-008** — Regression guard for the pure-DFT and HF paths, which must be bit-unchanged when `alpha_hf == 0`.
