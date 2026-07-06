---
slug: hf-ao2mo-mp2
criteria:
  - id: ac-001
    summary: H2O/STO-3G MP2 correlation energy matches pinned PySCF in default tier
    type: scientific
    pass_when: |
      test_cpu_scf_mp2 (default ctest tier, no ATV_RUN_E2E_TESTS) contains a
      case computing mp2_energy for H2O/STO-3G (O 0 0 0.1173 / H 0 ±0.7572
      −0.4692 Å, rhf_solve tol 1e-12) and asserts
      |e_corr − (−0.0355456516)| <= 1e-8 Ha and e_tot == scf.energy + e_corr;
      case passes.
    status: pending
  - id: ac-002
    summary: He/STO-3G with zero virtuals returns exactly 0.0, never NaN
    type: scientific
    pass_when: |
      test_cpu_scf_mp2 (default tier) contains a He/STO-3G case asserting
      e_corr == 0.0 exactly and std::isfinite(e_corr); case passes.
    status: pending
  - id: ac-003
    summary: ao2mo_ovov output satisfies (ia|jb) = (jb|ia) permutation symmetry
    type: code
    pass_when: |
      A test_cpu_scf_mp2 unit case builds the ovov tensor for H2O/STO-3G via
      ao2mo_ovov and asserts max_{iajb} |(ia|jb) − (jb|ia)| <= 1e-12; case passes.
    status: pending
  - id: ac-004
    summary: mp2_energy guards throw std::invalid_argument (no fallback values)
    type: code
    pass_when: |
      test_cpu_scf_mp2 asserts std::invalid_argument is thrown for each of:
      (1) IntegralPack with unmaterialized ERI (compute_eri=false),
      (2) SCFResult with converged=false,
      (3) unrestricted result (n_alpha != n_beta); all three cases pass.
    status: pending
  - id: ac-005
    summary: E2E tier — all 10 (system, basis) pairs match pinned PySCF to 1e-8 Ha
    type: scientific
    pass_when: |
      With -DATV_RUN_E2E_TESTS=ON, test_cpu_scf_mp2 E2E cases pass for
      {He, H2, H2O, NH3, CH4} × {STO-3G, 6-31G} against the pinned E_corr
      table in the spec's Domain basis (tolerance 1e-8 Ha), 10/10.
    status: pending
  - id: ac-006
    summary: Example 13_mp2 builds under ATV_BUILD_RHF and prints E_HF, E_corr, E_MP2
    type: runtime
    pass_when: |
      cmake --build build --target 13_mp2 succeeds with -DATV_BUILD_RHF=ON and
      ./build/examples/13_mp2 exits 0 printing labeled E_HF, E_corr, E_MP2
      lines with E_HF + E_corr equal to E_MP2 at printed precision.
    status: pending
  - id: ac-007
    summary: ch07 rewritten as real MP2 chapter; quoted numbers pinned by test
    type: docs
    pass_when: |
      docs/tutorial/ch07-correlation.md no longer contains "not implemented"
      / "status chapter" claims for MP2; it contains the MP2 energy expression,
      DOIs 10.1103/PhysRev.46.618 and 10.1016/0009-2614(88)85250-3, a 13_mp2
      walk-through, and every quoted E_corr digit-for-digit equals a constant
      asserted in tests/cpu/scf/test_mp2.cpp.
    status: pending
  - id: ac-008
    summary: Public symbols in mp2.h carry Doxygen with units
    type: docs
    pass_when: |
      Every public symbol in include/atomiverse/cpu/scf/mp2.h (ao2mo_ovov,
      MP2Result and its fields, mp2_energy) has @brief plus @param/@return
      stating units (Hartree; dimensionless for counts).
    status: pending
  - id: ac-009
    summary: Default ctest gate preserved; full suite green
    type: runtime
    pass_when: |
      Default-tier ctest (ATV_RUN_E2E_TESTS=OFF) passes with test_cpu_scf_mp2
      included and total wall time within the ~16 s budget documented in
      tests/CMakeLists.txt; the E2E sweep is absent from the default run.
    status: pending
---

# Acceptance criteria

ac-001/ac-002 pin the default-tier scientific bar (parity + the no-virtuals NaN
edge). ac-003/ac-004 are the code-level unit contract for the transform and the
no-fallback guards. ac-005 is the full E2E-tier validation sweep. ac-006 covers
the runnable example, ac-007/ac-008 the tutorial-regression and Doxygen
deliverables, ac-009 the repository test-gate invariant.
