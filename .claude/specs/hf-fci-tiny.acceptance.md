---
slug: hf-fci-tiny
criteria:
  - id: ac-001
    summary: He/STO-3G boundary — FCI equals RHF, e_corr exactly 0.0, never NaN
    type: scientific
    pass_when: |
      test_cpu_scf_fci (default ctest tier) contains a
      He/STO-3G case (He 0 0 0, rhf_solve tol 1e-12, compute_eri=true) asserting
      e_corr == 0.0 exactly, std::isfinite(e_fci), and
      |e_fci − scf.energy| <= 1e-10 Ha; case passes.
    status: pending
  - id: ac-002
    summary: H2/STO-3G FCI energy matches pinned PySCF fci.FCI to 1e-8 Ha (validation benchmark)
    type: scientific
    pass_when: |
      The val_fci_vs_pyscf validation benchmark
      (benchmarks/qc/fci/validation/test_fci_vs_pyscf.cpp, ATV_BUILD_BENCHMARKS)
      computes fci_energy for H2/STO-3G (H 0 0 0 / H 0 0 0.7414 Å, scf_solve
      rhf_config tol 1e-12) and asserts |e_fci − <pinned PySCF fci.FCI value>|
      <= 1e-8 Ha with the pinned constant regenerated from a PySCF oracle; case
      passes. (External-engine parity is a benchmark, not a unit test.)
    status: pending
  - id: ac-003
    summary: H_CI is Hermitian and determinant count equals the binomial product
    type: code
    pass_when: |
      test_cpu_scf_fci asserts, for H2/STO-3G, max |H_CI − H_CIᵀ| <= 1e-12 before
      diagonalization and that the enumerated determinant count equals
      C(n_mo,n_alpha)*C(n_mo,n_beta) (== 4 here); case passes.
    status: pending
  - id: ac-004
    summary: fci_energy guards throw std::invalid_argument (no fallback values)
    type: code
    pass_when: |
      test_cpu_scf_fci asserts std::invalid_argument is thrown for each of:
      (1) IntegralPack with unmaterialized ERI (compute_eri=false),
      (2) SCFResult with converged=false,
      (3) open-shell result (n_alpha != n_beta),
      (4) a system exceeding the TINY cap (n_mo > 14 or N_det > 2000);
      all four cases pass.
    status: pending
  - id: ac-005
    summary: Variational and correlation-sign invariants hold
    type: scientific
    pass_when: |
      test_cpu_scf_fci asserts e_fci <= scf.energy + 1e-10 Ha (variational bound)
      and e_corr <= 0 for the H2/STO-3G case; case passes.
    status: pending
  - id: ac-006
    summary: MP2 overbinds — e_corr_mp2 <= e_corr_fci < 0 on the 2-level H2/STO-3G pair
    type: scientific
    pass_when: |
      test_cpu_scf_fci computes both mp2_energy and fci_energy for H2/STO-3G and
      asserts e_corr_mp2 <= e_corr_fci and e_corr_fci < 0 (NOT |e_fci| >= |e_mp2|);
      case passes.
    status: pending
  - id: ac-007
    summary: ao2mo_full has 8-fold symmetry and mo_hcore equals Cᵀ·Hcore·C
    type: code
    pass_when: |
      A test_cpu_scf_fci unit case builds the full MO tensor via ao2mo_full for
      H2/STO-3G and asserts 8-fold permutation symmetry (e.g.
      max |g(pqrs) − g(qprs)|, |g(pqrs) − g(rspq)| <= 1e-12), and asserts
      mo_hcore output equals xt::linalg::dot(transpose(C), dot(Hcore, C)) and is
      symmetric to 1e-12; case passes.
    status: pending
  - id: ac-008
    summary: direct_spin1 fed engine MO integrals matches the mean-field FCI path (validation benchmark)
    type: scientific
    pass_when: |
      The val_fci_vs_pyscf validation benchmark cross-checks fci_energy against
      PySCF fci.direct_spin1.kernel(h1e, h2e, norb, nelec, ecore=V_nuc) fed the
      engine's own mo_hcore/ao2mo_full integrals to 1e-8 Ha, plus any larger tiny
      systems; runs under ATV_BUILD_BENCHMARKS, absent from the unit gate.
    status: pending
  - id: ac-009
    summary: Example 14_fci builds under ATV_BUILD_RHF and prints E_RHF, E_FCI, E_corr
    type: runtime
    pass_when: |
      cmake --build build --target 14_fci succeeds with -DATV_BUILD_RHF=ON and
      ./build/examples/14_fci exits 0 printing labeled E_RHF, E_FCI, E_corr lines
      for H2/STO-3G with E_RHF + E_corr equal to E_FCI at printed precision.
    status: pending
  - id: ac-010
    summary: ch07 gains a full-CI section pinned by test; new symbols carry Doxygen
    type: docs
    pass_when: |
      docs/tutorial/ch07-correlation.md contains a full-CI section (determinant
      expansion, exponential cost, reference-not-production framing) whose quoted
      He E_corr=0 boundary and H2/STO-3G number equal constants asserted in
      tests/cpu/scf/test_fci.cpp, and cite Slater/Condon + Knowles-Handy DOIs;
      every public symbol in fci.h (FCIResult + fields, fci_energy) and the two new
      mp2.h symbols (ao2mo_full, mo_hcore) has @brief plus @param/@return stating
      units (Hartree; dimensionless for counts).
    status: pending
  - id: ac-011
    summary: Default ctest gate preserved; full suite green
    type: runtime
    pass_when: |
      Default-tier ctest passes with test_cpu_scf_fci
      included and total wall time within the ~16 s budget documented in
      tests/CMakeLists.txt; the direct_spin1 E2E cross-check is absent from the
      default run.
    status: pending
---

# Acceptance criteria

ac-001/ac-002 pin the default-tier scientific bar: the He no-virtuals identity
(FCI≡RHF, exact zero) and the H2/STO-3G parity against pinned PySCF FCI.
ac-003/ac-004 are the code-level unit contract — CI-matrix Hermiticity plus the
binomial determinant count, and the four no-fallback guards (including the
load-bearing TINY cap). ac-005/ac-006 are the physics invariants: the variational
bound / correlation sign, and the MP2-overbinds cross-check that ties this
milestone to hf-ao2mo-mp2 (deliberately NOT the false "FCI recovers ≥ |MP2|"
claim). ac-007 is the transform contract for the two new ao2mo-family symbols.
ac-008 is the E2E direct_spin1 cross-check that isolates ao2mo bugs from FCI bugs.
ac-009 covers the runnable example, ac-010 the tutorial-regression + Doxygen
deliverables, ac-011 the repository test-gate invariant.
