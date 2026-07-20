---
slug: scf-slater-determinant-01-type
criteria:
  - id: ac-001
    summary: SlaterDeterminant is a validating public value type in cpu/scf
    type: code
    evaluator_hint: "inspect: include/atomiverse/cpu/scf/slater_determinant.h"
    pass_when: |
      include/atomiverse/cpu/scf/slater_determinant.h defines cpu::scf::SlaterDeterminant with
      per-spin-channel storage {RealMatrix C, RealVector epsilon, RealVector occupation}, a nested
      Spin{alpha,beta} enum, restricted single-channel form whose beta accessors forward to alpha,
      and constructors that throw std::invalid_argument with numeric context on any of: dimension
      mismatch, f_i < 0, f_i > 1, or |Sum f_i^sigma - N^sigma| > 1e-12*max(1,N). No make_*/factory
      function and no *Config type is introduced.
    status: verified
    last_checked: 2026-07-20
  - id: ac-002
    summary: density(Spin) delegates to linalg::build_density_matrix and matches GEMM
    type: code
    evaluator_hint: "target: test_cpu_scf_slater_determinant"
    pass_when: |
      density(Spin) is implemented via linalg::build_density_matrix (include/atomiverse/cpu/linalg/
      density.h, unrenamed), and the unit test asserts element-wise agreement <= 1e-12 with a
      hand-written C_occ * C_occ^T GEMM for Aufbau occupation, plus a non-Aufbau (1,0,1) and a
      fractional-occupation case against manual Sum f_i c_i c_i^T.
    status: verified
    last_checked: 2026-07-20
  - id: ac-003
    summary: density_total() correct in both spin forms
    type: code
    evaluator_hint: "target: test_cpu_scf_slater_determinant"
    pass_when: |
      Unit tests pass asserting density_total() == 2*density(alpha) (<= 1e-12) for the restricted
      form and density_total() == density(alpha) + density(beta) for the unrestricted form.
    status: verified
    last_checked: 2026-07-20
  - id: ac-004
    summary: expect(O) equals the manual trace and documents its one-electron limit
    type: code
    evaluator_hint: "target: test_cpu_scf_slater_determinant + header Doxygen"
    pass_when: |
      expect(O) matches a manually computed Tr(P_total * O) to <= 1e-12 in the unit test, and its
      Doxygen comment states both caveats: valid for one-electron operators only (two-electron
      needs the two-particle density matrix), and the mean-field total energy is not expect(F)
      (E_HF = 1/2 Tr[P(H_core+F)] + E_nn; Sum eps double-counts repulsion).
    status: verified
    last_checked: 2026-07-20
  - id: ac-005
    summary: homo/lumo/orbital_energy reproduce observer semantics
    type: code
    evaluator_hint: "target: test_cpu_scf_slater_determinant vs orbital_info_observer.cpp:43-60"
    pass_when: |
      On a known ascending epsilon sequence, homo(Spin) == eps(n_occ-1) and lumo(Spin) == eps(n_occ)
      in Hartree; with no virtual orbitals lumo() returns homo(); orbital_energy(Spin,i) returns
      eps(i) and throws on out-of-range i. Tests pass.
    status: verified
    last_checked: 2026-07-20
  - id: ac-006
    summary: with_occupation is a ref-qualified value-returning DeltaSCF/MOM entry
    type: code
    evaluator_hint: "inspect ref-qualifiers; target: test_cpu_scf_slater_determinant"
    pass_when: |
      with_occupation(Spin,from,to) exists as a const&/&& ref-qualified pair (with_gradient
      precedent, mean_field_method.h:180-187); the const& path returns a new re-validated value
      leaving the original unchanged, the && path chains; f_from != 1 or f_to != 0 throws
      std::invalid_argument. Tests pass.
    status: verified
    last_checked: 2026-07-20
  - id: ac-007
    summary: validate(S) enforces Tr(PS)=N and integer-occupation idempotency
    type: scientific
    evaluator_hint: "target: test_cpu_scf_slater_determinant (orthonormal toy basis)"
    pass_when: |
      validate(S, tol=1e-8) passes on an orthonormal toy basis (S=I) for integer occupations
      (both Tr(P^sigma S) = N^sigma and P^sigma S P^sigma = P^sigma per channel), skips the
      idempotency check for fractional occupations while still enforcing Tr(PS)=N, and throws
      with channel + deviation context for a non-orthonormal C against S=I — matching the
      invariants of Szabo & Ostlund Ch. 3 / the f-in-{0,1} single-determinant criterion.
    status: pending
  - id: ac-008
    summary: Zero-consumer wiring — no existing production code touched
    type: code
    evaluator_hint: "grep: SlaterDeterminant over include/ src/"
    pass_when: |
      A grep for SlaterDeterminant across include/ and src/ matches only
      include/atomiverse/cpu/scf/slater_determinant.h; git diff shows the only modified existing
      file is tests/CMakeLists.txt (test registration); scf_result.h, orbital_info_observer.*,
      mean_field_method.*, and linalg/density.h are unchanged.
    status: verified
    last_checked: 2026-07-20
  - id: ac-009
    summary: Full build and fast unit suite green with the new seconds-fast test
    type: runtime
    pass_when: |
      cmake --build build succeeds and ctest --test-dir build --output-on-failure is green with
      test_cpu_scf_slater_determinant registered and passing; the new test uses only toy matrices
      (n <= 4), no full SCF solve.
    status: verified
    last_checked: 2026-07-20
  - id: ac-010
    summary: Doxygen coverage with units on every public symbol
    type: docs
    evaluator_hint: "inspect: slater_determinant.h"
    pass_when: |
      Every public symbol in slater_determinant.h carries @brief, @param (with units), and @return;
      epsilon and homo/lumo/orbital_energy are documented in Hartree; C/f/P documented as
      dimensionless; expect() return units documented as those of O.
    status: pending
---

# Acceptance criteria

ac-001–ac-006 与 ac-010 为类型面（构造校验、密度/期望/访问器/ΔSCF 入口、文档）；ac-007 是 Domain basis 第 6 条不变量的科学验证；ac-008 是本环"零消费者接线"的硬门；ac-009 是 runtime 关口（快套件秒级，无 E2E）。
