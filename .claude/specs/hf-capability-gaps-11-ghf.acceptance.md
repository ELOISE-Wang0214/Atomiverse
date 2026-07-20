---
slug: hf-capability-gaps-11-ghf
criteria:
  - id: ac-001
    summary: Spinor SlaterDeterminant invariants validated
    type: code
    pass_when: |
      tests/cpu/scf/test_ghf.cpp spinor cases pass: generalized constructor
      stores C(2n x N_mo); density() is 2n x 2n; spin-block accessors return
      the four blocks; validate(S) enforces Tr(P·S⊗1₂)=N and idempotency;
      dimension/occupation violations throw std::invalid_argument.
    status: pending
  - id: ac-002
    summary: GHF collapses to UHF on unfrustrated H2/STO-3G
    type: code
    pass_when: |
      In tests/cpu/scf/test_ghf.cpp, converged GHF on H2/STO-3G satisfies
      |E_GHF - E_UHF| <= 1e-8 Ha against a UHF run on the same geometry.
    status: pending
  - id: ac-003
    summary: Non-collinear GHF below UHF on triangular H3, matches PySCF
    type: scientific
    pass_when: |
      In tests/cpu/scf/test_ghf.cpp, equilateral H3/STO-3G GHF with
      collinearity-breaking guess gives E_GHF <= E_UHF - 1e-6 Ha and
      |E_GHF - E_pinned(PySCF scf.GHF)| <= 1e-8 Ha (pinned value recorded
      in regressions/scf/h3_ghf_sto3g/REFERENCE.md).
    status: pending
  - id: ac-004
    summary: GHF registered as concrete CRTP method, no factory/config debt
    type: code
    pass_when: |
      methods.h includes ghf.h; src/cpu/scf/mean_field_method.cpp contains
      `template class MeanFieldMethod<GHF>;`; parameters are members
      (n_electrons, guess_spin_rotation); tests/architecture/
      no_factory_scrub.sh and test_stage_typing_scrub pass.
    status: pending
  - id: ac-005
    summary: Regression h3_ghf_sto3g passes against pinned oracle in <=60 s
    type: scientific
    pass_when: |
      `regressions` target runs reg_h3_ghf_sto3g as an out-of-library
      consumer (public headers + link atomiverse only), checks the
      REFERENCE.md pinned PySCF values, and completes in <= 60 s; it is
      not registered in default ctest.
    status: pending
  - id: ac-006
    summary: Tutorial ch03 gains the GHF spin-symmetry-breaking section
    type: docs
    pass_when: |
      docs/tutorial/ch03-open-shell.md contains a "GHF 与自旋对称破缺"
      section with the RHF→UHF→ROHF→GHF constraint-removal ladder, an
      honest Fukutome-classes mention with the 1981 DOI, the triangular-H3
      hands-on with molvis spin-vector rendering, and chapter-end PySCF
      comparison table rows (scalars in tables).
    status: pending
  - id: ac-007
    summary: New public symbols carry Doxygen with units
    type: docs
    pass_when: |
      Every public symbol added in ghf.h, the new SlaterDeterminant
      accessors, and SpinType::Generalized carries @brief, @param with
      units, and @return.
    status: pending
  - id: ac-008
    summary: Full unit suite stays green and fast
    type: runtime
    pass_when: |
      ctest --test-dir build passes with ATV_BUILD_HF=ON including
      test_cpu_scf_ghf; the whole tests/ suite still completes in seconds
      (no e2e added to ctest).
    status: pending
---

# Acceptance criteria

ac-001–ac-004、ac-008 由 /mol:impl 在实现回路中验证(code/runtime);ac-003、
ac-005 为领域验证,判据是外部钉定的 PySCF `scf.GHF`/`scf.UHF` 数值(oracle
依赖仅为数值,非构建依赖);ac-006、ac-007 为文档判据。ac-003 的塌缩/非共线
双测试共同编码 GHF 的核心物理陷阱:无自旋轨道项时 GHF 必须塌缩回 UHF,仅在
自旋阻挫体系上凭非共线性下降。SlaterDeterminant 旋量推广是本环架构决策,
实现前须与 scf-slater-determinant 链定案(扩展而非分叉)。
