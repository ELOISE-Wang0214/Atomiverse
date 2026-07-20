---
slug: hf-capability-gaps-13-cphf
criteria:
  - id: ac-001
    summary: Second/skeleton derivative integrals match finite differences
    type: scientific
    pass_when: |
      tests/cpu/scf/test_cphf.cpp passes: on HeH+/STO-3G, d2S/d2T/d2V agree
      with central FD of dS/dT/dV to 1e-7, MO-basis S^(x)/F^(x) agree with FD
      of S/F to 1e-6 [Hartree/Bohr], and B^x sums to zero under rigid
      translation.
    status: pending
  - id: ac-002
    summary: CPHF iterative solve equals dense solve using the shared kernel
    type: code
    pass_when: |
      Unit test asserts iterative U^x equals the dense ¹(A+B) direct solve to
      1e-10, and the CPHF LHS action is bit-identical to OrbitalHessian
      sigma_plus(v, singlet) — the link-08 kernel, no reimplementation.
    status: pending
  - id: ac-003
    summary: Analytic Hessian is symmetric, translation-clean, matches FD of forces
    type: scientific
    pass_when: |
      tests/cpu/scf/test_analytic_hessian.cpp: H2/STO-3G analytic Hessian
      matches central FD of analytic forces elementwise to 1e-6 eV/A^2,
      H == H^T, and row translational sums vanish (|sum_j H_ij| ≤ 1e-6).
    status: pending
  - id: ac-004
    summary: FIRST GATE - analytic vs in-house FD frequencies within 1 cm^-1
    type: scientific
    pass_when: |
      regressions/scf/h2o_rhf_freq: H2O HF/6-31G* harmonic frequencies from
      rhf.hessian() and from the in-house FD Hessian (both through the
      existing vib mass-weighting/diagonalizer) agree mode-by-mode to
      ≤ 1 cm^-1.
    status: pending
  - id: ac-005
    summary: Frequencies match externally pinned PySCF hessian.RHF
    type: scientific
    pass_when: |
      Same regression compares against PySCF hessian.RHF frequencies pinned in
      REFERENCE.md (external run, values committed): every mode ≤ 1 cm^-1;
      regression registered via atv_add_regression(scf h2o_rhf_freq) and runs
      ≤ 60 s.
    status: pending
  - id: ac-006
    summary: Surface is a method output, no new driver, no factories
    type: code
    pass_when: |
      RHF::hessian(Frame&) returns vib::HessianResult beside forces(); no new
      driver or public result type is introduced; vib::Hessian is unmodified;
      CphfPack is constructor-wired under internal/ and
      tests/architecture/no_factory_scrub.sh passes.
    status: pending
  - id: ac-007
    summary: Tutorial ch20 gains the analytic-Hessian section with PySCF table
    type: docs
    pass_when: |
      docs/tutorial/ch20-hf-properties-environment.md contains the
      "解析 Hessian 与振动频率" section: CPHF orbital-response narrative and
      third reuse of the stability kernel, H2O hands-on, molvis normal-mode
      visualization, chapter-end PySCF frequency table (scalars in a table).
    status: pending
  - id: ac-008
    summary: Full build and unit suite pass with new tests registered
    type: runtime
    pass_when: |
      cmake --build build succeeds and ctest --test-dir build passes,
      including test_cpu_scf_cphf and test_cpu_scf_analytic_hessian; existing
      SCF/gradient tests stay green.
    status: pending
---

# Acceptance criteria

ac-001 钉两大经典 bug 源的逐项 FD 门禁(skeleton 导数与 B^x);ac-002 钉共享
kernel 的第三配置(不重写 σ-向量);ac-003/004 钉"第一道闸"——解析对自家
FD 的自洽(微型与 H₂O 两级);ac-005 钉外部 PySCF oracle;ac-006 钉方法输出
/无新 driver/无工厂的架构裁定;ac-007/008 钉教程交付与全套件绿。
