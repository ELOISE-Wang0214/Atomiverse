---
slug: hf-capability-gaps-08-stability
criteria:
  - id: ac-001
    summary: OrbitalHessian sigma action matches explicit dense (A±B) build
    type: code
    pass_when: |
      tests/cpu/scf/test_orbital_hessian.cpp passes: for H2/STO-3G and random
      vectors v, |dense_plus(ch)·v − sigma_plus(v,ch)|_max ≤ 1e-12 (idem minus),
      for both singlet and triplet channels.
    status: pending
  - id: ac-002
    summary: Singlet A diagonal equals eps_a − eps_i from the SlaterDeterminant
    type: scientific
    pass_when: |
      Unit test asserts dense ¹A diagonal == ε_a − ε_i [Ha] taken from the
      converged SlaterDeterminant (H2/STO-3G) to 1e-12.
    status: pending
  - id: ac-003
    summary: Equilibrium H2 RHF is internally and externally stable
    type: scientific
    pass_when: |
      tests/cpu/scf/test_stability.cpp: at R = 0.735 Å (STO-3G) both
      internal_eigenvalue and external_eigenvalue ≥ 0 and both *_stable flags true.
    status: pending
  - id: ac-004
    summary: Stretched H2 shows Seeger-Pople split - internal stable, external not
    type: scientific
    pass_when: |
      At R = 1.5 Å (STO-3G): internal_eigenvalue ≥ 0 AND external_eigenvalue < 0;
      regression reg_h2_rhf_stability matches the PySCF mf.stability() lowest
      external eigenvalue pinned in REFERENCE.md to 1e-6 Ha.
    status: pending
  - id: ac-005
    summary: Broken-symmetry mode lowers UHF below RHF at stretched geometry
    type: scientific
    pass_when: |
      In regressions/scf/h2_rhf_stability: UHF seeded from StabilityResult
      rotated_C converges with E_UHF < E_RHF by ≥ 1e-4 Ha at R = 1.5 Å,
      consistent with Szabo-Ostlund §3.8.7; regression runs ≤ 60 s.
    status: pending
  - id: ac-006
    summary: Kernel is a class - no build_* free functions in cpu/scf
    type: code
    pass_when: |
      OrbitalHessian/StabilityAnalysis are constructor-wired classes; grep gate
      tests/architecture/no_factory_scrub.sh passes with the new files present.
    status: pending
  - id: ac-007
    summary: Tutorial ch03 gains the SCF stability section with PySCF table
    type: docs
    pass_when: |
      docs/tutorial/ch03-open-shell.md contains a stability-analysis section
      (internal vs external narrative, stretched-H2 hands-on, RHF-vs-UHF
      energy-vs-R curve, chapter-end PySCF lowest-eigenvalue table).
    status: pending
  - id: ac-008
    summary: Full build and unit suite pass with the new tests registered
    type: runtime
    pass_when: |
      cmake --build build succeeds and ctest --test-dir build passes, including
      test_cpu_scf_orbital_hessian and test_cpu_scf_stability.
    status: pending
---

# Acceptance criteria

ac-001/002 钉 kernel 的数学正确性(σ-向量与稠密一致、对角结构);ac-003/004 钉
Seeger-Pople 分类本身(含"internal 稳定但 external 不稳定"陷阱与 PySCF 外部
oracle);ac-005 钉破缺对称初猜到 UHF 的桥;ac-006 钉 OOP/无工厂架构规则;
ac-007/008 钉教程交付与全套件绿。
