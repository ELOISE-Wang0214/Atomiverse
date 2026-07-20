---
slug: hf-capability-gaps-07-guess
criteria:
  - id: ac-001
    summary: New guess/mixer headers are classes; no factory free functions
    type: code
    pass_when: |
      internal/hueckel_guess.h, sap_guess.h, sap_data.h, ediis_adiis.h exist with
      constructor-based classes (HueckelGuess/SapGuess/EDIISMixer/ADIISMixer);
      tests/architecture/no_factory_scrub.sh passes; no cpu/scf build_* free
      function added.
    status: pending
  - id: ac-002
    summary: Selection is method members; enums live in scf_enums.h
    type: code
    pass_when: |
      InitGuess in scf_enums.h has Hueckel and SAP; MixerKind{CDIIS,EDIIS_CDIIS,
      ADIIS_CDIIS} in scf_enums.h; MeanFieldMethod has members mixer (default
      CDIIS) and mixer_switch_tol; no new *Config type appears.
    status: pending
  - id: ac-003
    summary: Unit tests pin Hueckel/SAP guess matrix math
    type: runtime
    pass_when: |
      test_cpu_scf_init_guess passes: Hueckel F_uv = 0.5*1.75*(H_uu+H_vv)*S_uv
      off-diagonal and VSIP diagonal element-wise; SAP F symmetric and guess
      density satisfies Tr(PS)=N_elec within 1e-8; unsupported Z throws.
    status: pending
  - id: ac-004
    summary: Unit tests pin EDIIS/ADIIS simplex solutions
    type: runtime
    pass_when: |
      test_cpu_scf_ediis_adiis passes: two-point closed-form optimal c matches
      to 1e-10 for both mixers; coefficients satisfy sum=1, c_i>=0; degenerate
      history returns latest F without NaN.
    status: pending
  - id: ac-005
    summary: Unit suite stays fast and E2E-free; defaults bit-identical
    type: runtime
    pass_when: |
      Full ctest passes with no SCF-to-convergence run added under tests/;
      existing test_cpu_scf_solver and test_cpu_scf_accelerators pass unchanged
      with default mixer=CDIIS.
    status: pending
  - id: ac-006
    summary: Regression pins iteration counts and hard-case convergence
    type: scientific
    pass_when: |
      regressions/scf/h2o_scf_convergence (<=60 s, regressions target only):
      H2O/def2-SVP iterations Hueckel<Hcore and SAP<Hcore; hard case where
      CDIIS-alone exceeds budget but ADIIS->CDIIS converges to the PySCF energy
      pinned in REFERENCE.md within 1e-6 Ha.
    status: pending
  - id: ac-007
    summary: ch18 half-chapter with domain-of-use pitfall and PySCF table
    type: docs
    pass_when: |
      docs/tutorial/ch18-scf-techniques.md contains the 初猜与 DIIS 家族 half:
      EDIIS/ADIIS far-from-convergence-only pitfall with ADIIS/EDIIS->CDIIS
      recipe, a molplot convergence-curve figure, and a chapter-end PySCF
      iteration-count table (scalars in tables, curves as figures).
    status: pending
  - id: ac-008
    summary: Doxygen and DOI references on all new symbols
    type: docs
    pass_when: |
      Every new public/internal symbol carries @brief/@param(units)/@return;
      hueckel_guess.h and sap_guess.h cite DOI 10.1021/acs.jctc.8b01089;
      ediis_adiis.h cites DOI 10.1063/1.1470195 and DOI 10.1063/1.3304922.
    status: pending
---

# Acceptance criteria

ac-001/002 为结构闸(无工厂、成员即参数、枚举归位 scf_enums.h);ac-003/004 为
单次调用数学闸(秒级单元);ac-005 守护默认行为零变化与 tests/ 无 E2E;ac-006 是
唯一的端到端科学闸(regressions,PySCF 钉定值);ac-007/008 为教程与文档闸。
