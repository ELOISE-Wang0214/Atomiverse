---
slug: hf-capability-gaps-10-rohf
criteria:
  - id: ac-001
    summary: ROHF is a concrete CRTP method registered in methods.h, no factories
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/rohf.h declares `class ROHF : public
      MeanFieldMethod<ROHF>` with member `n_alpha`/`n_beta` occupation and no
      `build_*`/`make_*` free function anywhere in the diff; methods.h includes
      rohf.h; tests/architecture/no_factory_scrub.sh passes.
    status: pending
  - id: ac-002
    summary: O2 triplet ROHF/STO-3G total energy matches pinned PySCF to 1e-8 Ha
    type: scientific
    pass_when: |
      test_cpu_scf_rohf_solver asserts |E(ROHF, O2 R=1.2075 Å, (9,7), STO-3G) −
      E_pinned| ≤ 1e-8 Ha, where E_pinned is a PySCF 2.13.1 scf.ROHF value
      recorded verbatim (with geometry + command) in the test file header.
    status: pending
  - id: ac-003
    summary: CH3 doublet ROHF/STO-3G total energy matches pinned PySCF to 1e-8 Ha
    type: scientific
    pass_when: |
      test_cpu_scf_rohf_solver asserts |E(ROHF, planar CH3 R(C–H)=1.079 Å,
      (5,4), STO-3G) − E_pinned(PySCF 2.13.1)| ≤ 1e-8 Ha, provenance pinned in
      the test file header.
    status: pending
  - id: ac-004
    summary: ROHF is spin-pure; UHF contamination contrast holds
    type: scientific
    pass_when: |
      For converged O2 and CH3 ROHF runs, |spin contamination| ≤ 1e-10
      (⟨S²⟩ = 2.0 resp. 0.75 exactly); a same-geometry O2 UHF run in the same
      test file shows spin contamination > 1e-3.
    status: pending
  - id: ac-005
    summary: Orbital energies match PySCF only under Guest–Saunders convention
    type: scientific
    pass_when: |
      Occupied orbital energies of O2/CH3 match pinned PySCF canonical-ROHF
      values to 1e-6 Ha, and neither tests nor acceptance gate on any
      non-PySCF-convention orbital energies (comment in test cites the
      convention).
    status: pending
  - id: ac-006
    summary: Closed-shell degeneracy and guarded unsupported paths
    type: code
    pass_when: |
      Unit tests pass asserting: ROHF(1,1) on H2/STO-3G equals RHF energy to
      1e-8 Ha; ROHF with n_alpha < n_beta throws std::invalid_argument;
      forces()/mp2() on ROHF throw rather than silently taking a
      Restricted/Unrestricted branch.
    status: pending
  - id: ac-007
    summary: Regression o2_rohf_sto3g passes as out-of-library consumer ≤60 s
    type: runtime
    pass_when: |
      regressions/scf/o2_rohf_sto3g builds behind ATV_BUILD_REGRESSIONS using
      only public headers + link atomiverse, passes via the `regressions`
      target in ≤60 s, and REFERENCE.md cites Roothaan 1960 (DOI
      10.1103/RevModPhys.32.179) plus the pinned PySCF value and command.
    status: pending
  - id: ac-008
    summary: Tutorial ch03 ROHF section and Doxygen convention note shipped
    type: docs
    pass_when: |
      docs/tutorial/ch03-open-shell.md gains a ROHF section (spin-purity
      narrative, coupling-coefficient ambiguity stated honestly, O2/CH3
      hands-on, chapter-end scalar table rows E_total and ⟨S²⟩ vs PySCF, no
      bar charts) and every new public symbol in rohf.h /
      rohf_couple_stage.h carries @brief/@param(units)/@return, with the
      Guest–Saunders convention recorded in the ROHF class doc.
    status: pending
---

# Acceptance criteria

ac-002/003/005 的钉定值在实现期由外部 PySCF 2.13.1 一次性生成并连同几何、
命令行写入测试文件头与 REFERENCE.md;总能量是可移植观测量(1e-8 Ha 门槛),
轨道能仅在 PySCF canonical ROHF(Guest–Saunders)约定下比较(1e-6 Ha)。
ac-004 的自旋纯性是解析恒等式(Roothaan 1960),不依赖任何外部引擎。
