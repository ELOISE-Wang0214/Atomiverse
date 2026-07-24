---
slug: scf-review-07-posthf-perf
criteria:
  - id: ac-001
    summary: eri_dense assembles via native::eri_block + OpenMP, not scalar native::eri
    type: code
    pass_when: |
      src/cpu/scf/two_electron.cpp IntegralPack::eri_dense calls
      blocks_from_shells and native::eri_block, contains an
      `#pragma omp parallel for` over the block-pair loop, and no longer
      calls native::eri; the target builds with -DATV_BUILD_HF=ON.
    status: verified
    last_checked: 2026-07-24
  - id: ac-002
    summary: block-assembled dense ERI matches scalar reference ≤1e-12 per element
    type: code
    evaluator_hint: "ctest: test_cpu_scf_two_electron"
    pass_when: |
      Unit test test_cpu_scf_two_electron passes: on the STO-3G p-shell
      fixture max element-wise |eri_dense(block) − native::eri(reference)|
      ≤ 1e-12, including the pure-s and (pp|pp) cases.
    status: verified
    last_checked: 2026-07-24
  - id: ac-003
    summary: ri_mp2_energy contracts via a single GEMM, no scalar P-loop
    type: code
    pass_when: |
      src/cpu/scf/mp2.cpp ri_mp2_energy contains a single xt::linalg::dot
      over the (n_occ·n_vir)×n_aux B_ov tensor and no `for (int P = 0; P <
      n_aux` accumulation inside the O(no²nv²) amplitude loop, and no
      per-P element-wise M[μν] reconstruction; the target builds.
    status: verified
    last_checked: 2026-07-24
  - id: ac-004
    summary: RI-MP2 correlation energy unchanged ≤1e-10 Ha vs golden
    type: code
    evaluator_hint: "ctest: test_cpu_scf_ri_mp2"
    pass_when: |
      Unit test test_cpu_scf_ri_mp2 passes: ri_mp2_energy e_corr matches
      the embedded hard-coded golden within 1e-10 Ha.
    status: verified
    last_checked: 2026-07-24
  - id: ac-005
    summary: ao2mo entry N⁴ copy eliminated in mp2 + ccsd shared entry copies
    type: code
    pass_when: |
      Neither ao2mo_mixed_ (src/cpu/scf/mp2.cpp) nor ao2mo_uhf_block
      (src/cpu/scf/ccsd.cpp) contains the initializer
      `xt::xtensor<double, 4> cur = eri_ao;`; test_cpu_scf_mp2 passes and
      reg_ccsd-h2-sto3g exits 0.
    status: verified
    last_checked: 2026-07-24
  - id: ac-006
    summary: block-ERI MP2 regression reproduces PySCF-pinned H2O/STO-3G e_corr
    type: runtime
    evaluator_hint: "regressions: reg_posthf-eri-block-mp2"
    pass_when: |
      Under -DATV_BUILD_REGRESSIONS=ON the `regressions` target builds and
      runs reg_posthf-eri-block-mp2, which reproduces the REFERENCE.md
      PySCF-pinned H2O/STO-3G MP2 correlation energy (literal in the source,
      no live PySCF) within 1e-6 Ha and exits 0.
    status: verified
    last_checked: 2026-07-24
  - id: ac-007
    summary: all post-HF scf regressions reproduce cited values unchanged
    type: runtime
    evaluator_hint: "regressions: scf post-HF family"
    pass_when: |
      The `regressions` target builds and every scf post-HF case exits 0:
      reg_h2o_mp2_sto3g, reg_h2o_ri_mp2, reg_ump2-o2-sto3g,
      reg_romp2-o2-sto3g, reg_scs-mp2-h2o-sto3g, reg_sos-mp2-h2o-sto3g,
      reg_frozen-core-mp2-h2o, reg_mp3-h2-sto3g, reg_mp4-h2-sto3g,
      reg_ccsd-h2-sto3g, reg_ccsd-t-h2o-sto3g, reg_cc2-h2-sto3g,
      reg_cisd-h2co-sto3g, reg_adc2-h2co-sto3g, reg_adc2x-h2co-sto3g,
      reg_ip-adc2-h2o-sto3g, reg_ea-adc2-h2-sto3g, reg_gf2-h2-sto3g,
      reg_g0w0-h2-sto3g, reg_bse-h2-sto3g, reg_casscf-h2-sto3g,
      reg_casci-h2-sto3g, reg_nevpt2-h2-sto3g, reg_eom-ee-ccsd,
      reg_eom-ea-ccsd, reg_eom-ip-ccsd, reg_eom-sf-ccsd-o2, reg_eom-cc2-h2,
      reg_direct-rpa-h2-sto3g, reg_b2plyp-h2o, reg_delta-scf-n2.
    status: verified
    last_checked: 2026-07-24
  - id: ac-008
    summary: Doxygen with units on eri_dense and ri_mp2_energy
    type: docs
    pass_when: |
      IntegralPack::eri_dense (internal/integrals.h) and ri_mp2_energy
      (mp2.h) each carry @brief, @param with units (Hartree/Bohr), and
      @return in their declaring headers.
    status: verified
    last_checked: 2026-07-24
---

# Acceptance criteria

- **ac-001 / ac-002（item 17）**：稠密 ERI 走壳块引擎 + OpenMP，且块装配张量与标量参考逐元素 ≤1e-12——静态 grep 保证结构、单测保证数值等价。
- **ac-003 / ac-004（item 18）**：RI-MP2 收缩为单次 GEMM 且无标量 P-循环，能量对 golden ≤1e-10 Ha（浮点重结合界 N·ε，非机器精度）；拟合度量冻结不改。
- **ac-005（item 21）**：两处共享入口 N⁴ 拷贝消除，行为由既有 MP2 单测与 CCSD regression 守卫。
- **ac-006**：本 spec 的唯一 regression，经块装配稠密 ERI 端到端复现 PySCF 固定的 MP2 相关能（literal，无运行期第三方）。
- **ac-007**：全部 post-HF regression 家族复现引用值不变——item 17 张量等价的端到端总门禁。
- **ac-008**：公共/内部符号 Doxygen 带单位。
- 说明：性能墙钟前后对比为咨询性证据（impl 总结），**不**是 acceptance 门禁，故此契约不含 performance 型判据。
