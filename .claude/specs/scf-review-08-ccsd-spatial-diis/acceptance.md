---
slug: scf-review-08-ccsd-spatial-diis
criteria:
  - id: ac-001
    summary: CcAmplitudeDiis exists as an OOP class, no Config bag / factory
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/internal/cc_diis.h declares class
      CcAmplitudeDiis with ctor(int max_vectors=8), step(amp,err),
      reset(), size(); no *Config type and no make_/build_ free function;
      unit test test_cpu_scf_cc_diis passes (DIIS extrapolation drives the
      toy error sequence below the hard-coded threshold; single-vector and
      near-singular-B edge cases return finite values).
    status: pending
  - id: ac-002
    summary: DIIS wired into the spin-orbital CCSD amplitude loop
    type: code
    pass_when: |
      src/cpu/scf/ccsd.cpp so_ccsd_kernel constructs a CcAmplitudeDiis and
      feeds the Jacobi increment (e = Δt = R/D) as the error vector,
      extrapolating [t1;t2] each iteration; the file builds and the CCSD
      unit + regression suite stays green.
    status: pending
  - id: ac-003
    summary: DIIS reaches the same fixed point as plain Jacobi
    type: runtime
    pass_when: |
      regressions/scf/ccsd-spatial-so-parity/reg_ccsd-spatial-so-parity.cpp
      asserts |E_corr(DIIS) − E_corr(Jacobi)| ≤ 1e-9 Ha on H2O/STO-3G at
      equal residual tolerance (both converged max|Δt| < 1e-8).
    status: pending
  - id: ac-004
    summary: Closed-shell spatial CCSD kernel reproduces the pinned E_corr
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/internal/ccsd_spatial.h +
      src/cpu/scf/ccsd_spatial.cpp implement ccsd_energy_spatial with
      spatial n_o²·n_v² storage (no spin expansion); unit test
      test_cpu_scf_ccsd_spatial passes — H2/STO-3G E_corr =
      −0.020524527092106 Ha (±1e-9), He/STO-3G E_corr == 0.
    status: pending
  - id: ac-005
    summary: ccsd_energy routes closed-shell to spatial; SO consumers unchanged
    type: code
    pass_when: |
      src/cpu/scf/ccsd.cpp ccsd_energy dispatches closed-shell RHF to
      ccsd_energy_spatial; ccsd_solve/uccsd_solve signatures and the
      CCSDSolution SO blocks are unchanged (eom_*_ccsd, ccsd_gradient,
      triples_energy/ccsd_t_energy, cc2_energy still compile and consume
      CCSDSolution unmodified).
    status: pending
  - id: ac-006
    summary: Spatial vs spin-orbital oracle-equivalence regression
    type: runtime
    pass_when: |
      regressions/scf/ccsd-spatial-so-parity/reg_ccsd-spatial-so-parity.cpp
      runs BOTH paths on H2O/STO-3G, asserts
      |E_corr_spatial − E_corr_so| ≤ 1e-9 Ha (both converged max|Δt| < 1e-8)
      and E_corr_spatial matches the embedded PySCF cc.CCSD golden ≤ 1e-6 Ha
      (literal value, no live third-party at ctest time); exits non-zero on
      failure and is registered in regressions/CMakeLists.txt.
    status: pending
  - id: ac-007
    summary: Existing CCSD/EOM/CCSD-gradient/CC2 regressions unchanged
    type: runtime
    pass_when: |
      The regressions target passes ccsd-h2-sto3g, ccsd-t-h2o-sto3g,
      ccsd-grad-h2o-sto3g, cc2-h2-sto3g, eom-ee-ccsd, eom-ea-ccsd,
      eom-ip-ccsd, eom-sf-ccsd-o2, eom-cc2-h2 with their pre-existing pinned
      values (no numeric edits to those cases).
    status: pending
  - id: ac-008
    summary: New symbols carry Doxygen with units
    type: docs
    pass_when: |
      CcAmplitudeDiis and ccsd_energy_spatial have @brief, @param (with
      units where physical, e.g. e_scf_ha [Hartree], tol dimensionless) and
      @return per the project Doxygen rule.
    status: pending
  - id: ac-009
    summary: Iteration-count / wall-time evidence recorded (advisory)
    type: performance
    pass_when: |
      The impl summary reports before/after wall-clock time and CCSD
      iteration counts (with vs without DIIS; spatial vs spin-orbital) on
      the H2O/STO-3G parity case, with no regression in E_corr; no fixed
      speedup factor is encoded as a gate.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 (DIIS, deliverable A)** — 新增 `CcAmplitudeDiis`（OOP、构造参数、无 config 包/工厂），并接入自旋轨道 CCSD 振幅循环，立即惠及所有 SO 消费者。
- **ac-003 (DIIS 定点不变量)** — 由 parity 回归同时验证：含 DIIS 与纯 Jacobi 收敛到同一定点（≤1e-9 Ha）。
- **ac-004 / ac-005 (空间轨道 CCSD, deliverable B)** — 空间核复现 pinned `E_corr`，且闭壳层能量入口默认走空间路径，同时 `ccsd_solve`/EOM/(T)/梯度/CC2/UCCSD 的 SO 契约逐位不变。
- **ac-006 (oracle 交叉校验，永久回归)** — 空间 vs 自旋轨道 `E_corr` ≤1e-9 Ha，并对齐嵌入的 PySCF golden；这是捕捉自旋求和系数错误的门。
- **ac-007** — 既有 CCSD/EOM/梯度/CC2 回归维持原 pin 值不变。
- **ac-008 / ac-009** — Doxygen 文档；迭代数/墙钟证据仅作 impl summary 记录（advisory，不设数值门）。
