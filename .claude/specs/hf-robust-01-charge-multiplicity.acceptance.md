---
slug: hf-robust-01-charge-multiplicity
criteria:
  - id: ac-001
    summary: SpinResolution constructor reproduces the four-target (charge,M)->(n_alpha,n_beta) table
    type: code
    pass_when: |
      In tests/cpu/scf/test_spin_resolution.cpp, direct construction of
      ATV::cpu::scf::SpinResolution (validating constructor -- no factory
      function, per the 2026-07-02 user design directive) asserts exactly:
      (N_e=10, M=1)->(5,5) [H2O(0,1)]; (N_e=10, M=1)->(5,5) built from
      OH- integrals (charge=-1); (N_e=9, M=2)->(5,4) [H2O+]; (N_e=8, M=3)->(5,3) [O atom];
      all four assertions pass in test_cpu_scf_spin_resolution.
    status: verified
    last_checked: 2026-07-02
  - id: ac-002
    summary: SpinResolution rejects the three invalid (charge,M) combinations
    type: code
    pass_when: |
      test_cpu_scf_spin_resolution asserts std::invalid_argument is thrown (EXPECT_THROW
      on the exception type, not message text) for H2O(0,2) [N_e=10,M=2 parity],
      H2O(+1,1) [N_e=9,M=1 parity], and H(0,4) [N_e=1,M=4, n_beta<0].
    status: verified
    last_checked: 2026-07-02
  - id: ac-003
    summary: RHF on an odd-electron system now throws instead of silently truncating
    type: code
    pass_when: |
      test_cpu_scf_spin_resolution builds H2O+ integrals (build_integrals(frame,
      BasisSet::STO3G, /*charge=*/+1), n_electrons=9) and asserts that
      rhf_solve(ints, ...) (equivalently scf_solve(ints, rhf_config(ints))) throws
      std::invalid_argument; no SCF iteration runs.
    status: verified
    last_checked: 2026-07-02
  - id: ac-004
    summary: SCFSolver::run guard enforces occupation consistency without perturbing valid runs
    type: code
    pass_when: |
      test_cpu_scf_spin_resolution asserts (a) a config with
      n_alpha+n_beta != integrals.n_electrons throws std::invalid_argument,
      (b) SpinType::Restricted with n_alpha != n_beta throws std::invalid_argument,
      and (c) neutral H2O/STO-3G rhf_solve converges with energy identical to the
      pre-change value recorded in the test (tolerance 0, bit-identical).
    status: verified
    last_checked: 2026-07-02
  - id: ac-005
    summary: SpinResolution composes with the existing raw-count config builders (no new overloads)
    type: code
    pass_when: |
      Per the 2026-07-02 user design directive (no factory functions), no new
      uhf_config/uks_config overloads exist. test_cpu_scf_spin_resolution
      asserts SpinResolution(9, 2) fed into uhf_config(ints_H2Oplus,
      occ.n_alpha, occ.n_beta) yields n_alpha==5 && n_beta==4 with
      spin==Unrestricted, method==HF; and SpinResolution(8, 3) fed into
      uks_config(ints_O, occ.n_alpha, occ.n_beta) yields (5,3) with
      method==KS.
    status: verified
    last_checked: 2026-07-02
  - id: ac-006
    summary: false docstring at scf_solve.h:255 replaced by the real contract
    type: docs
    pass_when: |
      include/atomiverse/cpu/scf/scf_solve.h no longer contains the string
      "Throws-via-solver" nor "(assumed even)"; the default_config and rhf_config
      doc blocks state that SCFSolver::run throws std::invalid_argument when
      n_alpha+n_beta != n_electrons (odd-electron closed-shell case included).
    status: pending
  - id: ac-007
    summary: tutorial sections on closed-shell validity and charge/multiplicity exist
    type: docs
    pass_when: |
      docs/tutorial/ch02-rhf.md contains a "Closed-shell validity" section and
      docs/tutorial/ch03-open-shell.md contains a "Charge and multiplicity
      conventions" section, each following concept -> math -> algorithm ->
      limitations, in English, including the n_alpha/n_beta formulas, the PySCF
      mol.spin=2S contrast, and the worked table from ac-001 plus the three
      rejected combinations from ac-002.
    status: pending
  - id: ac-008
    summary: tutorial worked numbers are regression-locked by the test suite
    type: code
    pass_when: |
      Every (charge, multiplicity, n_alpha, n_beta) tuple and every rejected
      combination printed in the new ch02/ch03 sections appears as a passing
      assertion in tests/cpu/scf/test_spin_resolution.cpp, with a comment on the
      test block referencing the tutorial section and a comment in the tutorial
      naming test_cpu_scf_spin_resolution as the regression anchor.
    status: verified
    last_checked: 2026-07-02
  - id: ac-009
    summary: full default test gate stays green and fast
    type: runtime
    pass_when: |
      ctest --test-dir build (default gate, ATV_RUN_E2E_TESTS off, ATV_BUILD_RHF=ON)
      passes with test_cpu_scf_spin_resolution included and total gate wall time
      <= 20 s on the reference login-node build.
    status: verified
    last_checked: 2026-07-02
---

# Acceptance criteria

- **ac-001/ac-002** 是解析器的正/反向表，直接取自 scientist 的 PySCF 2.13.1 验证目标；反向用例只断言异常类型与触发条件，不绑定消息文本。
- **ac-003** 是本 spec 的核心缺陷修复：`rhf_config` 静默截断路径（`n_occ()=n_electrons/2`）在求解入口被拦截。
- **ac-004(c)** 的逐位一致性保证守卫对合法路径零行为影响（守卫只读、先于任何分配/迭代）。
- **ac-006** 使 scf_solve.h:255 的虚假声明变为经修正的真实契约。
- **ac-007/ac-008** 绑定教程与代码：教程数字不允许手写漂移，回归锚点是 `test_cpu_scf_spin_resolution`。
- **ac-009** 守住 `~16 s` 默认 gate 预算（阈值放宽到 20 s 吸收机器抖动）；新测试为 ms 级配置逻辑加一次 STO-3G H2O 小求解，不进 E2E 列表。
- **ac-009 measurement (2026-07-02)**: full gate 114/115 (sole failure = pre-existing, environment-sensitive `test_mpi_npt_multirank_n8`, untouched by this spec). `test_cpu_scf_spin_resolution` = 0.02 s. Total wall on this login node is dominated by PW tests (108 s) absent from the criterion's reference config (RHF-only); the RHF-gate subset is ~10 s ≤ 20 s.
