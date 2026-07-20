---
slug: hf-capability-gaps-06-ecp
criteria:
  - id: ac-001
    summary: ECP 表结构 + 查表类单元测试通过(钉 PySCF)
    type: code
    pass_when: |
      tests/cpu/scf/test_ecp.cpp 注册进 ctest 并通过:ecp_def2() 中 I/Au 的
      n_core、通道数与首末指数等于测试内钉住的 PySCF def2-ECP 值;
      EcpPotential 对表外元素与未知名字抛含支持项的 std::runtime_error;
      测试不跑任何完整 SCF,运行秒级。
    status: pending
  - id: ac-002
    summary: U 矩阵元(含高 l 投影道)钉 PySCF ECPscalar
    type: code
    pass_when: |
      EcpIntegrals 产出的若干 ⟨μ|U|ν⟩ 矩阵元(至少含一个 f 投影道高 l 项)
      与测试内钉住的 PySCF mol.intor('ECPscalar') 值差 ≤ 1e-10 Ha;
      实现为类成员,tests/architecture 与 04-methods 通配门下
      cpu/scf 无新增 build_* 自由函数。
    status: pending
  - id: ac-003
    summary: 价电子计数与 Z_val 接线正确,S/T/ERI 不动
    type: code
    pass_when: |
      单元测试断言:I₂ 的 IntegralPack(ecp="def2-ECP") n_electrons == 50、
      V_nuc 按 Z_val=25 计;同一几何 ecp 为空时 S/T 矩阵与现状逐元素相等;
      轻元素(全在表外)分子在 ecp="def2-ECP" 下 Hcore 与全电子路径一致。
    status: pending
  - id: ac-004
    summary: ECP↔x2c 逐元素互斥 THROW 门生效
    type: code
    pass_when: |
      设置 ecp="def2-ECP" 且 x2c=true 时,对含 ECP 元素(I)的 Frame 调 run()
      抛含元素名的 std::runtime_error;对全电子分子(H2O)同设置不抛
      (单元测试两向断言)。
    status: pending
  - id: ac-005
    summary: 生成器字节可复现(ecp_def2.h + def2_svp.h 重元素行)
    type: code
    pass_when: |
      按两表头 Regenerate 注释中的命令重跑 tools/gen_basis_data.py 后
      `git diff --stat include/atomiverse/cpu/scf/basis_data/` 为空,
      且 def2_svp.h 的 Z=37–86 行非空(test_ecp 抽查 I/Au 行)。
    status: pending
  - id: ac-006
    summary: I₂ 与 AuH RHF 能量对账 PySCF oracle
    type: scientific
    pass_when: |
      reg_ecp_heavy_atom 中 RHF/def2-SVP+def2-ECP 的 I₂ 与 AuH 总能量与
      源码注释里钉住的 PySCF(ecp='def2-svp',注明版本与命令)值差
      ≤ 1e-6 Ha,退出码 0,墙钟 ≤ 60 s,仅经公开 API(RHF 成员)驱动。
    status: pending
  - id: ac-007
    summary: 教程第 4 章 "有效核势" 一节落地
    type: docs
    pass_when: |
      docs/tutorial/ch04-basis-sets.md 含 "有效核势" 一节:价电子近似动机、
      半定域形式 U_L + Σ(U_l−U_L)|lm⟩⟨lm|、价电子计数说明、I₂ 上手代码,
      以及章末含 I₂/AuH 行的 Atomiverse-vs-PySCF 标量对照表(无柱状图)。
    status: pending
  - id: ac-008
    summary: 全量构建 + 单元测试套件绿
    type: runtime
    pass_when: |
      -DATV_BUILD_HF=ON 配置下 cmake --build 无错误,ctest 全绿且 tests/
      仍为秒级单元层;回归只挂 regressions target,不进默认 ctest。
    status: pending
---

# Acceptance criteria

ac-001/002/003/004 由 `/mol:impl` 跑 `test_ecp` 验证(PySCF 钉值在实现时
外部运行一次并写进测试源码注释);ac-005 由实现者按 Regenerate 命令重跑
生成器后查 `git diff`;ac-006 经 regressions 目标验证;ac-007 人工阅读;
ac-008 为收尾门。
