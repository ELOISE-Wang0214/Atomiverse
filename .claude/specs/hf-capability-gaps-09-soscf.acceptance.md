---
slug: hf-capability-gaps-09-soscf
criteria:
  - id: ac-001
    summary: Stepper 单元件全过:幺正性/梯度/信赖域
    type: code
    pass_when: |
      tests/cpu/scf/test_soscf.cpp 注册进 ctest 并通过:C·exp(κ) 保
      C^T S C = I(≤1e-12);g_ai 与微型体系旋转能量有限差分一致
      (≤1e-6 Ha);信赖域把 ‖κ‖ 截到 soscf_trust;运行秒级。
    status: pending
  - id: ac-002
    summary: CIAH 矩阵自由解 = 稠密增广矩阵参考解
    type: code
    pass_when: |
      同一测试内,小维度问题上 SOSCFStepper 的 (λ, κ) 与测试里显式构造的
      增广矩阵 [[0,g^T],[g,H]] 稠密本征解一致(κ 方向余弦 ≥ 1-1e-8);
      生产代码路径中 H 从不显式构造(σ-vector 调 08 环内核,代码审查确认)。
    status: pending
  - id: ac-003
    summary: 升级成员接线,默认关断路径逐位不变
    type: code
    pass_when: |
      RHF/UHF 经继承获得 soscf/soscf_start/soscf_max_micro/soscf_trust
      成员(默认 false/0.25/12/0.4,Doxygen 带单位);soscf=false 时
      H₂/STO-3G 微型 SCF 能量与迭代轨迹和本环之前逐位相同
      (test_soscf 断言),既有 scf 单元测试全绿。
    status: pending
  - id: ac-004
    summary: 难例:CDIIS 振荡而 SOSCF 收敛,能量对 PySCF newton
    type: scientific
    pass_when: |
      reg_soscf_hard_case 中同一难例、同 max_iter 预算:CDIIS-alone 不收敛
      (断言 converged==false 或能量振荡),soscf=true 收敛且终能量与源码
      注释钉住的 PySCF mf.newton() 值(注明版本与命令)差 ≤ 1e-8 Ha;
      退出码 0,墙钟 ≤ 60 s。
    status: pending
  - id: ac-005
    summary: 迭代数与 PySCF newton 钉值可比
    type: scientific
    pass_when: |
      同一 regression 断言 SOSCF 宏迭代数 ≤ 2× 源码注释里钉住的
      PySCF newton 宏迭代数(RHF 与 UHF 难例各一行钉值)。
    status: pending
  - id: ac-006
    summary: 教程 ch18 二阶方法半章落地
    type: docs
    pass_when: |
      docs/tutorial/ch18-scf-techniques.md 含 "二阶方法" 一节:轨道旋转
      Newton-Raphson 推导、二次收敛论证、DIIS-早期/SOSCF-收尾配方叙事衔接
      07 环上半章;molplot DIIS vs SOSCF 收敛曲线;章末 PySCF 对照表行
      (迭代数 + 终能量,标量只进表)。
    status: pending
  - id: ac-007
    summary: 全量构建 + 单元测试套件绿
    type: runtime
    pass_when: |
      -DATV_BUILD_HF=ON 配置下 cmake --build 无错误,ctest 全绿且 tests/
      仍为秒级单元层(难例回归只挂 regressions target,不进 ctest)。
    status: pending
---

# Acceptance criteria

ac-001/002/003 由 `/mol:impl` 跑单元测试 + 代码审查验证(ac-002 的"H 不显式
构造"看 internal/soscf.h 只经 08 环 σ-vector 接口触 H);ac-004/005 经
regressions 目标验证,oracle 在实现时用外部 PySCF(mf.newton())钉死并写进
回归源码注释;ac-006 人工阅读 + molplot 产物存在;ac-007 为收尾门。
