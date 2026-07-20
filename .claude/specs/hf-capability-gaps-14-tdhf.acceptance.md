---
slug: hf-capability-gaps-14-tdhf
criteria:
  - id: ac-001
    summary: TDHFPack 约化 RPA 与完整非对称 RPA 及稠密 TDA 一致
    type: code
    pass_when: |
      tests/cpu/scf/test_tdhf_pack.cpp 通过:H₂/STO-3G 上 TDA 根 == A 本征值;
      约化 (A+B)(A−B) 形式的 RPA 根与完整 [[A,B],[−B,−A]] 稠密本征值一致到
      1e-10 Ha;Davidson 最低根 == 稠密到 1e-8;每个最低根 ω_RPA ≤ ω_TDA。
    status: pending
  - id: ac-002
    summary: TDHFPack 为构造接线的 internal pack,matrix-free σ 复用 OrbitalHessian
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/internal/tdhf_pack.h 存在;TDHFPack 由构造函数
      接 (OrbitalHessian&, Channel);grep 无新增 build_*/make_* 自由函数,
      no_factory_scrub.sh 通过;公共头不含 IntegralPack/SCFConfig 引用。
    status: pending
  - id: ac-003
    summary: triplet-instability 陷阱:tdhf 前置稳定性检查 throw,tda 保持实数
    type: code
    pass_when: |
      tests/cpu/scf/test_tdhf_method.cpp 通过:拉伸 H₂ (R=1.5 Å, STO-3G) 上
      RHF::tdhf(…, triplet) 抛异常且消息含稳定性诊断;同一参考上 RHF::tda()
      返回全实激发能;合成不稳定 kernel 上 TDHFPack::rpa() 亦 throw。
    status: pending
  - id: ac-004
    summary: 公共薄表面:eV 激发能、跃迁偶极、振子强度、双通道、未 run throw
    type: code
    pass_when: |
      ExcitedStatesResult 含 energy_ev / transition_dipole_au /
      oscillator_strength / channel;未 run 先调 tda()/tdhf() throw;
      singlet 偶极经链环 03 内核装配、triplet 通道 μ=0 且 f=0;
      f == (2/3)·ω·|μ|² (a.u.) 到 1e-12(单测断言)。
    status: pending
  - id: ac-005
    summary: 甲醛回归:TDA 与 TDHF 激发能 vs 钉定 PySCF 至 1e-4 eV
    type: scientific
    pass_when: |
      regressions/scf/ch2o_tdhf 以库外消费者身份跑 HF/6-31G* 甲醛,TDA 与
      TDHF 各前 3 个 singlet 根与 REFERENCE.md 钉定的 PySCF tdscf 值逐根差
      < 1e-4 eV,n→π* 归属一致;运行 ≤ 60 s;REFERENCE.md 记录 PySCF
      版本/输入并引用 Hirata 1999 与 Dreuw 2005 DOI。
    status: pending
  - id: ac-006
    summary: 亮态振子强度与 PySCF 交叉核对
    type: scientific
    pass_when: |
      同一回归中 π→π* 亮态的长度规范振子强度与 REFERENCE.md 钉定 PySCF 值
      之差 < 1e-3(无量纲),且 n→π* 根 f ≈ 0(< 1e-3)。
    status: pending
  - id: ac-007
    summary: ch20 教程"激发态:TDHF/RPA 与 TDA"节完整
    type: docs
    pass_when: |
      docs/tutorial/ch20-hf-properties-environment.md 含该节:极化率极点叙事、
      A/B kernel 第四配置、TDA-vs-RPA 与 triplet-instability 回链 ch03、
      甲醛 hands-on(molvis 可视化)、章末 PySCF 激发能标量对照表(无条形图)。
    status: pending
  - id: ac-008
    summary: 全套构建与测试绿灯
    type: runtime
    pass_when: |
      -DATV_BUILD_HF=ON 配置下 cmake --build 无错;ctest 全绿(tests/ 秒级,
      新增两单测已注册);regressions 目标含 ch2o_tdhf 且通过。
    status: pending
---

# Acceptance criteria

ac-001/002 钉住第四配置的数学与架构形状(约化 RPA 正确性、pack 无工厂);
ac-003 是本链环的领域陷阱(RPA 虚根 ⇔ triplet 不稳定,TDA 免疫);ac-004 钉公共
表面契约;ac-005/006 是甲醛 oracle 的科学验收;ac-007/008 收尾教程与门禁。
