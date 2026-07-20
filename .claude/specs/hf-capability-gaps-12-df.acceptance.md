---
slug: hf-capability-gaps-12-df
criteria:
  - id: ac-001
    summary: def2-JKFIT 静态表 + AuxBasis 落地,表值抽查过、缺行元素抛异常
    type: runtime
    pass_when: |
      tests/cpu/scf/test_density_fitting.cpp 中表抽查用例通过:basis_data::def2_jkfit()
      的 H/C/O 抽查壳与 Weigend 2002 发表值逐数一致;对表中无行元素构造辅助
      ShellList 抛 std::runtime_error;ctest 全绿且整套单元仍秒级。
    status: pending
  - id: ac-002
    summary: g(l=4)壳 2c/3c 积分正确,根预算未涨
    type: runtime
    pass_when: |
      test_density_fitting 中 (g|g) 与含 g 辅助壳的 3-中心积分对测试注释里
      外部钉定的 PySCF 参考值 ≤1e-10 一致;kRysMaxRoots 仍为 7(静态断言编译过)。
    status: pending
  - id: ac-003
    summary: j_fitted/k_fitted 数学正确(小体系对稠密精确 J/K 收口)
    type: runtime
    pass_when: |
      H2/STO-3G 单次调用:V^{-1/2}·V_PQ·V^{-1/2}=I ≤1e-12;‖J_fit−J_exact‖∞ 与
      ‖K_fit−K_exact‖∞ 均落在测试内文档化的拟合误差上界内;用例通过。
    status: pending
  - id: ac-004
    summary: 接缝成员驱动:无 build_* 自由函数,RI 选路是方法成员
    type: code
    pass_when: |
      tests/architecture/no_factory_scrub.sh(含 04-methods 的 cpu/scf build_* 通配门)
      通过;RI 枚举在 scf_enums.h,选路仅经 MeanFieldMethod 成员 ri/aux_basis →
      FockStage 构造参 → HamiltonianBuilder 成员工作区,树内无 DF 相关自由函数工厂。
    status: pending
  - id: ac-005
    summary: 默认 RI::None 行为位级不变
    type: runtime
    pass_when: |
      未改动的现有 SCF 单元测试(test_scf_solver、test_fock_stage 等)在本环改动后
      全部通过,无期望值改动的 diff。
    status: pending
  - id: ac-006
    summary: H2O/def2-TZVP RI-JK 双门回归:≤50 μHa vs 自家精确 ERI 且 ≤1e-6 Ha vs PySCF
    type: scientific
    pass_when: |
      regressions/scf/h2o_ri_jk 通过:|E(RI-JK) − E(exact-ERI HF)| ≤ 5e-5 Ha,且
      |E(RI-JK) − REFERENCE.md 钉定的 PySCF mf.density_fit() 值| ≤ 1e-6 Ha;
      墙钟 ≤60 s;仅经公共头 + link atomiverse。
    status: pending
  - id: ac-007
    summary: glycine/def2-SVP 速度表钉定:ours-RI / ours-exact / PySCF-DF
    type: performance
    pass_when: |
      benchmarks/scf/validation/ri_speed/ 含可复跑脚本 + README(机器/编译旗标/
      重复次数方法学);产出的三列墙钟标量表钉入 ch19 章末,且该表中
      ours-RI-JK 墙钟 < ours-exact-ERI 墙钟。
    status: pending
  - id: ac-008
    summary: ch19"密度拟合 RI"一节成文,误差预算与取舍诚实
    type: docs
    pass_when: |
      docs/tutorial/ch19-hf-large-systems.md 含:4-指标→3-指标压缩叙事、μHa 误差
      预算表、单表 JKFIT(而非 JFIT)取舍说明、失配抛异常陷阱、hands-on 代码、
      章末精度+速度表;三篇 DOI 全部引用;新公共符号有 @brief/@param(带单位)。
    status: pending
---

# Acceptance criteria

ac-001/002/003 是 RED-first 单元层(表、积分、拟合数学);ac-004/005 锁架构接缝与
零回归;ac-006 是科学双门(二阶误差论断 + 外部 oracle);ac-007 把速度证据固化为
可复跑方法学 + 教程钉数(/mol:bench 或 benchmarks 流程执行);ac-008 收教程叙事。
