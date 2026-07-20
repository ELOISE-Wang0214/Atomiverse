---
slug: hf-capability-gaps-05-basis
criteria:
  - id: ac-001
    summary: 新表结构单元测试通过(壳数/AO 数钉 PySCF)
    type: code
    pass_when: |
      tests/cpu/scf/test_basis_tables.cpp 注册进 ctest 并通过:pople_631pgs()/
      aug_ccpvdz() 逐元素壳数、l 序列、AO 计数与首末指数抽查等于测试内钉住的
      PySCF mol._basis 值;测试不跑任何 SCF,运行秒级。
    status: pending
  - id: ac-002
    summary: 四套既有表补齐 Z=11–18 且 sto3g 不动
    type: code
    pass_when: |
      pople_631g/pople_631gs/def2_svp/def2_tzvp 四表 Z=11..18 行 shells 非空
      (test_basis_tables 断言);git diff 显示 sto3g.h 零改动。
    status: pending
  - id: ac-003
    summary: 生成器字节可复现
    type: code
    pass_when: |
      按各表头 Regenerate 注释中的命令重跑 tools/gen_basis_data.py 后
      `git diff --stat include/atomiverse/cpu/scf/basis_data/` 为空
      (六个生成头全部逐字节复现)。
    status: pending
  - id: ac-004
    summary: 枚举/分派/旗标注册完整,不支持 Z 抛错
    type: code
    pass_when: |
      BasisSet::Pople631pGs 与 BasisSet::AugCcPVDZ 可经 ShellList 构造;
      basis_set_name 返回 "6-31+G*"/"aug-cc-pVDZ";02 环旗标登记为
      cartesian-6d/spherical;对表外 Z 构造抛含基组名的 std::runtime_error
      (单元测试断言,无 default 回退分支)。
    status: pending
  - id: ac-005
    summary: 阴离子/弥散能量对账 PySCF oracle
    type: scientific
    pass_when: |
      reg_diffuse_anion 中 Cl⁻/6-31+G* 与 H2O/aug-cc-pVDZ 的 RHF 总能量与
      源码注释里钉住的 PySCF 值(注明版本与命令)差 ≤ 1e-7 Ha,退出码 0,
      墙钟 ≤ 60 s。
    status: pending
  - id: ac-006
    summary: OH⁻ 弥散束缚性翻转复现
    type: scientific
    pass_when: |
      同一 regression 内 OH⁻ 的 ε_HOMO 在 6-31G* 与 6-31+G* 下均与钉住的
      PySCF 同基组值一致(≤1e-6 Ha),且 6-31+G* 值严格更低并与 PySCF
      钉值同号(弥散恢复束缚)。
    status: pending
  - id: ac-007
    summary: 教程第 4 章弥散/第三周期扩写落地
    type: docs
    pass_when: |
      docs/tutorial/ch04-basis-sets.md 含 "极化、弥散与第三周期" 一节:
      基组命名法解码、OH⁻ 弥散演示、章末含 6-31+G*/aug-cc-pVDZ 行的
      Atomiverse-vs-PySCF 标量对照表,以及 molplot 基组收敛曲线(曲线可图,
      标量只进表)。
    status: pending
  - id: ac-008
    summary: 全量构建 + 单元测试套件绿
    type: runtime
    pass_when: |
      -DATV_BUILD_HF=ON 配置下 cmake --build 无错误,ctest 全绿且 tests/
      仍为秒级单元层(regression 只挂 regressions target,不进 ctest)。
    status: pending
---

# Acceptance criteria

ac-001/002/004 由 `/mol:impl` 跑单元测试与 grep 验证;ac-003 由实现者按
Regenerate 命令重跑生成器后查 `git diff`;ac-005/006 经 regressions 目标验证
(oracle 数值在实现时用外部 PySCF 钉死并写进回归源码注释);ac-007 人工阅读 +
molplot 产物存在;ac-008 为收尾门。
