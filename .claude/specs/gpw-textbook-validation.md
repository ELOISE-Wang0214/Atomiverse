---
title: GPW validation on minimal textbook models
slug: gpw-textbook-validation
status: approved
science_required: true
created: 2026-07-11
acceptance: gpw-textbook-validation.acceptance.md
supersedes: gpw-edl-aimd (2026-06-20)
---

# GPW 验证 — 教科书最小模型

## Summary

**2026-07-11 用户指令重写。** 前身 spec（gpw-edl-aimd）以 1214 原子 ZnO 电双层
生产体系为验收 capstone。该体系单点需数小时机时（13726 基函数、2888 占据轨道），
把"修复的正确性"与"生产规模收敛"混在一个判据里 —— 正确性验证被算力绑架。指令：

> 所有的模型都要用最小的教科书的典型模型，删除之前错误的几千个原子的垃圾模型
> 和代码。所有的模型测试要在 10 分钟内完成！

本 spec 的口径：**GPW 的每一条物理判据都在教科书最小体系上证明**（H₂、Ar、H₂O、
NaCl 双原子、纤锌矿 ZnO 小胞），每个验收用例 **≤ 10 分钟**。生产规模运行是下游
使用，不是代码交付物，永不出现在验收判据里。

已删除的垃圾模型与代码：`examples/10_gpw_znoedl.cpp`（1214 原子 EDL 单点）、
`tests/` 里的 1214 原子 `final_frame.xyz` 加载用例与 `ATV_FINAL_FRAME_XYZ`
CMake 机制、代码注释中的全部 EDL 品牌引用。

**同日追加的 API 指令（均已执行）：**
- 「全局应该只有一个Element」「Element应该直接用molrs定义的」→ molrs 经 CXX 桥导出
  周期表枚举 `molrs::Element`（H=1…Og=118，与核心 Rust 枚举逐值 parity 测试），
  `ATV::Element` 为其别名；`GpwElement` 删除，改为 CP2K 语义的
  `Kind{element, basis, pot}`（枚举不可继承，用隐式转换给出 is-a 用法）。
- 「不要用任何的工厂函数，全都是OOP风格」→ `make_gpw_atom`/`builtin_gth_*`/
  `parse_gth_*`/`load_dftd3_params`/`builtin_dftd3_params` 全部删除，换成
  `GpwAtom(pos, basis, pot, z)`、`GthPotential(Element, name)`、
  `GaussianBasis(Element, name)`、`D3Params{}` 构造函数。
- 「所有的参数表都要写成.h的形式，不要动态解析」→ CP2K 语法字符串解析器
  （cp2k_io/cp2k_basis）与内嵌原始字符串删除；数据一次性生成为静态 C++ 表
  `basis_data/gpw_tables.h`（%.17g，位级往返）与 `basis_data/dftd3_tables.h`；
  `GpwConfig.d3_param_file` 字段删除（D3 只有树内表一个来源）。

## 继承自前身 spec 的已验证成果（全部保留）

- 树内 GTH-PBE / DZVP-MOLOPT / DFT-D3 数据表（Krack 2005; VandeVondele-Hutter
  2007; Grimme 2010），`ATV_CP2K_DATA_DIR` 已删，构建与测试期零第三方读取。
- 三个真 bug 修复：力对错泛函求导（17 commits 潜伏）、缺正则正交化、ZnO 对角化
  SCF 小能隙不收敛（level shift）。
- Γ 点配对晶格和缺最外层镜像 → `pair_image_counts`（每轴 +1 层），PySCF 2.13.1
  外部真值验证到机器精度（s_min −1.92 → −6.8e-14；谱与参考一致）。

## Domain basis

同前身 spec：GTH（Goedecker-Teter-Hutter PRB 54, 1703 (1996)；Krack TCA 114,
145 (2005)）；MOLOPT（VandeVondele-Hutter JCP 127, 114105 (2007)）；解析 GPW 力
五项分解 vs 中心 FD oracle（唯一力判据，无外部引擎）；CSVR（Bussi-Donadio-
Parrinello JCP 126, 014101 (2007)）；D3(BJ)（Grimme JCP 132, 154104 (2010)）。

## Testing strategy

- `tests/`（单元，ms 级）：数据表 vs 文献值、周期重叠 PSD + PySCF 钉死谱
  （单 H 原子 → 3×3 H 网格 → 镜像计数算术回归）、积分/谐函数/collocate 纯数学。
- `benchmarks/qc/gpw/validation/`（验收时跑，ctest 永不跑）：教科书体系锚点
  单点（H₂/Ar/H₂O/ZnO/NaCl）、FD 力 oracle、geoopt、H₂ NVE + H₂ CSVR-NVT 短
  AIMD、OT-vs-diag。**每个用例 ≤ 10 分钟，全套一次跑完并记录耗时。**

## Out of scope

- 任何数百/数千原子体系的运行、收敛或性能画像（生产使用，非交付物）。
- GAPW、HFX、k 点、CUDA GPW、金属展宽、通用 CP2K `.inp` 解析器。
