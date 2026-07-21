---
title: GAFF 系数装配 — real 单位 → eV 与 DihedralCharmm 映射
slug: gaff-md-path-04-ff-compose
status: approved
chain: gaff-md-path
link: 04/07
created: 2026-07-21
depends_on: gaff-md-path-03-coul-long-specials
---

# GAFF 系数装配 — real 单位 → eV 与 DihedralCharmm 映射

## Summary

在 `cpu/io` 增加 **数据型** 装配辅助：把 LAMMPS real 单位下的 GAFF 风格系数（kcal/mol、Å、度）转为 Atomiverse 工程单位（eV、Å、rad），并为 PEO 所用的 **单 term `dihedral_style fourier`** 建立到 `DihedralCharmm` 的显式映射。产出可直接喂给既有 `PairLJCut` / `BondHarmonic` / `AngleHarmonic` / `DihedralCharmm` / `CoulLong` 构造函数的向量与 type 表 — 不引入新 Driver，不自动 `add_potential`。

## Domain basis

- LAMMPS `units real`：能量 kcal/mol，距离 Å，时间 fs，电荷 e。Atomiverse：能量 eV。转换 \(1\,\mathrm{kcal/mol} = 1/23.060547830619026\,\mathrm{eV}\)（与 `tools/debug_amorphous` 钉值一致；本环将该常量提升到 `core/units.h` 或 `cpu/io` 可见常量，禁止再散落 magic number）。
- Harmonic bond（LAMMPS）：\(E = K(r-r_0)^2\)，\(K\) 单位 energy/length²；Angle：\(E = K(\theta-\theta_0)^2\)，\(\theta\) 在 ATV 内用 **rad**。
- Fourier dihedral（LAMMPS `dihedral_style fourier`）：  
  \(E = \sum_{i=1}^{m} K_i\bigl(1+\cos(n_i\phi - d_i)\bigr)\)，\(d\) 为度。  
  PEO GAFF 夹具每 type 均为 **m=1** → 恰为 CHARMM 形式 \(K(1+\cos(n\phi-\delta))\)，与现有 `DihedralCharmm` 一致（δ 转 rad）。
- Lorentz–Berthelot：`PairLJCut` 已对对角 type 表做 arithmetic σ / geometric ε；GAFF `pair_modify mix arithmetic` 与此一致。
- 参考：Wang *et al.*, *J. Comput. Chem.* **25**, 1157 (2004); LAMMPS `dihedral_style fourier` 文档。

## Design

- **`cpu::io` 辅助（名称建议 `gaff_units.h` 或 `lammps_real_units.h`）**：
  - `kcal_per_mol_to_ev(double)`
  - bond/angle/dihedral/LJ 批量转换 free functions（输入 real，输出 eV/rad）。
- **`cpu::ForceField` 填充（可选薄封装）**：从手写表或测试内嵌 struct 填 `def_bondstyle("harmonic")` 等；**不**要求本环解析完整 `peo25.ff` 文本（完整 `.ff` 行解析若简单可做，但 acceptance 只钉 **结构化表 → 系数向量**）。
- **`FourierTerm` → Charmm**：`struct { double K_kcal; int n; double d_deg; }` → `(K_ev, n, delta_rad)`；m>1 时本环 **抛异常或返回错误**（显式 Out of scope 多 term 求和组件），避免静默丢 term。
- **特殊键默认**：装配文档/辅助返回推荐 Amber `(0,0,0.5)` 常量，由调用方 `set_special_bonds`。
- 单测用 PEO 子集：2–3 bond types、1 dihedral type（取自 peo25.ff 钉值）。

## Files to create or modify

- `include/atomiverse/core/units.h` — 增加 kcal/mol ↔ eV（或旁路 `cpu/io/lammps_real.h`）
- `include/atomiverse/cpu/io/lammps_real.h` (new) — 转换 + fourier→charmm
- `src/cpu/io/lammps_real.cpp` (new) — 若非 header-only
- `tests/cpu/io/test_lammps_real.cpp` (new)
- `tests/CMakeLists.txt`
- 若改 units.h：确认 amorphous 工具可后续改用（本环不强制改 tools）

## Tasks

- [ ] Write failing tests in `tests/cpu/io/test_lammps_real.cpp`：钉住 kcal→eV 因子、bond K 转换、angle θ 度→rad、单 term fourier→Charmm 向量；m>1 拒绝
- [ ] Implement conversion helpers in `include/atomiverse/cpu/io/lammps_real.h` (+ cpp if needed)
- [ ] Promote kcal/mol factor into `units.h` or document single source in `lammps_real.h` with constexpr and test pin 23.060547830619026
- [ ] Add Doxygen with units on every public converter
- [ ] Register test; Run full check + test suite

## Testing strategy

- **快乐路径**：`0.1078 kcal/mol` ε → eV 与手算一致（相对 1e-15）；`dihedral_coeff … 1 K n d` → Charmm 三向量长度 1。
- **边界**：m=0 或 n<1 抛 `invalid_argument`；空 type 表。
- **领域**：不跑 Evaluate；只保证单位与函数形式映射正确，05 负责能量对照。

## Out of scope

- 完整 `peo25.ff` 文本 parser / label-based type 名解析（可用 int type 索引表）。
- 多 term Fourier 专用 Component（需要时另开 spec）。
- Improper（PEO 本例为 0）。
- 自动安装到 `md::MD`。
