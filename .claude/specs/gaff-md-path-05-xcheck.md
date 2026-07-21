---
title: GAFF 单点能量交叉验证 — 小 PEO 对照钉值
slug: gaff-md-path-05-xcheck
status: approved
chain: gaff-md-path
link: 05/07
created: 2026-07-21
depends_on: gaff-md-path-04-ff-compose
---

# GAFF 单点能量交叉验证 — 小 PEO 对照钉值

## Summary

新增 **forcefield 回归**：在固定小体系（单链 PEO 片段或 2 链小盒，全原子 GAFF2 拓扑 + 电荷）上，用 `cpu::Evaluate` 组装 `PairLJCut`（多 type）+ `CoulLong` + `KSpacePME` + harmonic bond/angle + `DihedralCharmm` + Amber specials，输出分项能量，并与 **预钉的 LAMMPS 参考值**（写入仓库的 JSON/REFERENCE，不在 CI 调 LAMMPS）比对。这是「能跑 GAFF」的科学硬门：不是协议 30 ns，而是力场接线正确。

## Domain basis

- 能量分解（与 LAMMPS `thermo_style` 可对齐的项）：  
  \(E = E_{\mathrm{bond}} + E_{\mathrm{angle}} + E_{\mathrm{dihed}} + E_{\mathrm{LJ}} + E_{\mathrm{coul,real}} + E_{\mathrm{coul,recip}} + E_{\mathrm{coul,self}}\)。  
  ATV 经 `kEnergy*` meta / system scalars 暴露的项写入 REFERENCE 对照表。
- 容差（初始）：总分 \(E\) 相对误差 ≤ **1e-4**（或绝对 1e-3 kcal/mol 转 eV 后取严者）；分项 ≤ **5e-4** 相对（PME 网格/α 与 LAMMPS PPPM 非同一算法 — 文档记录 SPME vs PPPM 差异，α 与 cutoff 在 REFERENCE 钉死使总能量可复现）。
- 若逐位对齐困难：允许 **同一 ATV 实现的解析/FD 力检查** + 与钉值能量比对；钉值生成步骤写在 `tools/` 脚本（可选，CI 不跑 LAMMPS）。
- 夹具来源：精简自 peo25（DP 更短）或仓库内 `regressions/forcefield/gaff_peo_fragment/` 自带 `.data` + 系数表；**禁止** CI 依赖 `/home/jicli594/work/molpack-paper/...` 绝对路径。
- 文献：Wang GAFF；Gudla & Zhang, *J. Phys. Chem. B* **128** (2024) 仅作协议上下文，本环不做 anneal。

## Design

- **回归目录**：`regressions/forcefield/gaff_peo_fragment/`  
  - `reg_gaff_peo_fragment.cpp`  
  - `REFERENCE.md`（钉值、α、cutoff、kspace 设置、如何用 LAMMPS 再生）  
  - `fragment.data` + `coeffs.json`（或 cpp 内嵌表）
- **组装顺序**（Evaluate）：NeighborList → PairLJCut multi-type + specials(0,0,0.5) → CoulLong + specials → KSpacePME → Bond/Angle/Dihedral → run once → 读能量。
- **不**引入新 Driver；纯回归可执行文件。
- CMake：`regressions/CMakeLists.txt` 注册。

## Files to create or modify

- `regressions/forcefield/gaff_peo_fragment/reg_gaff_peo_fragment.cpp` (new)
- `regressions/forcefield/gaff_peo_fragment/REFERENCE.md` (new)
- `regressions/forcefield/gaff_peo_fragment/fragment.data` (new)
- `regressions/forcefield/gaff_peo_fragment/coeffs.json` (new) — 或等价内嵌
- `regressions/CMakeLists.txt`
- 可选：`tools/pin_gaff_lammps.py` (new) — 文档化钉值流程，CI 不强制

## Tasks

- [ ] Write failing regression binary skeleton + REFERENCE tolerances in `regressions/forcefield/gaff_peo_fragment/`
- [ ] Add fragment fixture (data + coeffs) vendored in-repo
- [ ] Implement Evaluate assembly (01–04 APIs) and compare total/per-term energies to pins
- [ ] Document SPME vs PPPM and pin α, real cutoff, grid/spacing in REFERENCE.md
- [ ] Register regression in `regressions/CMakeLists.txt`
- [ ] Run full check + test suite including the new regression target

## Testing strategy

- **快乐路径**：回归 exit 0；总能量在容差内。
- **边界**：缺电荷文件应在装载期失败（依赖 01）。
- **领域验证**：REFERENCE 含钉值表与再生命令；科学 criterion 绑定该回归。
- 单元测试层不重复完整 PME 盒；本环 = 科学回归。

## Out of scope

- 30 ns 协议 / 密度 / \(R_g\) 统计（→ 06 仅 smoke）。
- 在线调用本机 LAMMPS。
- CUDA 组装与 GPU 回归。
- 与 paper N=100 盒的定量重合。
