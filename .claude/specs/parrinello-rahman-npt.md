---
slug: parrinello-rahman-npt
status: approved
created: 2026-07-23
track: A-cell-md
grilled: false
---

# parrinello-rahman-npt

## Summary

在 CPU MD 上交付 **Parrinello–Rahman 全胞矩阵 NPT**（可倾斜/剪切），作为与对角-only `NoseHooverNPT` 并列的 `Integrator`（或明确扩展 MTK 到完整 \(\mathbf{h}\) 动力学）。公开 API 可驱动立方→正交→三斜平衡；回归钉住各向同性极限与 MTK 对角路径一致、以及非零剪切自由度可演化。

## Domain basis

Parrinello–Rahman (1981)：胞矩阵 \(\mathbf{h}\) 为动力学变量，拉格朗日含 \(W\,\mathrm{Tr}(\dot{\mathbf{h}}^T\dot{\mathbf{h}})\) 与应力约束。现代实现常与 MTK 扩展系统结合（Martyna et al.；Shinoda et al. 2004）。压力张量 \(\boldsymbol{\sigma}\) 与现有 `VirialStress` / LAMMPS 符号约定一致（正迹 = 压缩）。单位：Å、fs、eV、eV/Å³（压力）。

参考文献：Parrinello & Rahman, *J. Appl. Phys.* **52**, 7182 (1981), DOI 10.1063/1.328693；Martyna, Tobias, Klein, *J. Chem. Phys.* **101**, 4177 (1994), DOI 10.1063/1.467468。

## Design

- **类型**：`cpu::md::ParrinelloRahmanNPT`（`Integrator` kind），构造函数参数：\(T\) [K]、\(P\) 或应力目标、\(\tau_T\)、\(\tau_P\)、可选胞质量 \(W\)。**禁止** `make_pr_npt` factory；**禁止** 公共 `*Config` bag。
- **与 MTK 关系**：对角各向同性极限应数值对齐现有 `NoseHooverNPT`（同 \(T,P,\tau\)，立方盒）；实现可共享 NHC 链代码（generalize private helper）。
- **Phase**：沿用 `initial_integrate` / `final_integrate` / `end_of_step` 绑定模式（见 `nose_hoover_npt.h`）。
- **压力**：读取 `frame` / context 已有系统应力；需要时 auto 依赖 `VirialStress`（与现 barostat 一致）。
- **NPH**：本 spec **不做** 无 thermostat 路径（见 `nph-ensemble`）。

### Reuse decision

- `reuse` `NoseHooverNPT` Trotter/NHC 模式与压力符号约定。
- `reuse` `VirialStress` / `requires_system_stress` 注入路径。
- `reuse` `cpu::Box::set_cell` / 完整胞矩阵 API（triclinic cell relax 已用）。
- `generalize` NHC 链更新若与 MTK 重复 → private shared helper under `cpu/dynamics/`。
- `pattern` `Barostat` / `NoseHooverNPT` 成员参数 + Doxygen 单位。

## Files to create or modify

- `include/atomiverse/cpu/dynamics/parrinello_rahman_npt.h` (new)
- `src/cpu/dynamics/parrinello_rahman_npt.cpp` (new)
- `include/atomiverse/cpu/dynamics/` 导出/聚合头（若有）
- `tests/cpu/dynamics/test_parrinello_rahman_npt.cpp` (new)
- `regressions/md/pr_npt_lj/reg_pr_npt_lj.cpp` (new) + `REFERENCE.md` (new)
- CMake test/regression 注册文件（与现有 MD 目标一致）

## Tasks

- [ ] Write failing unit tests for PR cell update and isotropic limit (tests/cpu/dynamics/test_parrinello_rahman_npt.cpp)
- [ ] Implement cpu::md::ParrinelloRahmanNPT Integrator (header + src)
- [ ] Generalize shared NHC helper with NoseHooverNPT if duplication appears
- [ ] Document Doxygen parameters with units on public members
- [ ] Add regression regressions/md/pr_npt_lj (LJ fluid; hard-coded mean P or cell stats window)
- [ ] Verify isotropic PR vs MTK diagonal consistency on cubic box
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：单步胞更新有限、体积正、min-image 缓存刷新。
- **Unit**：立方各向同性目标下 off-diagonal \(\dot h\) 噪声有界。
- **Regression**：小 LJ NPT；硬编码平均压力/密度窗口（可对 MTK 同设定对照写入字面量）；≤ 60 s。
- 无 e2e 在 `tests/`。

## Out of scope

- NPH（无热浴）
- CUDA PR
- 表面张力半各向同性专用耦合协议（可用全矩阵表达则文档说明即可）
- Monte Carlo 体积移动

---

**Librarian note：** `cpu/dynamics/`；mirror NoseHooverNPT；confidence high。
