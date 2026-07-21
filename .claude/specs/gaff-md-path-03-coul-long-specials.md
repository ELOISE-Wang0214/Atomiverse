---
title: CoulLong real-space special_bonds — PME 半边也认 Amber
slug: gaff-md-path-03-coul-long-specials
status: approved
chain: gaff-md-path
link: 03/07
created: 2026-07-21
depends_on: gaff-md-path-02-pair-specials
---

# CoulLong real-space special_bonds — PME 半边也认 Amber

## Summary

为 `cpu::CoulLong`（Ewald/PME 实空间半边）增加与 `PairCoulCut` 同语义的 `set_special_bonds(s12,s13,s14)`，使 GAFF 在 `lj/cut/coul/long` + kspace 配置下 1–2/1–3 不进 real-space `erfc(αr)/r`、1–4 半权。缺此能力时，键合近邻会被长程静电实空间重复计数，聚合物熔体能量不可用。

## Domain basis

- Ewald 劈分：\(E_{\mathrm{coul}} = E_{\mathrm{real}}[\mathrm{erfc}(\alpha r)/r] + E_{\mathrm{recip}} + E_{\mathrm{self}}\)。LAMMPS `pair_style lj/cut/coul/long` 的 special_bonds 作用于 **real-space pair 循环**；reciprocal + self 仍对全部电荷求和（与 LAMMPS 一致）。
- Amber：\(s_{12}=s_{13}=0\)，\(s_{14}=0.5\) 应用于 real-space Coulomb。
- 现状：`coul_long.cpp` **无任何 special 处理**；所有邻居对全权 `erfc`。
- 参考实现语义对齐 02 的 PairCoulCut；α 与 cutoff 仍由构造参数给定，并与 `KSpacePME` 共享 α。
- 文献：Ewald, *Ann. Phys.* **369**, 253 (1921); Essmann *et al.*, *J. Chem. Phys.* **103**, 8577 (1995) (SPME).

## Design

- **API**：`CoulLong& set_special_bonds(double s12, double s13, double s14) &` 与 `&&` 重载，镜像 `PairCoulCut`。
- **实现**：复用 02 的 tier 缩放策略（抽取共享 helper 优先，若 02 已落地内部函数则直接调用；否则在 CoulLong 内平行实现，禁止行为漂移）。
- **`LongRangeParticipant` / neighbor_cutoff**：不变。
- **回归**：`reg_npt_pme_nacl` 为离子晶体、无 bonded specials — 默认 s=1 行为必须逐测例保持；不改 α 启发式。
- **CUDA `CoulLong`**：本环可选镜像；若时间紧可 CPU-only 并在 Out of scope 写明。默认 **CPU 必须**，CUDA 对齐列为加分但不进 acceptance 硬门。

## Files to create or modify

- `include/atomiverse/cpu/potentials/coul_long.h`
- `src/cpu/potentials/coul_long.cpp`
- `tests/cpu/potentials/test_coul_long_special.cpp` (new)
- `tests/CMakeLists.txt`
- 若抽取 helper：`include/atomiverse/cpu/potentials/special_scale.h` (new) + 对应 cpp（仅当 02 未抽时在本环抽）

## Tasks

- [ ] Write failing tests in `tests/cpu/potentials/test_coul_long_special.cpp`：有拓扑两电荷 + CoulLong（可配合大盒/小 α 使 recip 可忽略或固定 KSpace 关闭的 pure real 配置）验证 s14=0.5 半权；s12=0 时 1–2 对 real 能量为 0
- [ ] Implement `set_special_bonds` + compute-path scaling in `coul_long.h` / `coul_long.cpp`
- [ ] Share or mirror scaling helper with PairCoulCut so Amber weights cannot diverge
- [ ] Verify `reg_npt_pme_nacl` and existing ewald/pme unit tests remain green (no specials set)
- [ ] Register test; Doxygen units and Amber example on CoulLong
- [ ] Run full check + test suite

## Testing strategy

- **快乐路径**：1–4 对、固定构型、`CoulLong` only（或 real+固定 recip）：\(E_{\mathrm{real}}(0.5)=\tfrac12 E_{\mathrm{real}}(1)\) 到 1e-10 相对（erfc 路径）。
- **边界**：无 special 列表 ≡ 全权；`(0,0,0)` 排除 1–2/1–3/1–4。
- **领域守卫**：NaCl PME 回归不设 specials，能量门槛不变。

## Out of scope

- Reciprocal 空间 special 修正（LAMMPS 也不对 recip 做 1–4 scale）。
- 完整 GAFF 单点能量（→ 05）。
- CUDA CoulLong specials 硬门。
- `pair_style` 复合组件（LJ+Coul 一体 class）。
