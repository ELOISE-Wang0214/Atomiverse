---
title: Pair 势 special_bonds 真缩放 — Amber 1–4 = 1/2
slug: gaff-md-path-02-pair-specials
status: approved
chain: gaff-md-path
link: 02/07
created: 2026-07-21
depends_on: gaff-md-path-01-io-full
---

# Pair 势 special_bonds 真缩放 — Amber 1–4 = 1/2

## Summary

把 `cpu::PairLJCut` 与 `cpu::PairCoulCut` 的 `set_special_bonds(s12,s13,s14)` 从「仅当 scale=0 时硬排除」升级为 **LAMMPS 语义的按 tier 缩放**：1–2 / 1–3 / 1–4 对分别乘以 \(s_{12},s_{13},s_{14}\)，支持 Amber 权重 \((0,0,1/2)\)。GAFF/Amber 分子在有拓扑时，LJ 与 cut 库仑的 1–4 作用为半权，而不是全无或全算。

## Domain basis

- LAMMPS `special_bonds`：对非键 pair 循环中属于 1–2、1–3、1–4 的原子对，能量与力乘以权重；`special_bonds amber` ≡ LJ 与 Coulomb 均为 \(0,0,0.5\)（LAMMPS 文档 *special_bonds* 命令）。
- Amber/GAFF 惯例：1–2、1–3 非键完全关闭（由 bonded 项覆盖）；1–4 保留一半 LJ 与一半库仑。见 Cornell *et al.*, *J. Am. Chem. Soc.* **117**, 5179 (1995); Wang *et al.*, *J. Comput. Chem.* **25**, 1157 (2004) (GAFF)。
- 当前实现：`any_excl_ = (s12==0)||(s13==0)||(s14==0)` 且仅在 scale==0 时把对放入 exclusion 集并 **跳过**；\(s_{14}=0.5\) 时 1–4 仍全权计入 — 与 Amber 不符。
- 目标语义（与 LAMMPS 一致）：
  - scale = 0 → 跳过该对；
  - scale ∈ (0,1] → 能量与力乘 scale；
  - 默认未设置 specials 时行为与今日「全权、无 exclusion」一致。
- 单位：势能 eV，力 eV/Å；缩放无量纲。

## Design

- **共享逻辑**（优先抽到 `cpu/potentials/` 内匿名/内部 helper，或复制最小公共刷新路径，避免新顶层概念）：从 `kAtom_nspecial` / `kAtom_special` 构建 **按 tier 的 partner 表**（或单一 map + tier 标签），compute 热路径对每对查 1–2→1–3→1–4 优先级（LAMMPS：较近 tier 优先）。
- **`PairLJCut` / `PairCoulCut`**：`set_special_bonds` API 签名不变；实现改为真缩放。文档字符串写明 Amber 例子 `.set_special_bonds(0.0, 0.0, 0.5)`。
- **向后兼容**：现有测试使用 `(0,0,0)` 全排除 — 仍应全跳过，数值不变。
- **CUDA**：本环 **CPU only**。CUDA `PairLJCut` 已有部分 scale 但要求 s12==s13==s14；CUDA Amber 差异列入 Out of scope / 后续。
- 不改 neighbor list 构建；specials 只在 pair 核内缩放。

## Files to create or modify

- `include/atomiverse/cpu/potentials/pair_lj_cut.h` — Doxygen 语义
- `include/atomiverse/cpu/potentials/pair_coul_cut.h` — Doxygen 语义
- `src/cpu/potentials/pair_lj_cut.cpp` — 真缩放
- `src/cpu/potentials/pair_coul_cut.cpp` — 真缩放
- `tests/cpu/potentials/test_pair_special_scale.cpp` (new)
- `tests/CMakeLists.txt`

## Tasks

- [ ] Write failing tests in `tests/cpu/potentials/test_pair_special_scale.cpp`：两原子 1–4 拓扑下 LJ 与 CoulCut 在 s14=0.5 时能量 = 无 special 的一半；s14=0 时能量为 0；s14=1 时等于全权
- [ ] Implement tier-aware scaling in `src/cpu/potentials/pair_lj_cut.cpp`
- [ ] Implement tier-aware scaling in `src/cpu/potentials/pair_coul_cut.cpp`
- [ ] Update Doxygen on both headers with Amber (0,0,0.5) example and scale semantics
- [ ] Register test; verify existing tip3p/rigid_water exclusions with (0,0,0) still pass
- [ ] Run full check + test suite

## Testing strategy

- **快乐路径**：合成丁烷式最小拓扑（或两原子人为标为 1–4）：`s14=0.5` → \(E = \tfrac12 E_{\mathrm{full}}\)，力同比例（1e-12 相对，双精度）。
- **边界**：`s12=0,s13=0,s14=0` 全排除；无 special 列表时与旧行为一致；迁移后 refresh exclusion/scale map（既有 multi-rank 注释路径不回归）。
- **领域**：本环不跑完整 GAFF；缩放因子本身是无量纲约定，对照解析半权即可。

## Out of scope

- `CoulLong` specials（→ 03）。
- CUDA pair specials 完整 Amber 混合权重。
- 不同 lj/coul 权重（LAMMPS 可分设）；v1 与现 API 一样 **lj 与 coul 各自组件各自 set**，调用方对两者都传 Amber 三元组。
- 1–4 dihedral 内建 pair weight（CHARMM 可选）— 仍由 pair specials 表达。
