---
title: CUDA MD Review Cleanup — link 08 — file placement + conventions + alias header removal
slug: cuda-md-review-cleanup-08-placement-and-conventions
status: draft
created: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 10 + Cluster 11, 7 findings)
chain: cuda-md-review-cleanup-* (link 8/9)
---

# CUDA MD Review Cleanup — link 08 — file placement + conventions + alias header removal

## Summary

机械式清理(放在 chain 末尾以避免与 01-07 的 diff merge 冲突):
- 重定位 `src/cuda/md/neighborlist/` 整目录到 `src/cuda/neighborlist/`(私有 impl 应跟随公共组件 `cuda/neighbor_list.cuh` 在 cuda root)
- 重定位 `include/atomiverse/cuda/observer/stress.{cuh}` 和 `src/cuda/observer/stress.cu` 到 `cuda/md/virial_stress.{cuh,cu}`(CPU 同胞类是 `cpu/md/virial_stress.h`,CLAUDE.md 把 VirialStress 列在 MD 组件)
- 修 `angle_harmonic.cuh:18` docstring `[degrees]→[rad]`(MEMORY `[[atv-cuda-angle-theta0-double-deg2rad]]` 跟踪)
- 删除 4 个 backward-compat alias header(`lj.cuh`、`harmonic_bond.cuh`、`harmonic_angle.cuh`、`multi_harmonic_dihedral.cuh`)
- 跟踪多 GPU halo-exchange CUDA 端缺失到 `.claude/NOTES.md`(实现不在本 spec)
- 重跑 `/mol:map` 刷新 `.claude/architecture.md` 中 stale 的 `cuda::Frame embeds SCFWorkspace` 描述
- 确认 `cuda::Thermo` 命名空间模式有意(`ATV::cuda::` 而非 `ATV::cuda::recorder::`,与 CPU 一致)— 无代码改动,只是 documenting

## Domain basis

**目录归置规则**(`CLAUDE.md ## Architecture § Backend directory layering`):"Files at the backend root are **shared infrastructure** used by multiple domains. Domain subdirectories hold domain-specific code."

`cuda::NeighborList` 的公共组件在 `include/atomiverse/cuda/neighbor_list.cuh`(root,跨 md/observer/recorder 使用)。私有 impl(`nlist_impl.cuh` 等)应跟随,放在 `src/cuda/neighborlist/`,不应放在 `src/cuda/md/neighborlist/`(MD 域)。

`VirialStress` 的 CPU 同胞类在 `cpu/md/virial_stress.h`;CLAUDE.md `## Implemented components § MD components` 明确列 VirialStress。CUDA 端放在 `cuda/observer/` 是误归。

**Backward-compat alias headers**:`CLAUDE.md ## Architecture § Rules` 第二条:"No fallback values"扩展为"No `[[deprecated]]` shims"(notes 中 2026-05-19 LAMMPS-命名 refactor 是 no-backwards-compat refactor)。4 个 alias header 应该早就删除。

**multi-GPU halo-exchange**:CPU 端 `md.h` 在 `#ifdef ATV_BUILD_MPI` 时自动注入 halo;CUDA 端无对应。本 spec out-of-scope,但 NOTES 应该跟踪,避免下次审计再发现。

## Design

7 项发现见 inventory cluster 10 + cluster 11。增量:

1. **目录重定位**:
   - `mkdir src/cuda/neighborlist && git mv src/cuda/md/neighborlist/* src/cuda/neighborlist/`
   - 更新所有 `#include "md/neighborlist/..."` 改成 `#include "neighborlist/..."`(在 `src/cuda/neighbor_list.cu` 和其他 src 文件中)
   - 更新 `CMakeLists.txt` / `src/CMakeLists.txt` 中相关 glob 或显式列表
2. **VirialStress 重定位**:
   - `git mv include/atomiverse/cuda/observer/stress.cuh include/atomiverse/cuda/md/virial_stress.cuh`
   - `git mv src/cuda/observer/stress.cu src/cuda/md/virial_stress.cu`
   - 全仓库更新 `#include "cuda/observer/stress.cuh"` → `"cuda/md/virial_stress.cuh"`
   - 更新 `CMakeLists.txt`
3. **docstring 修复**:`angle_harmonic.cuh:18` `@param theta0 [degrees] (...)` → `@param theta0 [rad] (API contract — already converted by caller)`
4. **删 alias headers**:`rm include/atomiverse/cuda/potentials/{lj.cuh,harmonic_bond.cuh,harmonic_angle.cuh,multi_harmonic_dihedral.cuh}`;全仓库 grep 这些 include 改成 canonical 头
5. **NOTES 跟踪**:在 `.claude/NOTES.md` 加一段"Multi-GPU halo exchange CUDA backend gap"——描述 CPU 自动注入而 CUDA 不支持,要求未来 multi-GPU spec 处理
6. **`/mol:map`**:运行 `/mol:map`(交互式,需用户确认)刷新 `.claude/architecture.md` 中 `cuda::Frame embeds SCFWorkspace` 这一段
7. **`cuda::Thermo` 命名空间确认**:无代码改动,只在本 spec 标注"模式有意"

## Files to create or modify

- `src/cuda/md/neighborlist/` → `src/cuda/neighborlist/`(整目录 `git mv`)
- `include/atomiverse/cuda/observer/stress.cuh` → `cuda/md/virial_stress.cuh`
- `src/cuda/observer/stress.cu` → `src/cuda/md/virial_stress.cu`
- `include/atomiverse/cuda/potentials/angle_harmonic.cuh:18` — docstring
- 删:`lj.cuh` / `harmonic_bond.cuh` / `harmonic_angle.cuh` / `multi_harmonic_dihedral.cuh`
- 修改全仓库:include 路径
- `.claude/NOTES.md` — 加 multi-GPU CUDA gap 条目
- `.claude/architecture.md` — `/mol:map` 刷新
- `CMakeLists.txt` / `src/CMakeLists.txt`

## Tasks

- [ ] `git mv src/cuda/md/neighborlist/` → `src/cuda/neighborlist/`; update all `#include` paths; update CMakeLists.txt
- [ ] `git mv include/atomiverse/cuda/observer/stress.cuh` → `cuda/md/virial_stress.cuh`; `git mv src/cuda/observer/stress.cu` → `src/cuda/md/virial_stress.cu`; update all `#include` paths
- [ ] Fix `angle_harmonic.cuh:18` docstring: `[degrees]` → `[rad]` with API-contract note (`include/atomiverse/cuda/potentials/angle_harmonic.cuh:18`)
- [ ] Remove backward-compat alias headers: `include/atomiverse/cuda/potentials/{lj.cuh, harmonic_bond.cuh, harmonic_angle.cuh, multi_harmonic_dihedral.cuh}`
- [ ] Update all repo-wide `#include` of the alias headers to use canonical headers (`pair_lj_cut.cuh`, `bond_harmonic.cuh`, `angle_harmonic.cuh`, `dihedral_multi_harmonic.cuh`)
- [ ] Add multi-GPU halo-exchange CUDA gap entry to `.claude/NOTES.md` describing the divergence from CPU `#ifdef ATV_BUILD_MPI` behavior; mark "out of scope" for current chain
- [ ] Run `/mol:map` to refresh `.claude/architecture.md` (specifically remove the stale "`cuda::Frame` embeds `SCFWorkspace`" claim at line 91)
- [ ] (Documentation only) Note `cuda::Thermo` namespace pattern in NOTES if not already documented: `ATV::cuda::` for backward compat with CPU sibling, not `ATV::cuda::recorder::`
- [ ] Verify clean build after all renames + deletions
- [ ] Verify all tests pass

## Testing strategy

1. **Build** — `cmake --build build -j$(nproc)` 干净通过(`#include` 路径更新无遗漏)
2. **Regression** — 完整 `ctest --test-dir build --output-on-failure` 全过
3. **静态检查** — `grep -rn '#include.*cuda/observer/stress\|#include.*cuda/potentials/lj\.cuh\|#include.*cuda/potentials/harmonic_bond\.cuh\|#include.*cuda/potentials/harmonic_angle\.cuh\|#include.*cuda/potentials/multi_harmonic_dihedral\.cuh\|#include.*cuda/md/neighborlist/' include/ src/ tests/ benchmarks/` 返回 0 hits(全部更新到 canonical 路径)
4. **NOTES validation** — `.claude/NOTES.md` 中存在 "Multi-GPU halo exchange CUDA backend gap" 条目
5. **Architecture refresh** — `.claude/architecture.md` 不再含 `SCFWorkspace` 提及 in `cuda::Frame` 段

## Out of scope

- 多 GPU halo-exchange CUDA 端 **实现**(只跟踪,不实现)
- 其他 cluster 修复 — 已在 link 01-07 完成

## Notes

`/mol:map` 是交互式 skill,需要用户在本 link 实施时确认重写 `.claude/architecture.md`。`/mol:impl` 不应代为运行 `/mol:map` —— 该步骤标 manual checkpoint。
