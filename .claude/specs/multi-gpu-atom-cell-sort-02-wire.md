---
title: Multi-GPU atom cell-sort — wire into reneighbor
status: in-progress
created: 2026-06-08
chain: multi-gpu-atom-cell-sort
link: 2
---

# Multi-GPU atom cell-sort — wire into reneighbor (link 2/3)

## Summary
链条第 2 环。把 link-01 的置换原语接入 `MultiGpu::reneighbor`,**始终开启**(无门控),由完整 parity 测试套
作为唯一绑定门。目标:多卡力与单卡逐字节一致,同时让 build/force 访存因 cell 排序而合并。

## Domain basis
保持物理不变(置换搬运全部物理量)。绑定校验面 = parity 测试(力 vs 单卡逐字节一致),非方程。

## Design
- 插入点(行号确认):`reneighbor()` 中 `R_migrate`(~1324)→ **[插入 sort]** → `R_build_halo`(~1327)
  → `R_rebuild_bonded`(~1329)。
- **顺序是关键且 load-bearing**:
  - 在 `build_halo()` **之前**排序 owned —— HaloPlan 存的 src 原子索引是对已排序 owned 集构建的,保持一致。
  - `build_tag_maps`(~964)与 `rebuild_bonded_tuples`(~1013)每次 reneighbor 从头重建、按 global-tag→local-row
    解析,排序在它们**之前**发生即可自动吸收新顺序。
- owned-only 排序(ghost 紧随其后由 `build_halo` 重建);本环 sort 始终开启,不加启发式。
- migrate 已在排序前完成 owned 集 finalize。

## Files
- `src/cuda/md/multi_gpu.cu`:`reneighbor` 接入置换调用;确认 `build_halo`/`build_tag_maps`/
  `rebuild_bonded_tuples` 在其后运行。
- link-01 的置换原语(调用)。

## Tasks
- [x] Wire the permute call into `reneighbor()` between `migrate()` and `build_halo()` (R_cell_sort, multi_gpu.cu).
- [x] Verify ordering: sort precedes build_halo / build_tag_maps / rebuild_bonded_tuples.
- [x] GREEN: 18/18 `test_cuda_multi_gpu` (1 disabled) + 9/9 `test_cuda_md` + 2/2 `test_cuda_neighborlist_skin_gating` 通过(GH200 job 11185)。
- [x] Hygiene: simplify ran clean.

> NOTE: link-02 uses a per-reneighbor local `AtomCellSortBuffers` (correctness-first); link-03 moves to reused per-device buffers + the beneficial-or-skip gate to remove the alloc/sort overhead.

## Testing
GH200 gtest:`test_cuda_multi_gpu`(尤其 FullPotentialParityVsSingleGpu / BondForceParityVsSingleGpu /
NveTrajectoryParityNlist / NveMigrationParity)、`test_cuda_md`、`test_cuda_neighborlist_skin_gating`。
sbatch partition gpu, module GPU/buildenv-gcccuda/2026.03-cu13.0。

## Out of scope
beneficial-or-skip 门控 + 性能 benchmark(link 03)。
