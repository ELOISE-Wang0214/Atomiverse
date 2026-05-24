---
title: CUDA MD Review Cleanup — link 02 — NeighborList Verlet-skin gating + over-launch fixes
slug: cuda-md-review-cleanup-02-nlist-skin-gating
status: done
created: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 3 + Cluster 4, 8 findings)
chain: cuda-md-review-cleanup-* (link 2/9)
---

# CUDA MD Review Cleanup — link 02 — NeighborList Verlet-skin gating + over-launch fixes

## Summary

接通已经分配但 **从未读写** 的 Verlet-skin 重建机制(`launch_cluster_displacement` kernel、`ref_cluster_{x,y,z}` buffer、`rebuild_flag`),让 `cuda::NeighborList::launch()` 在原子位移未超过 half-skin 时跳过整条重建管线(Morton sort、CUB radix sort、AABB、grid、coarse count、fill、compact、csr_to_flat、prune、export)。同时修复 4 项过启动 / 浪费:`launch_setup_cell_device<<<1,1>>>` 每步重跑、`smart_rolling_prune` / `export_atom_pairs_kernel` 用 `pair_capacity` 而非真实 pair 数定 grid、`setup_cell_from_*<<<1,1>>>` 单线程 O(stencil³)、`cut_skin_sq` / `cutoff_sq` 每步重算、CUB temp `cudaMalloc` 风险。**这是 chain 中单项最大 perf 收益**(预计邻居列表成本 ~28× 加速;~30-50% step time 释放)。关联 user MEMORY `[[cuda-nbnxn-smart-prune-missing]]`。

## Domain basis

Verlet-skin 重建判据:cutoff_extended = cutoff + skin;只要任意原子相对上次构建的位移 < skin/2,邻居列表保证不遗漏 pair(标准 Verlet list 理论)。skin = 0.3 Å @ T = 300 K LJ 系统典型重建周期 ~30 步。当前实现每步全量重建,等价于 skin = 0。

无穷小位移上界:`max_displacement²_per_atom = (skin/2)²`;`launch_cluster_displacement` 在 cluster 粒度(8 atoms/cluster)上求 max,写 `rebuild_flag = 1` 当任一 cluster max > skin/2。1 步滞后(在 step N 末尾算 step N+1 的 flag)对正确性无影响 —— skin 设计本就是缓冲。

参考代码:`src/cuda/md/neighborlist/fused_search.cu:116-136` 已有 `launch_cluster_displacement` 完整实现,但 `src/cuda/neighbor_list.cu:177-320` 的 `NeighborList::launch()` 没调用它;`nlist_impl.cuh:54,68` 分配了 `ref_cluster_x/y/z` + `rebuild_flag` 但没人写。

## Design

完整 8 项发现见 inventory `cuda-md-review-cleanup.md` Cluster 3 + Cluster 4。本 link 增量设计:

1. **重建 gate**(cluster 3):
   - 重建后立即把 cluster centers 写到 `ref_cluster_{x,y,z}`(用 `update_ref_clusters` kernel)
   - 下一步入口先 `launch_cluster_displacement`,写 `rebuild_flag`(device-side single int)
   - host 端用 lagged `cudaMemcpyAsync` D→H + `cudaEventQuery` poll 决定下一步是否跳过(1-step 滞后接受)
   - 若 `rebuild_flag == 0`,`launch()` 早退,只更新现有 `interaction_mask` / `flat_pair_list` 缓存
2. **Over-launch 修复**(cluster 4):
   - `launch_setup_cell_device<<<1,1>>>`:host 缓存 dims / n_stencil,只在 box volume 变化超阈值时重启
   - `smart_rolling_prune` / `export_atom_pairs_kernel` grid:用滞后读的 `flat_pair_counter` 而非 `pair_capacity`
   - `setup_cell_from_*<<<1,1>>>`:host 在 `build()` 算一次 O(stencil³) 数据,`cudaMemcpyToSymbolAsync` 上传到 `__constant__`;box 变化时重传
3. **状态缓存**(cluster 4):
   - `cut_skin_sq = r_search * r_search` 和 `cutoff_sq` 缓存在 `NeighborList::State` 内,`build()` 时一次性算
   - CUB temp storage 在 `build()` 预分配到最大尺寸(查询 CUB API)

## Files to create or modify

- `src/cuda/neighbor_list.cu` — `NeighborList::launch()` 重组,加 skin gate;`State` 加 cache 字段
- `src/cuda/md/neighborlist/nlist_impl.cuh` — `State` 字段扩展(`ref_cluster_*` 已存在;加 `rebuild_flag_host`、`cut_skin_sq`、`cutoff_sq`、`dims_cached`、`n_stencil_cached`)
- `src/cuda/md/neighborlist/fused_search.cu` — `launch_cluster_displacement` 已存在;补一个 `update_ref_clusters` kernel + launcher
- `src/cuda/md/neighborlist/setup.cu` — `setup_cell_from_*` 改为可选(只在 box 变化时调);新增 `__constant__` 拷贝路径
- `src/cuda/md/neighborlist/sort_and_cluster.cu` — CUB temp 预分配在 `build()` 时执行
- `include/atomiverse/cuda/neighbor_list.cuh` — 公共契约不变,但 `launch()` 文档更新说明 skin gating 行为
- (新)`tests/cuda/md/test_neighborlist_skin_gating.cu` — 验证 30 步不动原子 → radix sort 调用 ≤ 2 次

## Tasks

- [ ] Wire `launch_cluster_displacement` into `NeighborList::launch()` entry; consume `rebuild_flag` to gate the rebuild pipeline (`src/cuda/neighbor_list.cu:177-320`, `fused_search.cu:116-136`)
- [ ] Add `update_ref_clusters` launcher: snapshot cluster centers into `ref_cluster_{x,y,z}` after every rebuild (`src/cuda/md/neighborlist/fused_search.cu`)
- [ ] Add lagged D2H + `cudaEventQuery` polling for `rebuild_flag` (1-step lag acceptable; do not block hot path)
- [ ] Cache `dims_x/y/z`, `n_stencil` in host-side `State`; only re-run `launch_setup_cell_device<<<1,1>>>` when box volume changes beyond threshold (`src/cuda/neighbor_list.cu:190-194`)
- [ ] Replace `pair_capacity`-sized grid with `flat_pair_counter`-sized for `smart_rolling_prune` (`src/cuda/md/neighborlist/fused_search.cu:233-235`) and `export_atom_pairs_kernel` (`src/cuda/neighbor_list.cu:312`)
- [ ] Move `setup_cell_from_*<<<1,1>>>` work to host computation in `build()` + `cudaMemcpyToSymbolAsync` upload to `__constant__`; re-upload only when box changes (`src/cuda/md/neighborlist/setup.cu:178,193`)
- [ ] Cache `cut_skin_sq` and `cutoff_sq` in `NeighborList::State`, computed once in `build()` (`src/cuda/neighbor_list.cu:233, 293`)
- [ ] Pre-allocate CUB temp storage at max in `build()` (query CUB API for max size) (`src/cuda/md/neighborlist/sort_and_cluster.cu:112-120`)
- [ ] Move `cudaMemsetAsync(interaction_mask, 0xFF, pair_capacity*8)` from per-step to `build()`, OR memset only the csr_to_flat-touched portion (`src/cuda/neighbor_list.cu:267-269`)
- [ ] Add test `tests/cuda/md/test_neighborlist_skin_gating.cu`: 30-step MD loop with skin=0.3·cutoff and atoms held stationary (zero force), assert radix-sort kernel invocation count ≤ 2 (initial build + 1 safety margin)
- [ ] Verify existing `tests/cuda/md/test_md.cu` and `tests/cuda/test_nbnxn.cu` (if it exists) still pass

## Testing strategy

1. **Build + regression** — `cmake --build build` + `ctest --test-dir build -R cuda_md` 全过
2. **Skin gating runtime test** — `test_neighborlist_skin_gating.cu`:30 步零力 MD,断言 CUB radix sort kernel 调用 ≤ 2
3. **PE-realistic perf check**(softgate)— 在 GH200 上跑 `bench_md_step` 的 PE 用例(由 link 09 添加),记录 baseline-vs-after 邻居列表 wall-clock;预计 ≥ 1.30× 加速(若 < 1.10×,说明 skin gate 没起作用)
4. **Correctness** — `test_md.cu` 全跑通(力和能量不变)

## Out of scope

- `launch()` 内的 D2H+sync 修复 → link 03(本 link 改后,`launch()` 内 `interaction_mask` memset 移到 `build()` 也减少了一次 device 写,但不直接修 sync)
- `cuda::Box` 按值传 → link 04
- atomic 争用 → link 07(实际上 cluster 7 在 link 05)
- 文件重定位 `src/cuda/md/neighborlist/` → link 08
