---
title: Multi-GPU atom cell-sort — permutation primitive
status: in-progress
created: 2026-06-08
chain: multi-gpu-atom-cell-sort
link: 1
---

# Multi-GPU atom cell-sort — permutation primitive (link 1/3)

## Summary
链条第 1 环。提供一个独立、单元测试覆盖的"按邻居网格 cell 索引对一个设备上的 owned 原子做物理重排"原语,
**暂不接入** `MultiGpu::reneighbor`(dead-but-tested,可独立合并)。这是后续修复 `build_verlet_direct_k`
非合并访存(多卡 LJ even-count 性能悬崖)的基础:把同一 cell 的原子在数组中排到一起,使 warp 内线程共享 cell。

## Domain basis
纯访存布局优化,**逐字节保持力不变**(置换是单射,物理量随原子一起搬运)。无新物理/方程;正确性由置换
往返恒等 + cell 单调性单元测试保证。前置已完成:力核 `pair_compute_verlet` 已用 AoS float4 修好
(ng2 +68%, ng4 +60%, 29/29 CUDA 测试通过);本链针对剩余的 build 核瓶颈。

## Design
- 新增一个原语(建议置于 `src/cuda/md/neighborlist/` 或 `src/cuda/md/`,与既有 `sort_and_cluster.cu`
  的 cub 依赖一致):对给定设备 frame 的 owned 原子计算 cell key(用 `domain_decomposition` 的 cell 几何,
  与 build 网格一致),`cub::DeviceRadixSort::SortPairs(key=cell, value=0..nl-1)` 得到置换 `perm`。
- 用 `perm` 对**所有** per-atom 数组做 gather 重排(写入临时再回拷,或 ping-pong 缓冲):
  - LJ:`x,y,z,vx,vy,vz,fx,fy,fz,mass,charge,type,tag,image`
  - bonded SoA:`num_bond, bond_atom[6], bond_type[6], num_angle, angle_atom[12], angle_type[6],
    num_dihedral, dihedral_atom[24], dihedral_type[8], nspecial[3], special[16]`
- 提供逆置换以支持往返测试。owned-only(占 97%);ghost 由 `build_halo` 重建,不在本原语职责内。
- **风险点**:漏掉任一 bonded 成员只在跨卡 bonded tuple 下静默破坏拓扑 —— 因此本环以单元测试先行隔离。
- 既有 `sort_and_cluster.cu` 的 Morton 排序是单卡 nlist 内部用途、**不**物理重排 owned frame;本原语是新的。

## Files
- `src/cuda/md/neighborlist/` 或 `src/cuda/md/`:新增置换原语(.cu/.cuh)。
- `include/atomiverse/cuda/frame.cuh`:per-atom + bonded SoA 字段(只读引用,确认全字段)。
- `include/atomiverse/cuda/domain_decomposition.cuh`:cell key 几何。
- `tests/cuda/md/`:新增单元测试。

## Tasks
- [x] Write failing unit test: known frame (positions + 全 bonded SoA) → permute → 断言每个数组按同一 perm 重排;perm∘perm⁻¹ 恒等;sorted cell 非降。 (tests/cuda/md/test_atom_cell_sort.cu)
- [x] Implement cell-key compute + `cub::DeviceRadixSort` permutation. (atom_cell_sort.cu)
- [x] Implement all-array gather (12 LJ + 全 bonded SoA) + inverse.
- [x] GREEN: 单元测试通过(GH200 job 11182:2/2 PASSED — AtomCellSort.ConsistentPermutationAndSortedByCell + PermuteThenInverseIsIdentity)。
- [x] Hygiene: simplify ran clean (cleaned one stale stage-sizing comment).

## Testing
单元测试 in `tests/cuda/md/`,GH200(sbatch,partition gpu,module GPU/buildenv-gcccuda/2026.03-cu13.0,
`LD_LIBRARY_PATH=$LIBRARY_PATH:...`)。build:`cmake --build build-cuda`。

## Out of scope
接入 `reneighbor`(link 02);性能门控 + benchmark(link 03)。
