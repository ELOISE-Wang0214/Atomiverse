---
title: CUDA MD Review Cleanup — link 05 — warp-shuffle reductions + cluster-staged force atomics
slug: cuda-md-review-cleanup-05-reductions-and-atomics
status: draft
created: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 5 + Cluster 7, 7 findings)
chain: cuda-md-review-cleanup-* (link 5/9)
---

# CUDA MD Review Cleanup — link 05 — warp-shuffle reductions + cluster-staged force atomics

## Summary

两组改造合并到一个 link:

1. **规约模式现代化(cluster 5,5 项)** — 用 sm_90 硬件加速的 `__reduce_add_sync` / `cooperative_groups::reduce` 替换 hot kernel 中 6 次顺序 shared-mem stress 规约(每次完整 `__syncthreads()` barrier);把 `kinetic_stress_k` 与 `KineticEnergy` 融合;`reduce_ke_axes_k` 改 warp-shuffle;删除 `if (s_energy[0] != 0) atomicAdd` 代码味。
2. **力 atomic 争用消除(cluster 7,2 项)** — LJ pair 每 pair 6 个 atomicAdd(稠密系统 ~10-30% step time)改为 cluster-staged 写回(shared-mem 累加 + 每 block 每 atom 一次 atomicAdd,~64× 争用降低);bonded kernel 按 dest atom 分组(LAMMPS-style per-atom owner list),寄存器累加,最后写一次。

两组都触碰相同 hot force kernel(`pair_lj_cut.cu`、`bond_harmonic.cu`、`angle_harmonic.cu`、`dihedral_multi_harmonic.cu`),合并避免相同文件被两轮 PR 改动。

## Domain basis

**Warp-shuffle 规约**:Hopper sm_90 提供 `__reduce_add_sync(mask, value)`,一个指令完成 warp 内 32-way sum,延迟 ~10 cycles。等价于 5 次 `__shfl_xor_sync` + 5 次 add 的 13-15 cycles。比 shared-mem 规约(块大小 256 → log2(256) = 8 iter × 5-10 cycles barrier = ~60-80 cycles)快 ~5×。

`cooperative_groups::reduce` 是 portable warp/block reduction,内部用 `__reduce_add_sync` 在支持的硬件上。

**Cluster-staged 力写回**:邻居列表已经把原子分组成 cluster(8 atoms/cluster)。让一个 thread block 处理一个 cluster pair 时,partial force 累加到 `__shared__ work_t fi_local[8][3]`(64 atoms × 3 = 192 floats = 768 bytes),最后只为本 cluster 的 8 atoms 各做一次 atomicAdd 到 global force。相比每 pair 6 个 atomicAdd,一个 16×16 cluster pair (256 pairs) 原来 1536 atomicAdds → 16 atomicAdds,~96× 减少。

**Bonded 按 dest 分组**:LAMMPS neighbor 模式把每个 atom 的"作为某个 bond/angle/dihedral 中心 atom"的所有 topology 实例列在该 atom 的列表里。每个 dest atom 的 contribution 在该 atom 的 thread 内累加到 register,最后写一次。代价:topology 数据结构需要按 dest 反向索引重组(目前是按 topology entry 顺序)。

## Design

完整 7 项发现见 inventory `cuda-md-review-cleanup.md` Cluster 5 + Cluster 7。本 link 增量:

1. **统一 warp-shuffle 规约 helper**:`include/atomiverse/cuda/cg_reduce.cuh`(新)— 提供 `warp_reduce_sum<T, N>(T values[N])`,内部用 `cooperative_groups::reduce` 或 `__reduce_add_sync`;打包 7 个值(1 E + 6 stress)一次性规约
2. **替换 force kernel 规约**:`pair_lj_cut.cu:140-155`、`bond_harmonic.cu:97-110`、`angle_harmonic.cu:145-158`、`dihedral_multi_harmonic.cu:237-250` 各自替换
3. **融合 kinetic_stress_k 与 KineticEnergy**:观察:都是 per-atom velocity 规约。融合为一个 kernel 写 `{ke_scalar, stress_xx, stress_yy, stress_zz, stress_xy, stress_xz, stress_yz}` 7 个值
4. **删 `if (s_energy[0] != 0)`**:atomicAdd(0) 在现代 GPU 上免费,条件判断反而引入分支
5. **Cluster-staged 力写回**(LJ pair):
   - 一个 block 处理一组 cluster pair
   - `__shared__ work_t fi_local[CLUSTER_SIZE][3]` + `fj_local[CLUSTER_SIZE][3]` 初始化 0
   - 每 thread 算自己的 pair contribution 累加到 shared
   - 块结束时,每 thread 对应一个 atom,做一次 atomicAdd 到 global force
6. **Bonded 按 dest 分组**:
   - 新增 `BondedOwnerList` 数据结构(在 `cuda::Frame` 之外的工具结构,build 时计算)
   - 每 dest atom 一个 thread,迭代该 atom 的 owner 列表,寄存器累加,最后写一次
   - 注:对 dihedral 来说"dest atom"概念需推广(4-atom topology),可选方案 — 仍用 per-topology thread,只是把 fi_local 加在 shared

## Files to create or modify

- (新)`include/atomiverse/cuda/cg_reduce.cuh` — warp-shuffle 规约 helper
- `src/cuda/potentials/pair_lj_cut.cu` — 规约替换 + cluster-staged 力写回
- `src/cuda/potentials/bond_harmonic.cu` — 规约替换 + (可选)按 dest 分组
- `src/cuda/potentials/angle_harmonic.cu` — 同上
- `src/cuda/potentials/dihedral_multi_harmonic.cu` — 规约替换;dihedral atomic 优化策略待评估
- `src/cuda/observer/stress.cu` — 融合 `kinetic_stress_k` 与 `KineticEnergy`(可能拆为 `src/cuda/md/kinetic_energy.cu` 接管)
- `src/cuda/md/kinetic_energy.cu` — 接受融合
- `src/cuda/md/nose_hoover_npt.cu` — `reduce_ke_axes_k` 改 warp-shuffle
- (可能新)`include/atomiverse/cuda/md/bonded_owner_list.cuh` — 按 dest 分组的反向索引

## Tasks

- [ ] Add `include/atomiverse/cuda/cg_reduce.cuh` warp-shuffle reduction helper using `cooperative_groups::reduce` or `__reduce_add_sync(0xffffffff, val)` for sm_90
- [ ] Replace 6 sequential shared-mem stress reductions in `pair_lj_cut.cu` with packed 7-value warp-shuffle (`src/cuda/potentials/pair_lj_cut.cu:140-155`)
- [ ] Same for `bond_harmonic.cu` (`src/cuda/potentials/bond_harmonic.cu:97-110`)
- [ ] Same for `angle_harmonic.cu` (`src/cuda/potentials/angle_harmonic.cu:145-158`)
- [ ] Same for `dihedral_multi_harmonic.cu` (`src/cuda/potentials/dihedral_multi_harmonic.cu:237-250`)
- [ ] Fuse `kinetic_stress_k` with `KineticEnergy` into a single kernel that writes 7 values `{ke_scalar, stress_xx..yz}` (`src/cuda/observer/stress.cu:103-105` + `src/cuda/md/kinetic_energy.cu`)
- [ ] `kinetic_stress_k` shared-mem reduction switched to warp-shuffle (~384 bytes shared instead of 12 KiB) (`src/cuda/observer/stress.cu:43-66`)
- [ ] `reduce_ke_axes_k` switched to warp-shuffle 3-axis reduction (`src/cuda/md/nose_hoover_npt.cu:129-167`)
- [ ] Remove `if (s_energy[0] != 0) atomicAdd(...)` pattern across hot kernels; rely on `atomicAdd(0)` being free
- [ ] Implement cluster-staged force write in `pair_lj_cut.cu`: shared-mem partials + one atomicAdd per atom per block (`src/cuda/potentials/pair_lj_cut.cu:102-107`)
- [ ] Implement per-atom owner-list reverse index for bonded topology; bond/angle kernels accumulate in register and write once per dest atom (`src/cuda/potentials/bond_harmonic.cu:69-74`, `angle_harmonic.cu:111-119`)
- [ ] Evaluate dihedral atomic strategy: either same owner-list approach (4-atom topology complicates) OR shared-mem partial accumulation per block (`src/cuda/potentials/dihedral_multi_harmonic.cu:197-208`)
- [ ] Verify correctness: force/energy/stress match pre-fix values within existing tolerances

## Testing strategy

1. **Build + regression** — clean, `ctest -R 'cuda_md|cuda_potentials|cuda_observer'` 全过
2. **Correctness** — 强检验:`tests/cuda/test_pair_lj.cu`、`test_bonded.cu`、`test_cuda_nose_hoover_npt.cu` 力/能量/stress 比对在原误差带内
3. **Perf check** — `bench_md_step` PE 用例,记录 baseline-vs-after;预计 LJ 力 kernel ≥ 1.5×,稠密 PE 系统总 step time ≥ 1.10×

## Out of scope

- `__launch_bounds__` 标注 → link 06
- 力流池 / 多 stream → link 06
- `inv_r = 1/r` FP64 div → link 06
