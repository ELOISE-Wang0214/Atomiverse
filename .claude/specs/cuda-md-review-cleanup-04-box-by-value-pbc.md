---
title: CUDA MD Review Cleanup — link 04 — pass cuda::Box by value (PBC kernel reform)
slug: cuda-md-review-cleanup-04-box-by-value-pbc
status: done
created: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 6, 5 findings)
chain: cuda-md-review-cleanup-* (link 4/9)
---

# CUDA MD Review Cleanup — link 04 — pass cuda::Box by value (PBC kernel reform)

## Summary

`cuda::Box`(`include/atomiverse/cuda/box.cuh:31-91`)已缓存了 `inv_l[3]`、`half_l[3]`、`ortho` 等用于 PBC 的所有量,但 integrator 和 bonded kernel **不用** —— 每个 kernel 各自 inline 重写 PBC wrap,每线程做 FP64 division `1/Lx`(~20-30 cycle × 3 axes × N atoms)+ 读 `box[0],box[4],box[8]` 三个 uncoalesced load;bonded kernel 用 `round(dx/bx)` 而非缓存 `inv_l`;triclinic `minimum_image` 在 ortho=true 时仍走完整三斜路径(对 ortho 多算 4-5× FLOPS)。本 link 把所有 integrator(`VelocityVerlet`、`NVELimit`、`NoseHooverNPT`)和 bonded(`HarmonicBond`、`HarmonicAngle`、`MultiHarmonicDihedral`)kernel 改为按值传 `cuda::Box`(`trivially_copyable`),全部走 `cuda::Box::min_image`(自动选 ortho 快路径);把 `smart_rolling_prune`、`export_atom_pairs_kernel`、`fill_pairs`、`min_aabb_dist_sq` 按 `bool ortho` 模板特化。

## Domain basis

`CLAUDE.md ## Architecture § Backend directory layering` 把 `Box`(`cpu::Box` / `cuda::Box`)指定为 backend root 的"canonical simulation-cell type",owns `inv_cell`、per-axis `inv_L`、`half_L`、`bool ortho`,并暴露 `min_image(dx, dy, dz)`、`wrap_axis(q, d)`。CPU 同胞类 `cpu::Box` 已被全 backend 使用;CUDA 端 `cuda::Box` 的 8-14 行注释说"deferred migration"——本 link 落地。

PBC 算法 trade-off:
- ortho path:`dx -= L * nearbyint(dx * inv_L)`(3 个 axis,共 9 FMA + 3 rint)
- triclinic path:`s = inv_cell * dx; s -= round(s); dx = cell * s`(12 mul + 6 add + 3 rint + 9 mul + 6 add)≈ 4-5× 更多 FLOPS

`bool ortho` 模板特化让编译器在编译期消除 triclinic 分支,小系统/大系统均有 ~10-15% inner-loop FLOPS 减少。

`trivially_copyable<cuda::Box>` 保证按值传 `__device__` 参数无开销 —— 编译器把它放入 register。

## Design

完整 5 项发现见 inventory `cuda-md-review-cleanup.md` Cluster 6。本 link 增量设计:

1. **Integrator kernels**:
   - `velocity_verlet.cu:42-47`、`nve_limit.cu:62-67`、`nose_hoover_npt.cu:101-110`
   - kernel 签名加 `cuda::Box box`(by value);删除所有 `box[0]`/`box[4]`/`box[8]` 直接读和 `inv_Lx = 1/Lx`
   - PBC wrap 用 `box.wrap_axis(x, 0)`(自动 ortho fast path)
2. **Bonded kernels**:
   - `bond_harmonic.cu:53-55`、`angle_harmonic.cu:62-67`、`dihedral_multi_harmonic.cu:72-80`
   - 同样按值传 `cuda::Box`;`round(dx/bx)` → `box.min_image(dx, dy, dz)`
3. **NeighborList ortho 模板特化**:
   - `smart_rolling_prune`、`export_atom_pairs_kernel`、`fill_pairs`、`min_aabb_dist_sq` 加 `template<bool Ortho>` 参数
   - launcher 在 build 时根据 `box.ortho` 选择特化
4. **cuda::Box 增量**(如果需要):
   - 确保 `cuda::Box::min_image` 和 `wrap_axis` 在 `__device__` 上下文中正确 inline
   - 添加 `__forceinline__` 注解(若编译器不自动 inline)

## Files to create or modify

- `src/cuda/md/velocity_verlet.cu` — 改 kernel 签名,接 `cuda::Box`
- `src/cuda/md/nve_limit.cu` — 同上
- `src/cuda/md/nose_hoover_npt.cu` — 同上
- `src/cuda/potentials/bond_harmonic.cu` — kernel 签名 + `min_image`
- `src/cuda/potentials/angle_harmonic.cu` — 同上
- `src/cuda/potentials/dihedral_multi_harmonic.cu` — 同上
- `src/cuda/md/neighborlist/fused_search.cu` — `smart_rolling_prune` 模板特化
- `src/cuda/neighbor_list.cu` — `export_atom_pairs_kernel` 模板特化 + launcher dispatch
- `src/cuda/md/neighborlist/pair_search.cu` — `fill_pairs` 模板特化
- `src/cuda/md/neighborlist/triclinic.cuh` — `min_aabb_dist_sq` 模板特化(若适用)
- `include/atomiverse/cuda/box.cuh` — 确保 `min_image` / `wrap_axis` 有 `__forceinline__`

## Tasks

- [ ] Change `VelocityVerlet` kernel signature to take `cuda::Box` by value; delete inline `box[0]/box[4]/box[8]` reads and `1/Lx` divisions (`src/cuda/md/velocity_verlet.cu:42-47`)
- [ ] Same for `NVELimit` (`src/cuda/md/nve_limit.cu:62-67`)
- [ ] Same for `NoseHooverNPT` per-atom kernels (`src/cuda/md/nose_hoover_npt.cu:101-110`)
- [ ] Replace `round(dx/bx)` with `cuda::Box::min_image` in `HarmonicBond` kernel (`src/cuda/potentials/bond_harmonic.cu:53-55`)
- [ ] Same for `HarmonicAngle` (`src/cuda/potentials/angle_harmonic.cu:62-67`)
- [ ] Same for `MultiHarmonicDihedral` (`src/cuda/potentials/dihedral_multi_harmonic.cu:72-80`)
- [ ] Template-specialize `smart_rolling_prune` on `bool Ortho` (`src/cuda/md/neighborlist/fused_search.cu:186-203, 233-235`)
- [ ] Template-specialize `export_atom_pairs_kernel` on `bool Ortho` (`src/cuda/neighbor_list.cu:312`)
- [ ] Template-specialize `fill_pairs` on `bool Ortho` (`src/cuda/md/neighborlist/pair_search.cu:124-214`)
- [ ] Add launcher dispatch: select ortho vs triclinic specialization at launch time based on `box.ortho`
- [ ] Verify `cuda::Box::min_image` and `wrap_axis` inline in `__device__` context (use `cuobjdump` or compile-output inspection); add `__forceinline__` if needed (`include/atomiverse/cuda/box.cuh:31-91`)
- [ ] Verify all existing tests pass (no force/energy drift)

## Testing strategy

1. **Build + regression** — clean, `ctest -R 'cuda_md|cuda_potentials|cuda_neighbor'` 全过
2. **Correctness** — `tests/cuda/md/test_md.cu` 力/能量值不变(triclinic 和 ortho box 都跑)
3. **Perf check** — `bench_md_step` PE 用例(link 09 加),记录 baseline-vs-after;预计 ortho box 上 ≥ 1.05× 加速(PBC kernel 端);triclinic box 上不退化

## Out of scope

- 力 atomic 争用(cluster 7)→ link 05
- 6 次顺序 stress 规约 → link 05
- `__launch_bounds__` → link 06
- 力流池 → link 06
