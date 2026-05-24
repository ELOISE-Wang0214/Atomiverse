---
title: CUDA MD Review Cleanup — link 07 — shared kernel header alignment + microfixes
slug: cuda-md-review-cleanup-07-kernel-header-and-microfixes
status: draft
created: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 9 + Cluster 12, 8 findings)
chain: cuda-md-review-cleanup-* (link 7/9)
---

# CUDA MD Review Cleanup — link 07 — shared kernel header alignment + microfixes

## Summary

低风险扫尾:
- 修复 2 项共享 kernel 头声明与定义签名不符(`cuda::kernels::lj_force` 14 vs 18 参,`harmonic_bond_force` 类似)—— 选择 (a) 头补 14 参 overload 转发到 18 参 或 (b) 头改 18 参 + 更新 torch 后端调用点。**这是 cluster 9 的 2 项 🔴**(silent torch 后端 phantom 符号)
- 6 项 mini-fixes(cluster 12):`box_muller` cos/sin → `sincos`、`berendsen_scale_k` 在 host 算 cbrt、dihedral `dv_dcos` 预算系数数组、device `invert_3x3` 缓存、`deform.cu` thread-0 box update 折叠成单 `<<<1,1>>>` follow-up、3 个 bonded launch 评估是否加 gather kernel(可选)

放在 link 07 是因为 hot kernel 已经在 link 04-06 稳定;此时 header 改动和微优化不会和 hot kernel diff 冲突。

## Domain basis

**共享 kernel 头契约**:`include/atomiverse/cuda/potentials/lj_kernel.cuh` 顶部注释声称"single source of truth for LJ force+energy computation ... can be called from any backend",但 14-参签名 vs 18-参 `.cu` 定义。torch 后端 `src/torch/potentials/lj.cu:96-105` 调用 14 参 phantom 符号 —— 实际链接时找不到,要么报错要么链到旧 build artifact。

要选哪个 fix 策略:
- (a) 头补 14 参 overload 转发:对 torch 后端无侵入,可立即生效
- (b) 头改 18 参 + 更新 torch:把 stress + special-bonds CSR 暴露给 torch 后端,长期一致性更好,但侵入 torch 后端

推荐 (b):"single source of truth" 意味着所有 backend 用同一签名;phantom 符号是 latent bug,torch 后端可能任何时候编译失败。

**`sincos`**:NVIDIA CUDA 提供 `sincos(double, double*, double*)` 和 `sincosf(float, float*, float*)`,单 PTX op 同时算 sin + cos,比顺序两次调用快 ~2×。`sincospi(x, ...)` 对 `x = θ/π` 时有更精的实现。

**`berendsen_scale_k` host cbrt**:`μ = pow(1 + dt/τ_P * (P - P_target)/P_factor, 1.0/3.0)`(NPT volume rescale)。计算在 host 上一次,作 kernel arg 传入,代替每个 block 的 thread 0 重算。

**`invert_3x3` device 缓存**:cuda triclinic 把 cell matrix 3×3 求逆;每步重算浪费(NVE/NVT 不变)。在 `NeighborList::State` 缓存 inv_cell,box 变化时重新求逆。

## Design

完整 8 项发现见 inventory `cuda-md-review-cleanup.md` Cluster 9 + Cluster 12。增量:

1. **lj_kernel.cuh / harmonic_bond_kernel.cuh 签名对齐**(选 b)
2. **box_muller `sincos`**(`fix_langevin.cu:43-49`):`double s, c; sincos(2*PI*u1, &s, &c); z1 = sqrt(-2*log(u2)) * c; z2 = sqrt(-2*log(u2)) * s;`
3. **berendsen `cbrt` host 端**(`barostat.cu:40-73`):`Berendsen::launch()` 内先 host 算 `μ = cbrt(...)`,作为 kernel arg
4. **dihedral 系数预算**(`dihedral_multi_harmonic.cu:115-119`):`State` 加 `device_vector<work_t> dv_dcos_coeffs_`,build 时填 `{C[0], C[1], 2*C[2], 3*C[3], 4*C[4]}`
5. **invert_3x3 缓存**(`triclinic.cuh:47-62`):`NeighborList::State` 加 `inv_cell_cache_[9]`,box 变化时由 host 重算上传到 `__constant__`(与 link 02 setup_cell 一致策略)
6. **deform.cu thread-0 折叠**(`deform.cu:113`):分离为独立 `<<<1,1>>>` `update_box_diagonal_k` kernel,在主 propagate_x_k 之后调用
7. **bonded gather kernel**(可选):如果 link 05 cluster-staged 已经合并 E/stress 写,则本项无需做

## Files to create or modify

- `include/atomiverse/cuda/potentials/lj_kernel.cuh` — 改 18 参签名
- `src/torch/potentials/lj.cu` — 更新调用点匹配 18 参
- `include/atomiverse/cuda/potentials/harmonic_bond_kernel.cuh` — 同 LJ
- `src/torch/potentials/harmonic_bond.cu` — 同 LJ torch
- `src/cuda/md/fix_langevin.cu` — `sincos`
- `src/cuda/md/barostat.cu` — host-side `cbrt`
- `src/cuda/potentials/dihedral_multi_harmonic.cu` — 预算系数
- `include/atomiverse/cuda/potentials/dihedral_multi_harmonic.cuh` — `State` 字段
- `src/cuda/md/neighborlist/triclinic.cuh` — inv_cell 缓存(可能并入 `NeighborList::State`)
- `src/cuda/neighbor_list.cu` — inv_cell 缓存 launcher
- `src/cuda/md/deform.cu` — 分离 box-update kernel

## Tasks

- [ ] Align `cuda::kernels::lj_force` signature: change header to 18-arg (with stress + special-bonds CSR); update torch backend call site (`include/atomiverse/cuda/potentials/lj_kernel.cuh:32-41`, `src/torch/potentials/lj.cu:96-105`)
- [ ] Align `cuda::kernels::harmonic_bond_force` signature similarly; verify torch call site (`include/atomiverse/cuda/potentials/harmonic_bond_kernel.cuh:29-37`, `src/torch/potentials/harmonic_bond.cu`)
- [ ] Replace cos/sin in box_muller with sincos (single PTX op); consider sincospi (`src/cuda/md/fix_langevin.cu:43-49`)
- [ ] Move cbrt(μ) in berendsen_scale_k to host; pass μ as kernel arg (`src/cuda/md/barostat.cu:40-73`)
- [ ] Pre-compute dv_dcos coefficients `{C[0], C[1], 2*C[2], 3*C[3], 4*C[4]}` in `MultiHarmonicDihedral::build()`; store device array (`src/cuda/potentials/dihedral_multi_harmonic.cu:115-119`)
- [ ] Cache `inv_cell` (3×3 matrix) on host; re-compute only on box change; upload to `__constant__` (`src/cuda/md/neighborlist/triclinic.cuh:47-62`)
- [ ] Separate `if (i==0) { box[0] = ...; }` in `deform.cu:113` into a follow-up `<<<1,1>>>` `update_box_diagonal_k` kernel
- [ ] (Optional) Add downstream gather kernel for bonded E/stress if cluster 5 cluster-staged doesn't fold them in
- [ ] Verify correctness across all tests

## Testing strategy

1. **Build + regression** — clean
2. **Torch backend build** — `cmake --build build -DATV_BUILD_TORCH=ON` (if enabled) compiles cleanly; phantom symbol issue resolved
3. **Correctness** — `test_md.cu`、`test_bonded.cu`、`test_pair_lj.cu`、`test_fix_langevin.cu` 全过
4. **Perf softgate** — `box_muller` 微优化预计 0.5-1% step time 节省(very small),只要不退化即可

## Out of scope

- 文件重定位 → link 08
- 文档/docstring 修复(`angle_harmonic.cuh` [degrees]→[rad])→ link 08
- Alias header 删除 → link 08
- `/mol:map` 刷新 → link 08 / 09
