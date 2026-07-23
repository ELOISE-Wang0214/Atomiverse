---
spec: cuda-pme-accuracy
status: done
created: 2026-07-23
chain: cuda-pme-accuracy
phase: single
grilled: false
---

# cuda-pme-accuracy — CUDA PME 倒空间能量偏差 3%

## Summary

`tests/cuda/longrange/test_pme_coul.cu` 的 `CudaPmeCoul.MadelungNaCl` 失败:算出的
NaCl 静电能比解析 Madelung **低约 3.0%**,而 CPU 同参数达 5e-4。

## Root cause (verified host-side, 2026-07-23)

**不是 PME 倒空间精度问题。** CUDA `KSpacePME` 的 recip+self 与 CPU 逐位一致:

```
e_real (CPU CoulLong)     =  -8.6907 eV
e_pme  recip+self (CPU=CUDA host port) = -276.8618 eV
e_tot                     = -285.5525 eV   (= Madelung)
GPU reported total        = -276.8618 eV   ← 恰好等于 recip+self alone
```

`alpha=3.5/cutoff` 时实空间约贡献 3%;`alpha=5.0`(Ewald 测试)时实空间仅
0.07%,所以 kspace-only 仍落在 Ewald 的 2e-3 容差内 —— 假绿。

**真正的 bug:** `cuda::CoulLong` 仍走 cluster-CSR (`cl_row_ptr`),而
`NeighborList` 在 `flat_pairs_requested=false`(全树默认,RDF 迁 molrs 后
无人置位)时 **不填充** cl_*,只填充 per-atom Verlet 列表。`PairLJCut` 已切到
Verlet;CoulLong 静默 early-return → 实空间能量 = 0。

**Fix:** `CoulLong::compute` 改用与 `PairLJCut` 相同的 Verlet +
`md::launch_pair_compute` + `potentials::CoulKernel` 路径(cuFFT PME 不动)。

## Tasks

- [x] Host 二分:CPU real / recip / self 拆分 + CUDA 算法 host 复现
- [x] 修 `src/cuda/potentials/coul_long.cu` → Verlet 路径
- [x] GPU job 572602: `test_cuda_pme_coul` Passed
- [x] GPU job 572602: `test_cuda_ewald_coul` Passed

## Testing strategy

GPU 节点 `sbatch`(配方见 memory `cuda-gpu-node-build-recipe`)。二分用
`ctest -R test_cuda_pme_coul --output-on-failure`;对照值取 CPU 的
`tests/cpu/longrange/test_pme_coul.cpp` 与 `regressions/md/npt_pme_nacl`
(REFERENCE.md 记录同网格/阶数下 CPU 达 5e-4)。

## Out of scope

- 改 `ATV_CUDA_REAL` 默认值(保持 `float`)。
- CPU PME/Ewald 任何改动。
- 放宽测试容差来"修"失败 —— 容差是物理门禁,CPU 同参数能达标,CUDA 也应该。
