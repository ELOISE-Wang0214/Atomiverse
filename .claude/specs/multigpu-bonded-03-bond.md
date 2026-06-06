---
title: Multi-GPU bonded — cross-GPU BondHarmonic + reverse_comm wiring
status: code-complete
created: 2026-06-06
chain: multigpu-bonded
link: 3
---

# Multi-GPU bonded — cross-GPU BondHarmonic + reverse_comm wiring

## Summary
链条第 3 环。用 -02 重展开的 per-GPU bond tuple(partner 已解析到 owned/ghost local)在每张卡跑单卡现成的
`harmonic_bond_k`(对 i_local 与 partner_local 各 atomicAdd,partner 可能是 ghost),随后调 `reverse_comm()`(-01)
把 ghost 上的键力累加回 owner;能量只在 owner 卡计一次。验证多卡 PE(仅键)能量/力与单卡、CPU 一致。

## Domain basis
键力是局域两体力;每条 bond 恰由 min-tag owner 卡持有(不重不漏)。force-on-owned-and-ghost + reverse-once ⇒
owned 力完整、无双计。能量分解 E_bond 跨卡求和与单卡一致(浮点/混合精度容差内)。

## Design
- `MultiGpu` 增 bond force engine(per-GPU `cuda::BondHarmonic` 或直接调核),在 force 阶段:
  ZeroForces over owned+ghost → (LJ owned 全力,ghost=0) → bond 核 atomicAdd owned+ghost → `reverse_comm()`。
- 能量/virial:owner 卡累加,跨卡 host-side reduce(沿用现有 gather 风格)或留待 -05。
- bench/test 入口:多卡 PE-bond-only 与单卡对照。

## Files to create or modify
- `src/cuda/md/multi_gpu.cu`、`.cuh`(bond force + reverse 接线)。
- `tests/cuda/md/test_cuda_multi_gpu_bonded.cu`(扩展:键力等价)。

## Tasks
- [x] **P3.1** Wire per-GPU BondHarmonic over re-materialized bond tuples (deposit owned+ghost)
- [x] **P3.2** Call reverse_comm after bonded force; owner-only energy accumulation + cross-GPU reduce
- [x] **P3.3** RED→GREEN: 2/4-GPU bond energy + per-atom force == single-GPU == CPU; ΣF≈0

## Testing
- `test_cuda_multi_gpu_bonded`: 2/4-GPU vs single-GPU bond energy + forces within tol; ΣF≈0 each step; LJ tests green.

## Out of scope
- Angle/dihedral (→ -04); special-bonds (→ -04); optimization (→ -05).
