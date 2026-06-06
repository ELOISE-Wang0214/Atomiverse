---
title: Multi-GPU bonded — optimization (P2P batching, launch bounds, atomics)
status: in-progress
created: 2026-06-06
chain: multigpu-bonded
link: 5
---

# Multi-GPU bonded — optimization

## Summary
链条第 5 环(收尾)。在 -02..-04 正确的跨卡键力上做性能优化:合并每次 P2P 传输(xyz/force 打包成单缓冲)、
键力核 `__launch_bounds__`、owner-grouped atomicAdd 降竞争、forward+reverse 共用同一次 halo、消除每步分配。
不改变物理(正确性测试保持绿),报告多卡 PE step-time before/after。

## Domain basis
无新物理;纯实现优化。可证伪:优化后多卡 PE 的能量/力与优化前逐位/容差内一致(正确性不回归),且 step-time 改善有实测。

## Design
- 合并 reverse/forward 的 3 个 cudaMemcpyPeer 为 1 个(打包 [cnt*3]);减少 P2P 启动次数。
- 键力核加 `__launch_bounds__`;tag-map 二分用共享内存/缓存;re-materialize 与 LJ nlist 重建共用 owned+ghost 扫描。
- 复用 HaloPlan 缓冲(已在 -01 预分配),确认无每步 malloc/free。
- before/after profile(nsys / 计时)在多卡 PE。

## Files to create or modify
- `src/cuda/md/multi_gpu.cu`、`bonded_materialize.cu`、`tag_map.cu`、键力核(launch bounds)。

## Tasks
- [~] **P5.1** Batch P2P transfers (single packed buffer for xyz / forces)
- [x] **P5.2** __launch_bounds__ + owner-grouped atomics on bonded kernels
- [ ] **P5.3** Share the owned+ghost scan between re-materialize and nlist rebuild; eliminate per-step alloc
- [ ] **P5.4** Before/after multi-GPU PE step-time profile; correctness unchanged

## Testing
- All -02..-04 correctness tests stay green; multi-GPU PE step-time before/after recorded.

## Out of scope
- New potentials; the benchmark campaign (→ benchmark-matrix).
