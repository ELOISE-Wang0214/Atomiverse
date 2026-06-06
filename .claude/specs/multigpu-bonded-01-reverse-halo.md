---
title: Multi-GPU bonded — reverse force halo (reverse_comm)
status: code-complete
created: 2026-06-06
chain: multigpu-bonded
link: 1
---

# Multi-GPU bonded — reverse force halo (reverse_comm)

## Summary
链条第 1 环(基础)。`MultiGpu::reverse_comm()` 是 `forward_comm` 的伴随:按逆 cascade(z→y→x)把每个 dst 的
ghost 行力 gather、经 NVLink P2P 送回 src、scatter-add 到 src 的 owned/前轴-ghost 行。这是跨卡键力的前置——键力核在
ghost 行 atomicAdd 后,reverse_comm 把它累加回真正的 owner。LJ 路径算 owned 全力(ghost 力恒 0),不调用本函数。
**已实现并在 GH200 验证**(commit bdc1645)。

## Domain basis
无新物理。reverse-once 不变量:每步在所有力核之后恰调用一次,镜像单次 forward halo,Σ F 守恒不被破坏。

## Design
`gather_ghost_forces`(读 dst ghost 行 [base,base+cnt))/ `scatter_add_forces`(atomicAdd 进 src `d_src_idx` 行)两个核;
`HaloPlan` 增 dst 侧 reverse buffer `d_f{x,y,z}_dst`,src 侧 `d_s{x,y,z}` 复用作接收缓冲;在 `build_halo` 分配、
`destroy_halo_plan` 释放。

## Files to create or modify
- `src/cuda/md/multi_gpu.cu`、`include/atomiverse/cuda/md/multi_gpu.cuh`(已改)。

## Tasks
- [x] **P1.0** reverse_comm + gather/scatter kernels + HaloPlan reverse buffers (committed bdc1645)

## Testing
- GH200 build green;`test_cuda_multi_gpu` 13/13(LJ 不变)。reverse_comm 的端到端正确性由 -03 的键力等价测试覆盖。

## Out of scope
- 键力核本身(→ -02/-03/-04)。
