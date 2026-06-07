---
title: Force infra — generic per-atom comm contract (charge + arbitrary fields), kill the Allgatherv hack
status: approved
created: 2026-06-07
chain: force-infra
link: 2
---

# Force infra — generic per-atom comm contract (charge + arbitrary fields), kill the Allgatherv hack

## Summary
链条第 2 环，**基石**。当前 halo 通信字段**硬编码**：CPU `src/cpu/comm/halo_exchange.cpp` 正向只
`{x,y,z,tag}`、反向只 `{fx,fy,fz}`；CUDA `src/cuda/md/multi_gpu.cu:forward_comm(462-496)` 只 `x,y,z`、
`reverse_comm(498-541)` 只 `fx,fy,fz`。这逼出两处现成的丑陋：
- CPU `CoulLong` 拿不到 ghost 电荷，只能每步用 `MPI_Allgatherv` 建全局 `tag→charge` 表
  （`coul_long.cpp:50-83,74-75`）——O(N) 全局通信，纯属补丁。
- CUDA `MultiGpu` 根本没传电荷，故多卡跑不了 `CoulLong`/`Ewald`。

本环照搬 **LAMMPS 通用通信契约**（`src/comm.h`、`src/comm_brick.cpp:910-982`、`src/pair.h:50-52,207-210`）：
`ComponentBase` 声明每原子通信宽度 + 实现 pack/unpack 回调，halo 服务遍历**注册的通信参与者**而非硬编码字段。
现有位置-正向 / 力-反向退化为两个内置参与者（行为不变）；**电荷**成为第三个内置正向字段——
正是 LAMMPS 里电荷随原子边界自动到 ghost 的等价物（`src/verlet.cpp` borders；Pair 实空间库仑读 `q[j]`，
j 可为 ghost，`pair_lj_cut_coul_long.cpp:118-141`）。`CoulLong` 改读 ghost 电荷、删掉 Allgatherv。

LAMMPS 契约（直接对标）：
```cpp
int  comm_forward; int comm_reverse;                                  // pair.h:50-51
virtual int  pack_forward_comm (int n,int*list,double*buf,int,int*);  // pair.h:207
virtual void unpack_forward_comm(int n,int first,double*buf);         // pair.h:208
virtual int  pack_reverse_comm (int n,int first,double*buf);          // pair.h:209
virtual void unpack_reverse_comm(int n,int*list,double*buf);          // pair.h:210
```

## Domain basis
**行为保持 + 一处净化**：位置/力内置参与者与改造前逐位一致（cascade 顺序、PBC shift、reverse-once 动量守恒不变）。
电荷作为内置正向字段后，`CoulLong` 多卡/多 rank 的实空间库仑能量、逐原子力与原 Allgatherv 路径在容差内一致
（同一物理，换通信通道）。契约对任意每原子字段通用（宽度由组件声明），为后续 kspace 电荷同步与多体中间场
（如 EAM 密度）预留通道——但本环不实现任何新势。

## Design
- `ComponentBase`（`include/atomiverse/cpu/component.h` + `include/atomiverse/cuda/component.cuh`）增通信契约
  （默认 size=0 / no-op，旧组件零改动）：标量 `comm_forward/comm_reverse` + 四个 pack/unpack 回调（贴合 LAMMPS 签名；
  CUDA 侧为设备 gather/scatter 包装）。
- CPU `halo_exchange`、CUDA `MultiGpu::forward_comm/reverse_comm` 改为遍历注册参与者，缓冲区按 Σ(width)·cnt。
  内置参与者：① position-forward ② force-reverse ③ **charge-forward**。
- `CoulLong`（CPU）注册 charge-forward 需求（或直接消费内置电荷字段），删除 `tag_charge_` Allgatherv 路径
  （`coul_long.h:55-71` 的 `use_charge_map_`/`refresh_charge_map`）。
- 不引入 Kokkos；不新增 potential。kspace 倒空间的 S(k) 归约 / PME 网格通信是**另一种机制**，不在本环
  （留 -04）；本环只管原子级 forward/reverse + 电荷。

## Files to create or modify
- `include/atomiverse/cpu/component.h`、`include/atomiverse/cuda/component.cuh`（契约）。
- `include/atomiverse/cpu/comm/halo_exchange.h`、`src/cpu/comm/halo_exchange.cpp`（遍历参与者 + 内置 pos/force/charge）。
- `include/atomiverse/cuda/md/multi_gpu.cuh`、`src/cuda/md/multi_gpu.cu`（同上，P2P 按 Σ(width)）。
- `src/cpu/potentials/coul_long.{h,cpp}`（改用 ghost 电荷，删 Allgatherv）。
- 测试：`test_comm_contract`（合成字段 round-trip）+ CoulLong 多 rank/多卡等价回归。

## Tasks
- [ ] **P2.1** Add `comm_forward`/`comm_reverse` + pack/unpack forward/reverse to CPU & CUDA `ComponentBase` (default no-op)
- [ ] **P2.2** Refactor CPU `halo_exchange` to iterate participants; built-ins position-forward + force-reverse + charge-forward (behavior-preserving for pos/force)
- [ ] **P2.3** Refactor CUDA `MultiGpu` forward/reverse to iterate participants; buffer sized by Σ(width); add charge-forward built-in
- [ ] **P2.4** Port CPU `CoulLong` onto ghost-charge from the contract; delete the MPI_Allgatherv tag→charge map
- [ ] **P2.5** RED→GREEN: synthetic 1-wide forward + 1-wide reverse-accumulate field round-trips exactly (CPU n8 + CUDA multi-GPU); CoulLong multi-rank == single == pre-change
- [ ] **P2.6** Verify position/force halo bit-identical (CPU n1/n8, CUDA 1/2/4-GPU)

## Testing
- `test_comm_contract` (CPU + CUDA): participant registers an extra per-atom field; forward broadcasts owner→ghost, reverse sums ghost→owner once.
- `test_ewald_coul` / coul multi-rank: real-space coulomb forces/energy match pre-change (charge via contract, not Allgatherv).
- Regression: `test_mpi_lj_energy_bruteforce` (n1/n8), `test_cuda_multi_gpu` (LJ/PE) bit-identical.

## Out of scope
- Reciprocal-space S(k) reduction / PME grid comm (a separate kspace mechanism → -04). Generic pair-kernel seam (→ -03). No new potentials. Migration unchanged.
