---
title: Force infra — generic per-atom comm contract (charge + arbitrary fields), kill the Allgatherv hack
status: in-progress
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
`ComponentBase` 声明每原子通信宽度 + 实现 pack/unpack 回调，halo 服务遍历**注册的额外字段参与者**。
**按 LAMMPS 真实架构**（见下 Design decision）：位置-正向 / 力-反向**保持内置、逐步不变**（proven path），
**电荷随 build-time border 消息**到 ghost（每次 reneighbor 一次，点电荷恒定）——正是 LAMMPS `borders()`/AtomVec
携带 `q` 的等价物（Pair 实空间库仑读 `q[j]`，j 可为 ghost，`pair_lj_cut_coul_long.cpp:118-141`）。参与者契约
专供**额外**每步字段（EAM `fp`/`rho`）。`CoulLong` 改读 ghost 电荷、删掉 Allgatherv。

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

## Design decision (2026-06-07, user: cleanest/most-maintainable/extensible/performant)
Re-examining LAMMPS settles the architecture: `comm->forward_comm()` carries **positions**
intrinsically, `reverse_comm()` carries **forces** intrinsically, **charge rides the
`borders()`/AtomVec path** (sent once per reneighbor — point charges are static), and
`comm->forward_comm(Pair*)` packs only **extra per-step fields** (e.g. EAM `fp`). So the
LAMMPS-faithful = cleanest = lowest-risk shape is **intrinsic core + participant extras**,
NOT "everything is a participant". pos/force stay on the proven per-step path; charge joins
the build-time border message; the participant pack/unpack contract is for EXTRA fields.

## Design
- `ComponentBase` (`cpu/component.h` + `cuda/component.cuh`): add the **extra-field**
  participant contract (default 0 / no-op, old components unchanged): `comm_forward()` /
  `comm_reverse()` widths + pack/unpack. CPU = LAMMPS-style host `double*` callbacks; CUDA =
  device field descriptor (per-atom device ptr + width) consumed by MultiGpu's gather/scatter.
- **Charge → build-time border set** (not a per-step participant): CPU `build_halo` ghost
  message becomes `x,y,z,tag,q`; CUDA `build_halo` ghost data carries `q`. Sent once per
  reneighbor (mirrors LAMMPS `borders()`/AtomVec carrying `q`). Per-step forward (positions) /
  reverse (forces) are **intrinsic and unchanged** (proven path, no churn).
- Halo iterates registered **extra-field** participants, batched into the per-swap buffer
  (Σ extra widths) after the intrinsic positions (forward) / before forces (reverse).
- `CoulLong` (CPU) reads ghost `q` from the border set; delete the `tag_charge_` Allgatherv
  path (`coul_long.h:55-71` `use_charge_map_`/`refresh_charge_map`).
- 不引入 Kokkos；不新增 potential。kspace 倒空间的 S(k) 归约 / PME 网格通信是**另一种机制**，留 -04。

## Files to create or modify
- `include/atomiverse/cpu/component.h`、`include/atomiverse/cuda/component.cuh`（契约）。
- `include/atomiverse/cpu/comm/halo_exchange.h`、`src/cpu/comm/halo_exchange.cpp`（遍历参与者 + 内置 pos/force/charge）。
- `include/atomiverse/cuda/md/multi_gpu.cuh`、`src/cuda/md/multi_gpu.cu`（同上，P2P 按 Σ(width)）。
- `src/cpu/potentials/coul_long.{h,cpp}`（改用 ghost 电荷，删 Allgatherv）。
- 测试：`test_comm_contract`（合成字段 round-trip）+ CoulLong 多 rank/多卡等价回归。

## Tasks
- [x] **P2.1** Add the extra-field participant contract to CPU & CUDA `ComponentBase` (default 0/no-op): `comm_forward()`/`comm_reverse()` widths + pack/unpack (CPU host `double*`; CUDA device field descriptor)
- [x] **P2.2** CPU: carry charge `q` in the `build_halo` border message (x,y,z,tag,q); per-step pos-forward/force-reverse unchanged. (The generic per-step extra-field participant *iteration* is deferred to -04, where EAM `fp`/`rho` is its first consumer — no consumer in link 02; YAGNI.)
- [ ] **P2.3** CUDA: carry `q` in `build_halo` ghost data; per-step pos/force unchanged. (Extra-field iteration likewise lands in -04 via MultiGpu querying its held components' `CommField`s.)
- [x] **P2.4** Port CPU `CoulLong` onto ghost `q` from the border set; delete the MPI_Allgatherv tag→charge map
- [x] **P2.5** RED→GREEN (CPU): `test_mpi_coul_energy_bruteforce` — distributed real-space Coulomb energy == O(N²) erfc brute force at n1 AND n8 (rel < 1e-9), proving ghost charge via the border. (The synthetic per-step participant round-trip → -04 with the iteration.)
- [x] **P2.6** Position/force halo bit-identical — CPU: `test_mpi_halo_exchange`/`test_mpi_lj_energy_bruteforce` n1/n2/n8 green (unchanged by the 4→5 border size). CUDA portion pending P2.3.

## Testing
- `test_comm_contract` (CPU + CUDA): participant registers an extra per-atom field; forward broadcasts owner→ghost, reverse sums ghost→owner once.
- `test_ewald_coul` / coul multi-rank: real-space coulomb forces/energy match pre-change (charge via contract, not Allgatherv).
- Regression: `test_mpi_lj_energy_bruteforce` (n1/n8), `test_cuda_multi_gpu` (LJ/PE) bit-identical.

## Out of scope
- Reciprocal-space S(k) reduction / PME grid comm (a separate kspace mechanism → -04). Generic pair-kernel seam (→ -03). No new potentials. Migration unchanged.
