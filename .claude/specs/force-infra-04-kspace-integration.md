---
title: Force infra — KSpace as a first-class family: real/recip split + multi-GPU reduction path
status: approved
created: 2026-06-07
chain: force-infra
link: 4
---

# Force infra — KSpace as a first-class family: real/recip split + multi-GPU reduction path

## Summary
链条第 4 环。把**现有**长程静电（`KSpaceEwald` 直接 k 求和、`KSpacePME` FFT，两后端均已存在：
`ewald.{h,cuh}`、`pme.{h,cuh}`）整合进通用力基础设施，确保抽象**照顾 kspace**——而非只为短程 pair 设计。
**不新增任何 potential**；本环是整合 + 打通多卡，不是新物理。

三件事：
1. **形式化实/倒空间拆分**：短程实空间库仑（erfc）是 -03 的一个 pair-kernel；倒空间是独立 KSpace 族
   （-01 已给基类槽位）。`g_ewald`/α 由引擎单一下发给两侧（-01），对标 LAMMPS `kspace->init()` 先于
   `pair->init()`（`src/force.cpp:187-188`）、combined `lj/cut/coul/long` 从 `force->kspace->g_ewald` 取值
   （`pair_lj_cut_coul_long.cpp:674`）。
2. **承认 KSpace 自有的通信机制**（与 -02 的原子 forward/reverse **不同**）：Ewald 需跨域**结构因子 S(k) 归约**
   （CPU 已用 `MPI_Allreduce`，`ewald.cpp:130-133`）；PME 需**网格通信**（对标 LAMMPS `Grid3d` 的 brick↔FFT，
   `pppm.cpp:703-737`，与原子 comm 正交）。把这一"kspace 归约/网格 comm"作为 KSpace 族的接口契约固化，
   而非塞进原子 halo。
3. **打通多卡 kspace**：当前 CUDA kspace **仅单卡**（S(k) 只对 `n_local` 求和，`ewald.cuh:6`），且 `MultiGpu`
   不传电荷、不持 kspace。靠 -01（kspace 槽位）+ -02（电荷成为 forward 字段，喂实空间库仑）+ 本环的
   **跨 GPU S(k) 归约（NCCL）**，使多卡 Ewald 跑通。PME 多卡网格分解作为后续，明确标注不在本环。

## Domain basis
**行为保持**：单卡/CPU 现有 Ewald+Coul 总能量、逐原子力、应力与改造前在容差内一致（实/倒拆分只是重组接线，
物理不变；α 单一来源消除"双份手传 α 不一致"的隐患）。多卡 Ewald 正确性判据：N-GPU 的总静电能、逐原子力
== 单卡 == CPU（容差内）——证明电荷-forward（-02）+ 跨 GPU S(k) 归约共同正确。
Ewald 理论：E = Σ_i F_i(实空间 erfc) + 倒空间 Σ_k A(k)|S(k)|² + 自能 + 净电荷项；S(k)=Σ_i q_i e^{ik·r_i}
必须是**全局**和，故多卡须对各卡 partial S(k) 做求和归约（NCCL allreduce），与 CPU `MPI_Allreduce` 同构。

## Design
- **实/倒接线**：实空间库仑走 -03 pair-kernel（已在 -03 移植）；KSpace 族（-01 基类）只管倒空间。引擎按
  kspace→pair 顺序下发 g_ewald/α（-01 契约），单卡/CPU 行为不变。
- **KSpace 通信契约**：在 KSpace 基类上定义"全局归约/网格 comm"钩子（区别于 -02 的原子 comm）：
  - Ewald：`reduce_structure_factor()` —— CPU=`MPI_Allreduce`，CUDA 多卡=NCCL allreduce（单卡=no-op）。
  - PME：网格 forward/reverse 钩子（对标 Grid3d），**多卡网格分解明确留后续**，本环仅定义接口 + 保持单卡/复制网格现状。
- **多卡 Ewald**：`MultiGpu` 经 -01 注册 KSpaceEwald；电荷经 -02 内置 forward 字段到位；力阶段
  pair(LJ+实空间库仑) → bonded → kspace(倒空间)，kspace 用跨 GPU S(k) 归约。补上 CUDA Ewald 现缺的
  净电荷修正（`ewald.cu:167-168` 注明遗漏）。
- 不引入 Kokkos；不新增 potential。

## Files to create or modify
- `include/atomiverse/{cpu,cuda}/kspace/kspace.*`（KSpace 基类 + 归约/网格通信钩子）。
- `include/atomiverse/cuda/kspace/ewald.cuh`、`src/cuda/kspace/ewald.cu`(跨 GPU S(k) NCCL 归约 + 净电荷修正)。
- `src/cuda/md/multi_gpu.cu`、`.cuh`(注册 kspace + 力阶段顺序 + 电荷已由 -02 到位)。
- 引擎 g_ewald/α 单一下发（与 -01 协同）。
- 测试:`test_ewald_coul`(单卡/CPU 回归)+ `test_cuda_multi_gpu_ewald`(多卡 == 单卡)。

## Tasks
- [ ] **P4.1** Formalize real/reciprocal split: real-space coulomb is the -03 pair-kernel; KSpace family owns reciprocal; engine is the single source of g_ewald/α (kspace-before-pair)
- [ ] **P4.1b** (carried from -02) CUDA: carry `q` to ghosts in `MultiGpu::build_halo` ghost data — the per-atom charge prerequisite for multi-GPU real-space Coulomb (the -02 atom comm contract is in place on `cuda::ComponentBase`; this wires its first CUDA consumer)
- [ ] **P4.1c** (carried from -02) Generic per-step extra-field participant ITERATION: CPU `halo_exchange` iterates registered `comm_forward`/`comm_reverse` participants (batched buffer); `MultiGpu` forward/reverse iterate held components' `CommField`s into the combined Σ(width) staging. RED→GREEN synthetic 1-wide forward + 1-wide reverse-accumulate round-trip (CPU n8 + CUDA multi-GPU) — the deferred -02 ac-004. Consumed by EAM (`fp` forward / `rho` reverse) in -04 proper.
- [ ] **P4.2** Define the KSpace comm contract (global structure-factor reduction for Ewald; grid forward/reverse hooks for PME) — distinct from the -02 atom comm
- [ ] **P4.3** CUDA Ewald cross-GPU S(k) reduction via NCCL (single-GPU = no-op); add the missing net-charge correction
- [ ] **P4.4** Register KSpaceEwald on `MultiGpu` (charge supplied by P4.1b); run force phase pair→bonded→kspace
- [ ] **P4.5** RED→GREEN: single-GPU/CPU Ewald+Coul unchanged; multi-GPU Ewald total electrostatic energy + per-atom force == single-GPU == CPU within tol
- [ ] **P4.6** PME: define grid-comm interface + keep current single-GPU/replicated-grid behavior; document multi-GPU distributed-grid PME as a follow-up (no silent gap)

## Testing
- `test_ewald_coul` (CPU + single-GPU): total electrostatic energy, per-atom force, stress unchanged vs baseline within tol.
- `test_cuda_multi_gpu_ewald`: 2/4-GPU Ewald total electrostatic energy + per-atom forces == single-GPU == CPU within tol (proves charge-forward + cross-GPU S(k) reduction correct).
- Full CPU + CUDA ctest green.

## Out of scope
- New potentials (none; Ewald/PME/Coul all pre-existing). Multi-GPU distributed-grid PME (interface defined, impl deferred — explicitly logged). Dispersion-PME (g_ewald_6) / MSM.
