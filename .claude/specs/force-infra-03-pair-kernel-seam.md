---
title: Force infra — generic multi-term pair-kernel seam + shared per-backend traversal
status: in-progress
created: 2026-06-07
chain: force-infra
link: 3
---

# Force infra — generic multi-term pair-kernel seam + shared per-backend traversal

## Summary
链条第 3 环。当前每个短程势都重写整套"遍历邻居表 + min-image + 截断 + 归约 + 力累加"：
`PairLJCut` 与 `CoulLong` 的内层循环**逐字重复**，仅"逐对力数学"不同（LJ 的 r⁻¹² vs 库仑 erfc，
见 `src/cpu/potentials/pair_lj_cut.cpp:229+` 与 `coul_long.cpp:153-220`，CUDA 同理 cluster 走两遍）。
加任何短程势 = 复制粘贴。

本环引入一个**通用、与具体势无关的 pair-kernel 缝隙**——**不是 `LJExpr`**。接口是一个泛型概念
`PairKernel`：给定 `(rsq, qi, qj, itype, jtype)`，**产出一个或多个力/能量贡献**。LJ、库仑实空间各是一个
实例；二者**可在同一个邻居循环里叠加**——正是 LAMMPS combined `lj/cut/coul/long`（一次遍历算 LJ + erfc 两项，
`src/KSPACE/pair_lj_cut_coul_long.cpp:133-159`）。共享遍历/归约 launcher 各后端各写一份（CPU xsimd ket-loop、
CUDA 列主序 Verlet），对标 LAMMPS-KOKKOS `PairComputeFunctor` + 每势仅 `compute_fpair`（`src/KOKKOS/pair_kokkos.h:129-197`）。
**抄模式不抄 Kokkos**：AV 手写内核已实测快过 LAMMPS-KOKKOS，故保留手写 launcher，只抽出力数学。

验证用**现有**两个不同函数（LJ + 库仑实空间），证明缝隙不绑 LJ；**不新增任何 potential**。

## Domain basis
**行为/性能保持**：`PairLJCut`、`CoulLong` 移植到缝隙后，逐原子力、能量、virial 与改造前逐位一致；
单卡/多卡/CPU 每核吞吐相对改造前无回归（容差内，对照已记录 AV 基线）。多贡献接口的正确性判据：
组合 kernel(LJ+coul) 的合力 == 分别跑 LJ 与 CoulLong 之和（同一邻居表，叠加律）。kernel 用模板标量 `T`
实例化（CPU `double`/`xsimd::batch`、CUDA device `real_t`），力数学单一真值源、三处复用。

## Design
- 通用 pair-kernel 概念（后端无关、CUDA-free 头，**命名与势无关**），形如：
  ```cpp
  struct PairContribution { work_t f_over_r; work_t energy; };
  // 任意短程势实现该 concept；可返回多项之和或通过 out 累加多贡献
  template<class T> AV_HD PairContribution operator()(T rsq, T qi, T qj, int it, int jt) const;
  ```
  LJ kernel、Coulomb-real kernel 各是一个实现；组合 kernel = 两者叠加（一次遍历）。
- CUDA launcher `pair_compute<Kernel>(ctx, kernel, pairlist, box, periodic)`：复用现有逐原子 Verlet kernel，
  内层"算 rsq → 调 kernel → 累加 → block_reduce"模板化。
- CPU launcher `pair_compute<Kernel>(...)`：复用现有 OMP + xsimd ket-loop（gather/批量 rsq/掩码/scatter/
  per-thread 归约 + Voigt6），逐对数学换成 kernel 的 batch 实例化；标量尾同样调 kernel。
- 移植**现有** `PairLJCut` 与 `CoulLong`（两后端）到 launcher；CPU `CoulLong` 当前是 `ComponentBase`、
  CUDA 是 `Pair`——本环顺带消解此不对称（统一为同一 pair-kernel 缝隙的两个实例）。
- 与 -01 协同：缝隙化的势仍是 `cuda::Pair` 子类，经 -01 工厂注册。不引入 Kokkos / 不新增 potential。

## Files to create or modify
- `include/atomiverse/potentials/pair_kernel.h`（通用 concept + 贡献类型，新建，命名与势无关）。
- `include/atomiverse/potentials/{lj,coul_long}_kernel.h`（把**现有** LJ / 库仑实空间数学搬成 kernel 实例）。
- `include/atomiverse/cuda/potentials/pair_compute.cuh`、`include/atomiverse/cpu/md/pair_compute.h`（共享遍历 launcher，新建）。
- `src/cuda/potentials/pair_lj_cut.cu`、`coul_long.cu`；`src/cpu/potentials/pair_lj_cut.cpp`、`coul_long.cpp`（移植到 launcher，瘦身）。
- 测试：`test_pair_kernel`（kernel vs 解析）+ 组合叠加等价 + LJ/Coul 移植 parity。

## Tasks
- [x] **P3.1** Define the backend-agnostic `PairKernel` concept + multi-contribution type (templated scalar `T`); no LJ in the interface name/shape
- [x] **P3.2** Express existing LJ and Coulomb-real-space math as kernel instances
- [~] **P3.3** CUDA `pair_compute<Kernel>` launcher (`pair_compute.cuh`) over the Verlet list; **`PairLJCut` ported** (hand `lj_force_verlet` deleted), GH200 27/27 + racecheck clean. CoulLong CUDA port still pending (it's a `Pair` cluster style; port with the CPU CoulLong in P3.4).
- [ ] **P3.4** CPU `pair_compute<Kernel>` launcher over the xsimd ket-loop; port `PairLJCut` + `CoulLong`; unify the CPU/CUDA CoulLong base asymmetry
- [x] **P3.5** kernel vs analytic (LJ + coul) + combined-kernel == LJ+Coul superposition + batch==scalar: `test_cpu_pair_kernel` green. Migrated CUDA `PairLJCut` bit-identical (GH200 27/27 incl. multi_gpu/PE). CPU-migration parity pends P3.4.
- [x] **P3.6** Perf: single-GPU 1M PE = **331 step/s vs ~340 hand-kernel baseline (−2.7%); still 1.15× over LAMMPS (288)**. User-accepted (2026-06-07): the generality/maintainability win is the chain's purpose; small regression documented, perf follow-up (functor param-cast hoist via `kernel.prepare<real_t>()`) deferred. CPU per-core LJ pends P3.4.

## Testing
- `test_pair_kernel`: LJ and coul kernels return correct f_over_r/energy vs analytic.
- Superposition: a combined LJ+coul kernel in one neighbor walk == sum of separate LJ and CoulLong runs (forces/energy within tol).
- Parity: migrated PairLJCut + CoulLong bit-identical to baseline (CPU n1/n8, CUDA 1/2/4-GPU).
- Perf gate (operator GH200 + CPU node): no throughput regression vs baseline.

## Out of scope
- New physics potentials (none added; only existing LJ + Coul migrated). Manybody (EAM/Tersoff). Reciprocal kspace (→ -04). Cross-backend single code path beyond the kernel (traversal stays per-backend by design).
