---
title: Force infra — engine composes force families by base pointer (Pair/Bond/…/KSpace), de-hardcode
status: in-progress
created: 2026-06-07
chain: force-infra
link: 1
---

# Force infra — engine composes force families by base pointer (Pair/Bond/…/KSpace), de-hardcode

## Summary
链条第 1 环。把力的"组合"做成通用基础设施，对标 LAMMPS `Force` 并列持有
`Pair* / KSpace* / Bond* / Angle* / Dihedral*`（`src/force.h:79,95`）。当前 `cuda::md::MultiGpu` 按
**具体类型**硬编码（`std::vector<std::unique_ptr<cuda::PairLJCut>> pair_;`，
`include/atomiverse/cuda/md/multi_gpu.cuh`；构造点 `src/cuda/md/multi_gpu.cu:870-876` 写死 `PairLJCut`），
且**完全不持有 KSpace**——多卡路径连 `CoulLong`/`Ewald` 都进不来。

本环：力引擎（单卡驱动 + `MultiGpu`）改为按**抽象基类指针 + 工厂/注册表**组合所有力族——
`cuda::Pair`、`cuda::Bond/Angle/Dihedral`，**以及 KSpace 族**（`cuda::KSpaceEwald`/`KSpacePME`，
现为 `ComponentBase` 子类，`ewald.cuh:22`/`pme.cuh:20`）。**没有任何一个力族被特殊对待**；LJ 只是
`Pair` 的一个实例，Ewald 只是 KSpace 的一个实例。默认组合 = 现有 LJ + harmonic（行为零变化）。
CPU 驱动已经是模板化 `set_pair<P>` 持基类 + `on(force).add(KSpace…)`，本环对齐其语义并补上 MultiGpu 的缺口。

## Domain basis
纯结构改造，**行为保持**：默认仍组合 LJ + harmonic，数值逐位不变（多态派发只把 `static_cast` 换成虚表）。
引入 KSpace 槽位不改变现有单卡 Ewald/Coul 结果（它们已能在单卡驱动里 `on(force).add()`）。
对标 LAMMPS 的力族初始化顺序：`kspace->init()` 必须先于 `pair->init()` 以共享 `g_ewald`
（`src/force.cpp:187-188`）——本环在引擎组合层固化这一顺序契约，取代现状"α 分别手传给 Coul 和 Ewald"
的易错做法（`coul_long` 与 `ewald` 各自构造参数 α，无单一真值源）。

## Design
- 力族基类指针 + 工厂（不含任何具体势名）：
  - `MultiGpu` 成员改为 `std::unique_ptr<cuda::Pair>` / `cuda::Bond` / `…` / `cuda::KSpace`
    （若无 `cuda::KSpace` 公共基类，则本环顺带提取一个：`KSpaceEwald`/`KSpacePME` 的公共
    `virtual compute/launch` 接口，与 LAMMPS `KSpace` 基类对应，`src/kspace.h`）。
  - `setup_force_engine` 接收调用方注入的工厂（`PairFactory`/`KSpaceFactory`/…）；缺省工厂 = LJ + harmonic、
    无 kspace（与今日行为一致）。
- 组合顺序契约：引擎按 `kspace.init/build → pair.init/build` 顺序构建，`g_ewald`/α 由引擎单一来源下发给
  Pair 实空间项与 KSpace（取代双份手传 α）。力阶段执行顺序对齐 LAMMPS：pair → bonded → kspace
  （`src/verlet.cpp:316-329`）。
- `compute_forces` / `destroy_force_engine` 一律走基类虚函数。
- 不引入 Kokkos；不新增任何 potential（LJ/Coul/Ewald/PME 均已存在）。

## Files to create or modify
- `include/atomiverse/cuda/md/multi_gpu.cuh`、`src/cuda/md/multi_gpu.cu`（成员→基类指针、工厂、kspace 槽位、组合顺序）。
- `include/atomiverse/cuda/kspace/kspace.cuh`（若需：提取 `cuda::KSpace` 公共基类；`KSpaceEwald/PME` 改继承之）。
- 单卡驱动侧（`cpu_driver`/`cuda_driver` 力族组合处）对齐 α/g_ewald 单一来源 + init 顺序。
- 测试：自定义工厂注入非默认 Pair/KSpace 桩，验证多态组合；现有 LJ/PE/Ewald 测试保持绿。

## Tasks
- [x] **P1.1** Introduce/confirm a `KSpace` base (cpu+cuda) with virtual build/launch; `KSpaceEwald`/`KSpacePME` inherit it
- [x] **P1.2** `MultiGpu` holds `Pair/Bond/Angle/Dihedral/KSpace` base pointers; `compute_forces`/`destroy_force_engine` dispatch virtually
- [x] **P1.3** Add factory seam to `setup_force_engine` (Pair/KSpace/bonded factories); default = LJ + harmonic, no kspace (behavior-preserving)
- [ ] **P1.4** Single source for g_ewald/α: engine sets it on KSpace then Pair (LAMMPS init order); remove the duplicated α constructor args as the authoritative source
- [ ] **P1.5** RED→GREEN: custom factory injects a stub `Pair` + stub `KSpace`, asserts both dispatched; existing LJ/PE bit-identical; existing single-GPU/CPU Ewald+Coul unchanged

## Testing
- Custom-factory test: stub `Pair` and stub `KSpace` are launched on each subdomain (polymorphic compose).
- Regression: `test_cuda_multi_gpu` (LJ/PE) bit-identical; `test_ewald_coul` (CPU + single-GPU) unchanged.
- Full CPU + CUDA ctest green.

## Out of scope
- Generic per-atom comm contract incl. charge (→ -02). Generic pair-kernel seam (→ -03). Multi-GPU kspace reduction (→ -04). No new potentials.
