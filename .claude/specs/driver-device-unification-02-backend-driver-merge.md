---
title: Driver device unification — link 02 — merge CPUDriver + CUDADriver into BackendDriver<P,DT>
slug: driver-device-unification-02-backend-driver-merge
status: draft
created: 2026-06-01
parent: driver-device-unification (design)
chain: driver-device-unification-* (link 2/4)
---

# Driver device unification — link 02 — merge CPUDriver + CUDADriver into BackendDriver<P,DT>

## Summary

链 2/4：删除手写两遍的 backend 胶水层 `cpu::CPUDriver<P>` + `cuda::CUDADriver<P>`，合并为单一 facade 模板 `ATV::BackendDriver<PhaseEnum, DeviceTag>`。它经 link 01 的 `stage_in/stage_out/drain` 统一接口驱动 setup/teardown，承载共享的 `build()`/循环辅助/`install_defaults()`/`scan_pending_infra_()` 与 8 个共享 categorical setter，并统一采用 typed-query 注入（删 CPU 的 `name()=="..."` 字符串探测与 `dynamic_cast<PairEvaluator>`）。`backend_traits<DT>` 扩展 glue 类型供 setter/注入用。

本 link **保留 7 个 category 类**（`cpu::/cuda::` × Evaluate/Dynamics/Optimize + CPU Response），仅把它们的父类从 `cpu::CPUDriver`/`cuda::CUDADriver` 改挂到 `BackendDriver<其 Phase, DT>`；category 的 `run()`/`validate_components()` 暂仍各自存在（link 03 才合并为模板）。即：本 link 消灭 **backend 胶水层**的重复，link 03 再消灭 **category 层**的重复。

## Domain basis

契约源：`CLAUDE.md ## Architecture § Canonical facade`（`Driver<PhaseEnum, DeviceTag>` 底层原语，用户经 category 基类使用）+ `§ Categorical setter API`（`set_pair` 等按 capability 自动注入）。

前提（link 01 已落地）：两个 `Backend` 都有 `stage_in/stage_out/drain`；`cpu::ComponentBase` 有 `provides_neighbor_list()`。

## Design

### `backend_traits<DT>` 扩展 glue 类型

```cpp
template<> struct backend_traits<cpu::device_tag> {
  using component_type      = cpu::ComponentBase;    // 既有
  using backend_type        = cpu::Backend;          // 既有
  using pair_evaluator_base = cpu::PairEvaluator;    // set_pair static_assert + cutoff query
  using thermo_type         = cpu::Thermo;           // set_thermo 注入
  using virial_stress_type  = cpu::VirialStress;     // install_defaults 注入
};
// cuda::device_tag 同构（cuda::PairEvaluator / cuda::Thermo / cuda::VirialStress）
```

### `ATV::BackendDriver<P,DT>`（facade，CUDA-free 头 + `_impl.h`）

- `build()`：`install_defaults()` → `validate_components()` → `Base::build()`（保持现 CPU 顺序）。
- 循环辅助（**经 `stage_*`，设备无关**，非模板成员 body 进 `_impl.h` 按 `(P,DT)` 显式实例化）：`setup_pass(frame)`（init/load/`stage_in`/build_components/reset outer_status）、`step_once(i)`、`teardown_pass(frame)`（destroy_components/`stage_out`/shutdown）、`outer_status()`/`mutable_outer_status()`/`for_each_pending`/`throw_missing`。
- `validate_components()`：纯虚（category 覆写，本 link 仍在各 category 类里）。
- `install_defaults()`：虚，默认做**设备无关** typed 注入——`reads()&Stress` 命中且无 stress-provider 时注入 `traits::virial_stress_type{}`；`needs_neighbor_list()` 命中且无 `provides_neighbor_list()` 时抛错。
- **8 个共享 setter**（模板成员，header-inline，不碰 backend）：`set_pair/set_bond/set_angle/set_dihedral/set_neighbor/set_observer/set_record/set_thermo`，经 `CategoryTraits<P>`（eval_phase/end_phase）+ traits glue 类型参数化。

### typed-query 统一

`scan_pending_infra_()` 只用 `reads()`/`capabilities()`/`needs_neighbor_list()`/`provides_neighbor_list()`；pair cutoff 经 `traits::pair_evaluator_base` typed 访问器取（`dynamic_cast<traits::pair_evaluator_base*>` 至多一处）。删全部 `name()=="virial_stress"`/`name()=="neighbor_list"`。

### 删除 + reparent

- 删 `include/atomiverse/cpu/cpu_driver.h`、`include/atomiverse/cuda/cuda_driver.cuh`。
- `cpu::{Evaluate,Dynamics,Optimize,Response}Driver` 与 `cuda::{Evaluate,Dynamics,Optimize}Driver` 父类改为 `ATV::BackendDriver<其 Phase, DT>`；它们原先从 `CPUDriver`/`CUDADriver` 继承的成员（setup_pass/setter/...）现来自 `BackendDriver`。
- `cuda::md::MD` 原先自带的 8 个共享 setter + `install_infra_` 删除（改继承 `BackendDriver` 版本）；保留其 Dynamics 专属 `set_integrator/set_thermostat/set_barostat`。CPU-only MPI halo 注入：`cpu::DynamicsDriver`(或 `cpu::md::MD`) 覆写 `install_defaults()` 先调 `Base::install_defaults()` 再追加。

## Files to create or modify

**Create:**
- `include/atomiverse/facade/backend_driver.h` — `BackendDriver<P,DT>` 声明 + header-inline 共享 setters + `CategoryTraits`（从 cpu_driver.h 迁入）
- `include/atomiverse/facade/backend_driver_impl.h` — 非模板 body（setup/teardown/build/install_defaults/scan）

**Modify:**
- `include/atomiverse/cpu/backend_traits.h` / `include/atomiverse/cuda/backend_traits.cuh` — 加 glue 类型
- `include/atomiverse/cpu/{evaluate,dynamics,optimize,response}_driver.h` — 父类改 `BackendDriver`
- `include/atomiverse/cuda/{evaluate,dynamics,optimize}_driver.cuh` — 父类改 `BackendDriver`
- `include/atomiverse/cuda/md/md.cuh` — 删共享 setter + `install_infra_`（保留 Dynamics 专属 setter）
- `include/atomiverse/cpu/md/md.h` — MPI-halo 注入移入 `install_defaults()` override
- `src/cpu/driver.cpp` / `src/cuda/driver.cu` — 显式实例化加 `BackendDriver<P,DT>`
- `CMakeLists.txt`（如需）— 新 facade 头无需编译单元，但确认 include 路径

**Delete:**
- `include/atomiverse/cpu/cpu_driver.h`
- `include/atomiverse/cuda/cuda_driver.cuh`

## Tasks

- [ ] Write failing test `tests/test_backend_driver.cpp`: a `BackendDriver`-derived stub builds a Plan via the 8 shared setters and auto-injects VirialStress via typed `reads()==Stress` (no name strings)
- [ ] Extend `backend_traits<cpu::device_tag>` + `<cuda::device_tag>` with `pair_evaluator_base`, `thermo_type`, `virial_stress_type`
- [ ] Create `facade/backend_driver.h` + `backend_driver_impl.h` (build / setup_pass / step_once / teardown_pass via stage_*; outer_status; for_each_pending; throw_missing; install_defaults typed-inject; scan_pending_infra_ typed-only; 8 shared setters + CategoryTraits)
- [ ] Reparent the 7 category classes onto `BackendDriver<Phase, DT>`; delete `cpu_driver.h` + `cuda_driver.cuh`; fix includes
- [ ] Delete the duplicated shared setters + `install_infra_` from `cuda::md::MD`; move CPU MPI-halo into `install_defaults()` override
- [ ] `grep` gate: no `name() == "virial_stress"` / `name() == "neighbor_list"` and ≤1 `dynamic_cast<*PairEvaluator*>` across driver/injection code
- [ ] Update explicit instantiations; build CPU+CUDA clean; full ctest green

## Testing strategy

1. **Build** CPU+CUDA clean.
2. **New** `test_backend_driver` — shared-setter Plan build + typed VirialStress injection.
3. **Regression** — full ctest green (no behavior change: same Plans, same throws).
4. **Grep gate** — string proxies gone; ≤1 pair-evaluator dynamic_cast.

## Out of scope

- Collapsing the 7 category classes into `*Driver<DT>` templates → link 03
- Removing old Backend public names → link 03
- Parity tests + docs → link 04
