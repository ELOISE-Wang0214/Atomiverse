---
title: Driver device unification — link 03 — collapse category classes into *Driver<DT> templates
slug: driver-device-unification-03-category-templates
status: draft
created: 2026-06-01
parent: driver-device-unification (design)
chain: driver-device-unification-* (link 3/4)
---

# Driver device unification — link 03 — collapse category classes into *Driver<DT> templates

## Summary

链 3/4：消灭 **category 层**的重复。把 link 02 仍并存的 7 个手写 category 类合并为 4 个设备泛型模板 `ATV::EvaluateDriver<DT>` / `DynamicsDriver<DT>` / `OptimizeDriver<DT>` / `ResponseDriver<DT>`（均继承 `BackendDriver<其 Phase, DT>`），每个的 `run()` 与 `validate_components()` **各只写一次**。`molrs::FrameRef` 重载上移到模板，CUDA 因此免费获得（之前仅 CPU 有）。具体 driver（`cpu::/cuda::md::MD`、`Evaluate`、`Minimize`、`cpu::vib::Hessian`）reparent 到模板。最后移除 link 01 引入时保留的旧 Backend 公开名（`set_frame`/`materialize`/`retrieve`/`drain_analysis`/`sync`），调用点全部走 `stage_*`。

## Domain basis

契约源：`CLAUDE.md ## Architecture § Driver taxonomy`（三/四个 category 由终止语义定义，CUDA 镜像同一套）。

前提（link 02 已落地）：`BackendDriver<P,DT>` 承载共享胶水 + setter + 注入；7 个 category 类已挂在其下。

各 category 的设备无关 `run()`/`validate` 语义（已逐字核对当前实现，两后端仅差 `stage_*` 名）：
- Evaluate：setup / step_once(0) / teardown；validate = ≥1 `Evaluator` 且 `capabilities()` 非空。
- Dynamics：续跑（`live_`/`frame_ptr_`/`step_counter_`），段末 `stage_out`+`drain`，full teardown 延迟到 dtor；validate = ≥1 `Integrator` + ≥1 `Capability::Force`。
- Optimize：`max_iter` 循环，`converged||stop_requested` 早停；validate = ≥1 `Optimizer` + ≥1 `Evaluator(Energy|Force)`。
- Response：6N + finalize，无 outer_status 轮询；validate = ≥1 `Force` + ≥1 `Hessian`。

## Design

- `facade/driver_categories.h` + `_impl.h`：四个模板，`run()`/`validate_components()` 各一份；Dynamics 模板加 `set_integrator/set_thermostat/set_barostat`；各模板带 `molrs::FrameRef` 重载（materialize cpu::Frame → 调 `cpu::Frame&` 重载 → write_back），对两 DT 都实例化。
- reparent：`cpu::md::MD : DynamicsDriver<cpu::device_tag>`、`cuda::md::MD : DynamicsDriver<cuda::device_tag>`；`Evaluate`→`EvaluateDriver<DT>`；`Minimize`→`OptimizeDriver<DT>`；`cpu::vib::Hessian`→`ResponseDriver<cpu::device_tag>`。
- 删 7 个 category 头（cpu evaluate/dynamics/optimize/response + cuda evaluate/dynamics/optimize）。
- 移除旧 Backend 公开名，`stage_*` 成唯一公开暂存接口（private 实现可保留旧名）。
- 显式实例化：`src/cpu/driver.cpp` 加 `Evaluate/Dynamics/Optimize/Response Driver<cpu::device_tag>`；`src/cuda/driver.cu` 加 `Evaluate/Dynamics/Optimize Driver<cuda::device_tag>`（Response 暂不实例化 CUDA）。

## Files to create or modify

**Create:**
- `include/atomiverse/facade/driver_categories.h` — 4 个模板声明 + Dynamics 专属 setter + molrs 重载
- `include/atomiverse/facade/driver_categories_impl.h` — `run()`/`validate_components()` body

**Modify:**
- `include/atomiverse/cpu/md/md.h` / `cuda/md/md.cuh` — reparent 到 `DynamicsDriver<DT>`
- `include/atomiverse/cpu/evaluate.h` / `cuda/evaluate.cuh` — reparent
- `include/atomiverse/cpu/minimizer/minimize.h` / `cuda/minimizer/minimize.cuh` — reparent
- `include/atomiverse/cpu/vib/hessian.h` — reparent `ResponseDriver<cpu::device_tag>`
- `include/atomiverse/cpu/backend.h` + `src/cpu/backend.cpp` — 移除旧公开名
- `include/atomiverse/cuda/backend.cuh` + `src/cuda/driver.cu` — 移除旧公开名
- `src/cpu/driver.cpp` / `src/cuda/driver.cu` — 显式实例化 category 模板

**Delete:**
- `include/atomiverse/cpu/{evaluate,dynamics,optimize,response}_driver.h`
- `include/atomiverse/cuda/{evaluate,dynamics,optimize}_driver.cuh`

## Tasks

- [ ] Write failing test `tests/test_category_templates.cpp`: a `DynamicsDriver<cpu::device_tag>` and (smoke) `<cuda::device_tag>` run a tiny LJ NVE trajectory; continuation (`run(5); run(5)`) keeps a cumulative step counter
- [ ] Create `facade/driver_categories.h` + `_impl.h` with the 4 device-generic templates (one run()/validate each), Dynamics setters, and molrs::FrameRef overloads
- [ ] Reparent concrete drivers (md::MD, Evaluate, Minimize on both backends; cpu::vib::Hessian) onto the templates
- [ ] Delete the 7 category headers; fix every include site
- [ ] Remove old Backend public names (set_frame/materialize/retrieve/drain_analysis/sync); route all callers through stage_in/stage_out/drain
- [ ] Update explicit instantiations in src/cpu/driver.cpp + src/cuda/driver.cu
- [ ] Build CPU+CUDA clean; full ctest green (continuation/NHC/Langevin state preserved)
- [ ] Hygiene cleanup via `/mol:simplify`

## Testing strategy

1. **Build** CPU+CUDA clean (all category instantiations).
2. **New** `test_category_templates` — Dynamics NVE + continuation counter parity.
3. **Regression** — full ctest green; specifically `test_md`, `cuda_md*`, minimizer, evaluate, vib/hessian, hot-path no-sync guards (unified teardown adds no per-step sync).
4. **Grep gate** — old Backend public names gone; 7 category headers deleted, no residual include.

## Out of scope

- Plan-equivalence + validation-parity cross-backend tests + CLAUDE.md taxonomy diagram → link 04
- Merging concrete drivers into single `md::MD<DT>` templates → future (needs all wired components traits-typed)
- CUDA Response/Hessian instantiation → future
