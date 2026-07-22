---
spec: backend-mirror-05-dissolve
status: approved
created: 2026-07-22
chain: backend-mirror
phase: 05-dissolve
depends_on: backend-mirror-04-driver
grilled: true
---

# backend-mirror-05-dissolve — 拆解 Backend 上帝类

## Summary

删除公共 `cpu::Backend` / `cuda::Backend` 上帝类，把职责拆进已有完整词 infra 目录：`communication/`（域分解与 MPI）、`device/`（线程/stream/timer）、`asynchronous/`（分析执行器）、`driver/`（薄 Plan 分发与组件生命周期）。**禁止**新建 `runtime/` 或 `loop/`。调用方与测试全部改写，无 Backend 别名。

## Domain basis

无新物理。零拷贝 CPU Frame 绑定与 CUDA materialize/retrieve 边界行为保持。

## Design

| 职责 | 目标 |
|------|------|
| DomainDecomposition 拥有/bind | `communication/` |
| halo/migrate/reduce/gather 编排入口 | `communication/` |
| OpenMP / CUDA device+streams | `device/` |
| TimerRegistry | `device/` 或 core timer + 后端包装 |
| AnalysisExecutor | `asynchronous/` |
| Plan 绑定、component build/destroy、step/run_stage | `driver/` 内最小调度类型（可由 `Driver<P>` 直接持有，不必再叫 Backend） |
| Frame materialize/transfer | **仅** `frame/` 既有 API |

Driver 预设的 `build()`/`run()` 改为组装上述碎片，而不是持有一个全能 Backend。

### Reuse decision

| Candidate | Decision |
|-----------|----------|
| Backend methods/bodies | **generalize（dissolve）** into communication/device/asynchronous/driver |
| `frame/materialize`, `transfer` | **reuse** |
| `communication` domain decomp types | **reuse** / **generalize** ownership out of Backend |
| `device_manager` (cuda) | **reuse** — Backend stops owning duplicate lifecycle |

## Files to create or modify

- Delete or gut: `include/atomiverse/cpu/driver/backend.h`, `src` counterparts; cuda `backend.cuh`
- `include/atomiverse/cpu/communication/*`, `device/*`, `asynchronous/*`, `driver/*`
- `include/atomiverse/cuda/device/*`, `communication/*`, `driver/*`
- All Driver::build/run call sites, tests (`test_backend.cu`, cpu backend tests), examples

## Tasks

- [ ] Write failing tests: no public `class Backend` in cpu/cuda driver headers; DD ownership test targets communication types; stream/thread setup targets device
- [ ] Move CPU DD/MPI lifecycle out of Backend into `communication/`; update Context publishing
- [ ] Move CPU thread/timer and CUDA stream ownership into `device/`; analysis executor into `asynchronous/`
- [ ] Implement thin plan dispatch + component lifecycle in `driver/` (on Driver or a small helper type); delete Backend
- [ ] Rewrite all tests/examples that constructed Backend; cuda `test_backend.cu` becomes device/driver-focused tests under new paths
- [ ] Add regression `regressions/layout/no_backend/` — MD or Evaluate run without naming Backend
- [ ] Run full check + test suite

## Testing strategy

- Unit tests for communication bind and driver step without Backend type.
- Prior no-alloc / role-validation tests still pass under new wiring.
- Regression is public preset API only.

## Out of scope

- New MPI algorithms
- Directory renames
- component kinds/ layout tidy
- Docs (07)
