---
title: Driver device unification — link 04 — cross-backend parity tests + taxonomy docs
slug: driver-device-unification-04-parity-tests-and-docs
status: draft
created: 2026-06-01
parent: driver-device-unification (design)
chain: driver-device-unification-* (link 4/4)
---

# Driver device unification — link 04 — cross-backend parity tests + taxonomy docs

## Summary

链 4/4：用测试锁死「设备无关 driver 骨架」这一核心命题，并更新文档。新增 **Plan-equivalence parity**（同一组 categorical setter 调用，CPU 与 CUDA driver `build()` 后产出的 Plan 的 (phase, order, role, capabilities) 序列一致）与 **validation parity**（非法配置在两后端抛同类异常）；为 CUDA category 模板加 `molrs::FrameRef` 往返 smoke；更新 `CLAUDE.md` 的 driver-taxonomy 继承图为单模板栈；跑 umbrella `/mol:review` 收口。

## Domain basis

契约源：`CLAUDE.md ## Architecture § Driver taxonomy` + `§ Conventions`（Driver 模板范式）。本 link 的 parity 测试是 chain 的验收闸——证明 link 02/03 的合并确实做到「写一次、双设备等价」，而非引入隐性漂移。

前提（link 03 已落地）：四个 `*Driver<DT>` 模板 + reparent 完成；`stage_*` 为唯一公开暂存接口。

## Design

- `tests/test_driver_device_parity.cpp`（+ `.cu` 镜像，门控 `ATV_BUILD_CUDA`/`-DATV_BUILD_RHF`）：
  - **Plan-equivalence**：用相同的 `set_pair/set_neighbor/set_integrator/...` 序列分别构造 CPU 与 CUDA `md::MD`（及 `Evaluate`/`Minimize`），`build()` 后比较 `plan()` 的 (phase, order, role, capabilities) 列表逐项一致。
  - **Validation parity**：缺 Integrator / 缺 Force / 缺 Evaluator(+capabilities) / 缺 NeighborList 四种非法配置，两后端均在 `build()` 抛同类异常（同一份 `validate_components()`/`install_defaults()`）。
- `tests/cuda/test_cuda_molrs_overload.cu`：materialize 一个 molrs.Frame → `cuda::Evaluate.run(molrs::FrameRef&)` → write_back，验证 CUDA 经模板新得的 molrs 重载可编译/链接/往返。
- `CLAUDE.md` taxonomy 图改为：`Driver<P,DT>` → `BackendDriver<P,DT>` → `{Evaluate,Dynamics,Optimize,Response}Driver<DT>`，删 `cpu::CPUDriver`/`cuda::CUDADriver` 描述。

## Files to create or modify

**Create:**
- `tests/test_driver_device_parity.cpp`
- `tests/cuda/test_driver_device_parity.cu`
- `tests/cuda/test_cuda_molrs_overload.cu`

**Modify:**
- `CLAUDE.md` — ## Architecture § Driver taxonomy 继承图 + § Implemented components 表
- `.claude/NOTES.md`（如适用）— 记录单模板栈的迁移决策

## Tasks

- [ ] Write `tests/test_driver_device_parity.cpp` (+ `.cu`): Plan-equivalence across CPU/CUDA from identical setter calls; assert run RED first (before any fix needed — these are new guards, expected GREEN once link 03 is in)
- [ ] Add validation-parity assertions (missing Integrator/Force/Evaluator/NeighborList throw the same category on both backends)
- [ ] Add `tests/cuda/test_cuda_molrs_overload.cu` — CUDA molrs.Frame round-trip via the template overload
- [ ] Update `CLAUDE.md` driver-taxonomy inheritance diagram to the single-template stack; refresh the Implemented-components driver table
- [ ] Build CPU+CUDA clean; full ctest green including the new parity tests
- [ ] Run umbrella `/mol:review` on the full chain diff; record verdict

## Testing strategy

1. **Build** CPU+CUDA clean.
2. **Parity** — `test_driver_device_parity` green: identical Plans + identical validation throws across backends.
3. **molrs CUDA** — `test_cuda_molrs_overload` green.
4. **Full regression** — entire ctest suite green.
5. **Review gate** — `/mol:review` on the chain: 0 🚨, minimal 🔴.

## Out of scope

- Performance benchmarking of the unified stack (no hot-path change expected) — defer to existing bench suite.
- Merging concrete drivers into `md::MD<DT>` single templates → future.
- CUDA Response/Hessian → future.
