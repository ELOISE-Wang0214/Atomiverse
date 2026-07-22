---
spec: backend-mirror-04-driver
status: approved
created: 2026-07-22
chain: backend-mirror
phase: 04-driver
depends_on: backend-mirror-03-paths
grilled: true
---

# backend-mirror-04-driver — Driver 命名与 Criterion 统一

## Summary

去掉类名上的后端前缀：循环节点为 `cpu::Driver` / `cuda::Driver`；原非模板可嵌套擦除基类改名为 `DriverNode` 以消除重名。CUDA 终止策略与 CPU 对齐为 **Criterion** + `add_criterion`，删除公共 `LoopHandler` / `add_handler`。**无 typedef/using 兼容别名**；测试、绑定、示例按新符号重写。本环 **不** 拆 Backend（见 05）。

## Domain basis

无新物理。Verdict 合成顺序保持 Stop > Converged > Continue。

## Design

| 现状 | 目标 |
|------|------|
| `cpu::CPUDriver<P>` | `cpu::Driver<P>` |
| `cuda::CUDADriver<P>` | `cuda::Driver<P>` |
| `cpu::Driver`（擦除基类） | `cpu::DriverNode` |
| CUDA 对等擦除基类（若有） | `cuda::DriverNode` |
| `cuda::LoopHandler`, `add_handler` | `cuda::Criterion`, `add_criterion`，与 `cpu::Criterion` 同构的 `verdict`/`name` 面（CUDA 可仍非 Component 若工程代价过高，但 **用户可见 API 名称必须一致**） |
| 头文件 `cpu_driver.h` / `cuda_driver.cuh` / `loop_handler.cuh` | 改为 `driver.h` / `criterion.cuh`（或 `component/criterion.cuh` 与 CPU 对称） |

**Backend 类本环可暂留**，但不得再扩展；05 删除。

### Reuse decision

| Candidate | Decision |
|-----------|----------|
| CPUDriver/CUDADriver 实现体 | **generalize（rename）** → `Driver` |
| 擦除基类 `cpu::Driver` | **generalize（rename）** → `DriverNode` |
| `cpu::Criterion` | **reuse** — CUDA 向其对齐 |
| `cuda::LoopHandler` | **generalize** → Criterion API |
| `criteria.h` / ForceConverged | **reuse** — CUDA 侧同名类 |

## Files to create or modify

- `include/atomiverse/cpu/driver/cpu_driver.h` → `driver` 实现头（命名以 Design 为准）
- `include/atomiverse/cpu/driver/driver.h`（DriverNode）
- `include/atomiverse/cuda/driver/cuda_driver.cuh` → `Driver`
- `include/atomiverse/cuda/driver/loop_handler.cuh` → Criterion
- All presets (dynamics MD, minimize, montecarlo, evaluate, path, vibration Hessian)
- tests/cpu/driver/*, tests/cuda/*, python bindings, examples, regressions

## Tasks

- [ ] Write failing scrubs/tests: product code has zero `CPUDriver`, `CUDADriver`, `LoopHandler`, `add_handler`; `DriverNode` and `add_criterion` appear on both backends
- [ ] Rename erased base to `DriverNode` and phase loop type to `Driver` on CPU; update presets/tests
- [ ] Rename CUDA loop type to `Driver`; replace LoopHandler with Criterion/add_criterion; update cuda tests
- [ ] Rewrite examples, regressions, python bindings for new symbols (no alias headers)
- [ ] Add regression `regressions/layout/driver_names/` — uses `cpu::Driver`-based evaluate/minimize surface + `add_criterion` only
- [ ] Run full check + test suite

## Testing strategy

- Grep-based architecture tests under `tests/architecture/` or `tests/core/`.
- Existing role-validation and minimize criterion tests rewritten to new names; behavior unchanged.
- No e2e under tests/.

## Out of scope

- Backend dissolve (05)
- Directory renames (03)
- Docs harness full rewrite (07)
