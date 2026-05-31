---
title: Driver device unification — link 01 — unified Backend staging interface
slug: driver-device-unification-01-backend-staging-interface
status: approved
created: 2026-06-01
parent: driver-device-unification (design)
chain: driver-device-unification-* (link 1/4)
---

# Driver device unification — link 01 — unified Backend staging interface

## Summary

链 1/4 的**纯加法地基**：给 `cpu::Backend` 与 `cuda::Backend` 各加一组语义统一的设备暂存方法 `stage_in(cpu::Frame&)` / `stage_out(cpu::Frame&)` / `drain()`，把两个 backend 之间唯一真正的设备差异（CPU `set_frame`/无/`drain_analysis` vs CUDA `materialize`/`retrieve`/`sync`）收敛到同一组接口名后面；并给 `cpu::ComponentBase` 补上 `provides_neighbor_list()` 虚函数（CUDA 端早已有，spec `cuda-md-review-cleanup-01` 落地），`cpu::md::NeighborList` 覆写返回 `true`。

本 link **不动任何 driver、不删任何旧名**：`stage_*` 是对现有方法的薄 wrapper，旧公开名保留；driver 层照旧调用旧名。零行为变化，全测原样绿。后续 link 02（`BackendDriver<P,DT>`）才开始消费 `stage_*` 与 `provides_neighbor_list()`，link 03 再移除旧公开名。

## Domain basis

契约源：`CLAUDE.md ## Architecture § Data residency`（CPU Frame 零拷贝、CUDA Frame 经 `cpu::to_cuda`/`cuda::to_cpu` 显式 H2D/D2H）+ `§ Component tagging`（typed `role()`/`capabilities()`/`reads()` 查询）。

既有事实（已核对当前 worktree HEAD）：
- `cpu::Backend`（`include/atomiverse/cpu/backend.h`）：`set_frame(cpu::Frame&)`（零拷贝绑定 + DomainDecomposition bind）、无帧出、`drain_analysis()`、`shutdown()`。
- `cuda::Backend`（`include/atomiverse/cuda/backend.cuh`，实现在 `src/cuda/driver.cu`）：`materialize(const cpu::Frame&)`（H2D）、`retrieve(cpu::Frame&)`（D2H）、`sync()`。
- `cuda::ComponentBase` 已有 `reads()` + `needs_neighbor_list()` + `provides_neighbor_list()`；`cpu::ComponentBase`（`include/atomiverse/cpu/component.h:88,99`）有前两个但**缺** `provides_neighbor_list()`。
- CPU `scan_pending_infra_`（`cpu/cpu_driver.h`）仍用 `name()=="neighbor_list"` 探测 NeighborList —— 本 link 只补虚函数，**不改 scan**（改 scan 属 link 02 的 typed-query 统一）。

## Design

### 统一暂存接口（语义映射，不改行为）

| 统一方法 | `cpu::Backend` 实现 | `cuda::Backend` 实现 |
|---|---|---|
| `void stage_in(cpu::Frame&)` | 调 `set_frame(frame)` | 调 `materialize(frame)` |
| `void stage_out(cpu::Frame&)` | no-op（frame 即存储，零拷贝） | 调 `retrieve(frame)` |
| `void drain()` | 调 `drain_analysis()` | 调 `sync()` |

- `cpu::Backend::stage_in` 取非 const `cpu::Frame&`（CPU 需要绑定 DD，写入 frame）；`cuda::Backend::stage_in` 也取非 const `cpu::Frame&`（内部按 const 处理，签名对齐以便 link 02 的统一模板能用同一调用）。
- 旧名 `set_frame`/`materialize`/`retrieve`/`drain_analysis`/`sync` **保持 public 不变**（link 03 才移除），`stage_*` 仅委托。

### ComponentBase typed-query 对齐

`include/atomiverse/cpu/component.h`：
```cpp
/// @brief Does this component supply the neighbor pair list?
/// @return true for NeighborList; default false.
virtual bool provides_neighbor_list() const { return false; }
```
`cpu::md::NeighborList`（`include/atomiverse/cpu/md/neighbor_list.h`）覆写返回 `true`，与 `cuda::NeighborList` 对齐。

## Files to create or modify

**Modify:**
- `include/atomiverse/cpu/backend.h` — 声明 `stage_in/stage_out/drain`
- `src/cpu/backend.cpp` — 实现（委托 `set_frame` / no-op / `drain_analysis`）
- `include/atomiverse/cuda/backend.cuh` — 声明 `stage_in/stage_out/drain`
- `src/cuda/driver.cu` — 实现（委托 `materialize` / `retrieve` / `sync`）
- `include/atomiverse/cpu/component.h` — 加 `provides_neighbor_list()` 虚函数
- `include/atomiverse/cpu/md/neighbor_list.h` — 覆写 `provides_neighbor_list()==true`

**Create:**
- `tests/cpu/test_backend_staging.cpp` — 验证 `stage_in/stage_out/drain` 与 `provides_neighbor_list()`（CPU 侧；CUDA 侧 staging 已被现有 cuda_md 端到端测试覆盖，本 link 不新增 .cu）

## Tasks

- [ ] Write failing test `tests/cpu/test_backend_staging.cpp`: (a) `cpu::Backend` exposes `stage_in/stage_out/drain` and a build→stage_in→step→stage_out→drain sequence on a tiny LJ frame matches the pre-existing `set_frame`-based path bit-for-bit; (b) `cpu::md::NeighborList{}.provides_neighbor_list()==true` and a bare `cpu::ComponentBase`-derived stub returns `false`
- [ ] Add `void provides_neighbor_list()` virtual (default `false`) to `cpu::ComponentBase` (`include/atomiverse/cpu/component.h`)
- [ ] Override `provides_neighbor_list()==true` on `cpu::md::NeighborList` (`include/atomiverse/cpu/md/neighbor_list.h`)
- [ ] Add `stage_in/stage_out/drain` to `cpu::Backend` (`backend.h` + `src/cpu/backend.cpp`) delegating to `set_frame` / no-op / `drain_analysis`
- [ ] Add `stage_in/stage_out/drain` to `cuda::Backend` (`backend.cuh` + `src/cuda/driver.cu`) delegating to `materialize` / `retrieve` / `sync`
- [ ] Build CPU + CUDA (`cmake --build build -j$(nproc)`) clean; run full `ctest` — all green, zero behavior change
- [ ] Hygiene cleanup: confirm `/mol:simplify` runs clean on touched files

## Testing strategy

1. **Build** — `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON && cmake --build build -j$(nproc)` 干净（CPU + CUDA TU 都编译新方法）。
2. **New staging test** — `ctest -R test_backend_staging` 绿：CPU `stage_*` 路径与旧 `set_frame` 路径等价；`provides_neighbor_list()` 真值正确。
3. **Full regression** — `ctest --test-dir build --output-on-failure` 全绿（`test_md`、`cuda_md*`、minimizer、evaluate、vib/hessian 等无行为变化）。
4. **Additive guarantee** — 旧公开名仍在、driver 调用点未改，确认本 link 不破坏任何现有调用。

## Out of scope

- `BackendDriver<P,DT>` 与 `backend_traits` glue 类型 → link 02
- 移除旧 backend 公开名（`set_frame`/`materialize`/`retrieve`/`drain_analysis`/`sync`）→ link 03
- 重写 `scan_pending_infra_` 用 typed-query / 删 CPU 字符串探测 → link 02
- category 模板合并、删 driver 头、reparent → link 03
- parity 测试 + CLAUDE.md taxonomy 图 → link 04
