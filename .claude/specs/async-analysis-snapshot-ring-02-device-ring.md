---
title: CUDA device-side snapshot ring + compute-stream freeze (async-analysis link 2/3)
status: approved
created: 2026-06-06
---

# CUDA device-side snapshot ring + compute-stream freeze (async-analysis link 2/3)

## Summary

修复 CUDA 快照撕裂（tearing）缺陷：当前 `cuda::AsyncBridge` 仅做单向跨流定序（`cudaEventRecord(main)` → `cudaStreamWaitEvent(defer)`），随后在 defer 流上从**活跃设备指针**做 D2H，而主流的第 N+1 步不会等待该 D2H 完成——main 与 defer 并发执行，构成跨流的 read-after-write 竞争，产生被错误标记为第 N 步的撕裂快照。本链路新增 `cuda::DevicePublish`：在**会改写 frame 的计算/主流上**用一个 D2D pack kernel 把静态联合字段冻结进设备环 `device_ring[slot]`（在第 N 步之后、第 N+1 步之前完成——这是修复的核心），记录 `pack_done_event`，再在 defer 流上把该 slot 异步 D2H 进 link 01 已有的共享 CPU 主机 `FrameSnapshot` 环，并提交到共享 `AnalysisExecutor` 工作线程池。同时为 `AsyncRecordBase` 增加设备环（单 buffer → `device_ring[K]`），并把四个 `AsyncBridge` 子类改为消费冻结后的设备环 slot 而非活跃设备指针。终态与 CPU 路径对称：联合冻结进环 + drop_oldest + 全局 K + 共享 CPU 工作池。

## Domain basis

无新增物理（声明 `science.required=true`）。本链路的领域依据是 **GPU 路径上的快照一致性不变量**：对任一交付到主机消费者的设备环快照 S，存在唯一步号 n，使 S 中**全部**原子坐标满足 `x_i = f(tag_i, n)`，且 `S.step == n`。该不变量可证伪：构造解析轨迹 `x_i(step) = f(tag_i, step)`（例如 `x = tag_i + step`，每步在计算流上以确定 kernel 推进），则任何撕裂（部分原子为第 N 步、部分为第 N+1 步）或错标（坐标属第 N 步而 `step` 标为 N-1）的设备环快照都会在主机消费者侧被逐元素断言捕获。该判据对 recorder 的坐标路径与 thermo 的 `scalars->energy`（次步 `pre_force` 的 `ZeroEnergy` 会清零，构成竞争）分别成立。单位沿用项目约定：坐标 Å、能量 eV、步号无量纲。

## Design

- **新组件 `cuda::DevicePublish`**（镜像 `cpu::Publish`），绑定于 `md::Phase::end_of_step @ +10`（紧随每步归约、先于订阅者）：
  - `build()`：从 `ctx.runtime.executor`（共享 `AnalysisExecutor`，由 link 01 引入）解析主机环；`cudaMalloc` 设备环 `device_ring_`（K 个 `Snapshot` 槽）；创建 `pack_done_event_`、`copy_done_event_`（`cudaEventDisableTiming`，owning）；从 `ctx.runtime.stream` / `ctx.runtime.defer_stream` 捕获计算流与 defer 流（non-owning）。
  - `launch()`：① 在计算/主流上启动 D2D pack kernel，把静态联合 mask 字段冻结进 `device_ring_[slot]`（slot 由共享 executor 的 `acquire_slot()` 决定，drop_oldest，全局 K）；② 在计算流上 `cudaEventRecord(pack_done_event_)`；③ defer 流 `cudaStreamWaitEvent(pack_done_event_)` 后，把 `device_ring_[slot]` 异步 D2H 进主机环对应 slot；④ `cudaEventRecord(copy_done_event_, defer)`，提交 (订阅者, slot) 给共享工作池。`launch()` 不在 caller 线程做任何 `cudaStreamSynchronize`。
  - `destroy()`：释放 `device_ring_`、销毁两个 event、置空句柄。worker 池由共享 executor 负责 drain/join，不在此处。
- **设备环所有权与生命周期**：`device_ring_` 由 `DevicePublish` 拥有（`build` 分配 / `destroy` 释放）；主机环、worker 池复用 link 01 的共享设施，GPU 侧只新增设备环。环大小 K 取共享 executor 的 `ring_capacity()`，保证设备环与主机环一一对应、深度一致。
- **`cuda::Context` 新增 `executor` 句柄**：当前 `cuda::Context` 无 `executor` 字段；按 `cpu::Context` 镜像，新增 `cpu::AnalysisExecutor* executor = nullptr`，由后端在 init 时发布，使 `DevicePublish` / 改造后的订阅者可从任意 `BuildContext` 解析共享池。
- **`AsyncRecordBase<Snapshot>` 改造**：`device_snapshot_`（单 buffer，行 ~255）→ `device_ring_[K]`；保留其计算流 pack + `pack_done_event_` 定序；`submit_async` 按 round-robin 写入 `device_ring_[slot]`，避免 back-to-back 捕获覆盖仍在 D2H 中的 in-flight 快照；仅 drop_oldest、全局 K（无 per-component 背压）。
- **四个 `AsyncBridge` 子类改造**（`recorder.cuh` / `thermo.cuh` / `observer/rdf.cuh` / `observer/mulliken.cuh`）：`issue_d2h` 不再从活跃设备指针（`ctx.frame.atoms.x` / `&ctx.frame.scalars->energy` / `pairs->...` / `workspace().P`）读取，改为消费 `DevicePublish` 冻结后的设备环 slot；折叠到 `DevicePublish` + `AnalysisExecutor`。CPU 侧消费逻辑（`cpu::thermo_print` / `cpu::trajectory_append` / `cpu::rdf_accumulate` / `cpu::mulliken_charge`）已共享，保持不变。

## Files to create or modify

- `include/atomiverse/cuda/device_publish.cuh` (new) — `cuda::DevicePublish` 组件声明，设备环成员，Doxygen（含单位）。
- `src/cuda/device_publish.cu` (new) — D2D pack kernel + `build`/`launch`/`destroy` 实现。
- `include/atomiverse/cuda/context.cuh` — 新增 `cpu::AnalysisExecutor* executor` 字段（镜像 `cpu::Context`）。
- `include/atomiverse/cuda/async/record_base.cuh` — `device_snapshot_` → `device_ring_[K]`；`submit_async` round-robin。
- `src/cuda/async/record_base.cu` — 环分配/释放与 round-robin 实现配合。
- `include/atomiverse/cuda/recorder/recorder.cuh` — 消费设备环 slot 取代活跃 `frame.atoms.x/y/z`。
- `include/atomiverse/cuda/recorder/thermo.cuh` — 消费设备环 slot 取代活跃 `frame.scalars->energy/kinetic`。
- `include/atomiverse/cuda/observer/rdf.cuh` — 消费设备环 slot 取代活跃 `pairs->*` 指针。
- `include/atomiverse/cuda/observer/mulliken.cuh` — 消费设备环 slot 取代活跃 `workspace().P` / `scalars->energy`。
- `tests/cuda/async/test_async_record_correctness.cu` — 新增设备环撕裂回归用例（解析轨迹，覆盖坐标路径与 thermo 能量竞争）。
- `tests/cuda/async/test_hot_path_no_sync.cu` — 扩展覆盖 `DevicePublish`/设备环路径 0 sync。

## Tasks

- [ ] Write failing CUDA tearing regression for DevicePublish device-ring over analytic trajectory, covering recorder position path AND thermo scalars->energy race (tests/cuda/async/test_async_record_correctness.cu)
- [ ] Write failing back-to-back-capture test asserting AsyncRecordBase device_ring[K] never overwrites an in-flight snapshot (tests/cuda/async/test_async_record_correctness.cu)
- [ ] Extend failing hot-path no-sync test to assert DevicePublish device-ring path adds ZERO cuda::stream_synchronize on the step loop (tests/cuda/async/test_hot_path_no_sync.cu)
- [ ] Implement cuda::DevicePublish (compute-stream D2D pack kernel into device_ring + pack_done_event + defer-stream D2H into shared host ring + executor submit) in include/atomiverse/cuda/device_publish.cuh and src/cuda/device_publish.cu
- [ ] Add executor handle to cuda::Context in include/atomiverse/cuda/context.cuh (mirror cpu::Context)
- [ ] Retrofit AsyncRecordBase device_snapshot_ to device_ring_[K] with round-robin submit_async in include/atomiverse/cuda/async/record_base.cuh and src/cuda/async/record_base.cu
- [ ] Retrofit recorder.cuh, thermo.cuh, observer/rdf.cuh, observer/mulliken.cuh to consume the frozen device-ring slot instead of live device pointers
- [ ] Add Doxygen with units to DevicePublish, device_ring, and the pack kernel per project doc style
- [ ] Verify analytic-trajectory snapshot-consistency invariant + zero-sync on the GH200 node (operator-action; ATV_BUILD_CUDA=ON)
- [ ] Run full check + CUDA ctest suite on the GH200 node (operator-action)

## Testing strategy

- **Happy path**：每步捕获、解析轨迹 `x_i(step)=f(tag_i,step)`，慢速主机 worker；断言每个交付的设备环快照均为单步且 `step` 标号正确。
- **Edge cases**：(a) back-to-back 每步捕获时设备环 `device_ring[K]` 不覆盖 in-flight 快照；(b) thermo `scalars->energy` 在次步 `ZeroEnergy` 清零前已被冻结进环（专项覆盖该竞争）；(c) recorder 坐标路径与 thermo 标量路径各自独立断言。
- **Hot-path 不变量**：`DevicePublish`/设备环路径在 step loop 上 `cuda::stream_synchronize` 计数为 0（沿用 `sync_guard.cuh` 计数器）。
- **结构性检查**：`issue_d2h` 不再出现对 `ctx.frame.atoms.*` / `frame.scalars->energy` / `pairs->*` / `workspace().P` 等活跃指针的读取（grep 断言）。
- **Domain validation**：快照一致性不变量在 GH200 上经 CUDA ctest 验证；CUDA 构建与测试为 GH200 节点上的 operator-action（`ATV_BUILD_CUDA=ON`，见 CLAUDE.md GPU 节点说明），登录节点无法执行。

## Out of scope

- link 01（简化版 `AnalysisExecutor`、静态 mask、共享主机 `FrameSnapshot` 环）与 link 03（最终统一/收尾）的内容；本链路假定 link 01 已落地。
- `AsyncBridge` 基类本身的废弃/删除（终态收敛由 link 03 处理）；本链路仅把四个子类折叠到 `DevicePublish` + executor。
- 任何 `block` / `grow` 背压策略或 per-component 背压——本链路仅 drop_oldest + 全局 K（已与备选方案对比后排除，理由是与 CPU 路径对称且 hot-path 最廉价）。
- CPU 侧消费函数（`cpu::thermo_print` / `trajectory_append` / `rdf_accumulate` / `mulliken_charge`）的任何改动——它们已共享且保持不变。
