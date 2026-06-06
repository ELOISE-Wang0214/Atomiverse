---
title: CUDA device-side snapshot ring + compute-stream freeze (async-analysis link 2/3)
status: code-complete
created: 2026-06-06
---

# CUDA device-side snapshot ring + compute-stream freeze (async-analysis link 2/3)

## Summary

修复 CUDA 快照撕裂（tearing）缺陷：原 `cuda::AsyncBridge` 仅做单向跨流定序（`cudaEventRecord(main)` → `cudaStreamWaitEvent(defer)`），随后在 defer 流上从**活跃设备指针**做 D2H，而主流的第 N+1 步不会等待该 D2H 完成——main 与 defer 并发执行，构成跨流的 read-after-write 竞争，产生被错误标记为第 N 步的撕裂快照。修复把"在改写 frame 的计算/主流上把 payload D2D 冻结进设备环、再从冻结槽 D2H"直接做进 `AsyncBridge` 基类（freeze 排在第 N 步之后、第 N+1 步之前——这是修复核心），覆盖全部五个 `AsyncBridge` 子类；并给 `AsyncRecordBase` 加内部设备环 + drop_oldest 背压。**设计修正**（见 Design）：原计划的 `cuda::DevicePublish` + 共享 `cpu::AnalysisExecutor` + 主机 `FrameSnapshot` 环不适配 Thermo/RDF/Mulliken/Density（它们捕获标量/pair 距离/密度矩阵，非 FrameSnapshot 字段），故落点改为 `AsyncBridge` 基类，无需 `DevicePublish`/`Context.executor`/后端集成。冻结-再分发不变量与 link 01 一致。

## Domain basis

无新增物理（声明 `science.required=true`）。本链路的领域依据是 **GPU 路径上的快照一致性不变量**：对任一交付到主机消费者的设备环快照 S，存在唯一步号 n，使 S 中**全部**原子坐标满足 `x_i = f(tag_i, n)`，且 `S.step == n`。该不变量可证伪：构造解析轨迹 `x_i(step) = f(tag_i, step)`（例如 `x = tag_i + step`，每步在计算流上以确定 kernel 推进），则任何撕裂（部分原子为第 N 步、部分为第 N+1 步）或错标（坐标属第 N 步而 `step` 标为 N-1）的设备环快照都会在主机消费者侧被逐元素断言捕获。该判据对 recorder 的坐标路径与 thermo 的 `scalars->energy`（次步 `pre_force` 的 `ZeroEnergy` 会清零，构成竞争）分别成立。单位沿用项目约定：坐标 Å、能量 eV、步号无量纲。

## Design

> **设计修正（实现期，已在 GH200 验证）。** 原方案拟新增 `cuda::DevicePublish` 把所有组件经 `cpu::AnalysisExecutor` + 主机 `FrameSnapshot` 环统一。实现时发现该机制只适配 `Recorder`（位置）——而 `Thermo` 捕获 `scalars->energy/kinetic`、`RDF` 捕获 pair 距离、`Mulliken`/`DensityRecorder` 捕获 SCF 密度矩阵，这些都**不是** `FrameSnapshot` 字段，强行塞入需给 `FrameSnapshot` 挂任意 payload 或改变各自的计算口径。由于这五个组件都派生自 `AsyncBridge`，正确且最小的修复是把"计算流冻结进环"直接做进 `AsyncBridge` 基类——于是**无需** `DevicePublish`、**无需** `cuda::Context.executor`、**无需** CUDA 后端集成。冻结-再分发的不变量与 link 01 完全一致，只是落点在 `AsyncBridge` 而非新组件。

- **`AsyncBridge` 基类新增设备环冻结（撕裂修复核心）。** 每个子类自带一个 K 深设备环（K=`kRingDepth`=3，与 CPU 主机环默认一致）。`launch()` 改为：
  1. `slot = ring_head_++ % K`；
  2. `freeze_slot(ctx, slot)`：子类在**计算/主流**（`ctx.runtime.stream`）上把所需 payload D2D 冻结进自己的 `device_ring[slot]`——排在本步物理之后、下一步改写 frame 之前，这一步关闭跨流 read-after-write；
  3. `cudaEventRecord(pack_events_[slot], main)` → `cudaStreamWaitEvent(defer, pack_events_[slot])`；
  4. `d2h_slot(ctx, slot)`：子类在 defer 流把**冻结后的** `device_ring[slot]` 拷进 pinned host slot（绝不读活跃 frame）；
  5. `cudaLaunchHostFunc(defer, cuda_callback_)`，回调在 CUDA worker 线程上以 slot 为参调用 `on_data_ready(ctx, slot)`。
  `launch()` 不做任何 `cudaStreamSynchronize`。子类钩子由 `issue_d2h`/`on_data_ready(ctx)` 改为 `freeze_slot` + `d2h_slot` + `on_data_ready(ctx, slot)`。
- **五个 `AsyncBridge` 子类改造**（`recorder.cuh` 位置 / `thermo.cuh` 标量 / `observer/rdf.cuh` pair 距离 / `observer/mulliken.cuh` 密度 / `recorder/density_recorder.cuh` 密度）：各自分配 K 个设备环槽 + K 个 pinned host 槽，实现 `freeze_slot`（D2D 活跃→环槽，main 流）+ `d2h_slot`（环槽→pinned，defer 流）+ `on_data_ready(slot)`。CPU 侧消费逻辑（`cpu::thermo_print` / `cpu::trajectory_append` / `cpu::rdf_accumulate` / `cpu::mulliken_charge` / `cpu::eval_density_grid`）已共享，保持不变。
- **`AsyncRecordBase<Snapshot>` 改造**：`device_snapshot_`（单 buffer）保留为子类暂存目标（**子类 API 不变**），基类内部新增 `device_ring_[K]` + `host_ring_[K]` + 每槽 event；`submit_async` 把 `device_snapshot_` D2D 冻结进环槽再 D2H。**关键**：`HostWorker` 队列无界，热循环可领先消费者远超 K 步，故加 **drop_oldest 背压**——每槽 `slot_busy_` 标志，若目标槽仍在途则丢弃本次捕获（计入 `n_dropped()`），绝不阻塞热循环、绝不投递被覆盖的快照。（`AsyncBridge` 经 `cudaLaunchHostFunc` 天然回压主机入队线程，热循环不会领先 defer 超过 K 步，故无需显式丢弃。）

## Files to create or modify

- `include/atomiverse/cuda/async_bridge.cuh` — 基类设备环冻结 API：`kRingDepth`、`pack_events_[K]`、`ring_head_`，钩子 `freeze_slot`/`d2h_slot`/`on_data_ready(ctx,slot)`，Doxygen（含单位）。
- `src/cuda/async_bridge.cu` — `launch()` 重排为 freeze→event→wait→d2h→`cudaLaunchHostFunc(slot)`；`CallbackPayload` 带 slot。
- `include/atomiverse/cuda/recorder/recorder.cuh` — K 设备环 + pinned 槽，冻结/拷回位置。
- `include/atomiverse/cuda/recorder/thermo.cuh` — 冻结 `scalars->energy/kinetic`（先于次步 `ZeroEnergy`）。
- `include/atomiverse/cuda/recorder/density_recorder.cuh` — 冻结 SCF 密度矩阵（RHF-gated）。
- `include/atomiverse/cuda/observer/rdf.cuh` — 冻结 pair 距离/计数。
- `include/atomiverse/cuda/observer/mulliken.cuh` — 冻结密度矩阵 + 能量（RHF-gated）。
- `include/atomiverse/cuda/async/record_base.cuh` — 内部 `device_ring_[K]`+`host_ring_[K]`+每槽 event + `slot_busy_` drop_oldest + `n_dropped()`（子类 API 不变）。
- `tests/cuda/async/test_async_bridge_tearing.cu` (new) — 解析轨迹撕裂回归（AsyncBridge）+ AsyncRecordBase no-in-flight-overwrite + 0-sync 断言。
- `tests/cuda/async/CMakeLists.txt` — 注册 benchmark-independent 测试（hot_path + tearing）无条件，benchmark-dependent 测试内部 gate。
- `tests/CMakeLists.txt` — `add_subdirectory(cuda/async)` 不再被 `ATV_BUILD_BENCHMARKS` 全量 gate。

## Tasks

- [x] Write CUDA tearing regression over analytic trajectory x_i=tag_i+step covering recorder position path (tests/cuda/async/test_async_bridge_tearing.cu)
- [x] Write AsyncRecordBase no-in-flight-overwrite test (analytic per-step value + lagging worker) (tests/cuda/async/test_async_bridge_tearing.cu)
- [x] Assert the device-ring freeze path adds ZERO cuda::stream_synchronize on the step loop (tearing test + test_hot_path_no_sync.cu)
- [x] Implement AsyncBridge device-ring freeze: freeze_slot/d2h_slot/on_data_ready(slot) + per-slot pack events + launch reorder (include/atomiverse/cuda/async_bridge.cuh, src/cuda/async_bridge.cu)
- [x] Retrofit the 5 AsyncBridge subclasses to freeze on the compute stream and D2H the frozen slot (recorder.cuh, thermo.cuh, density_recorder.cuh, observer/rdf.cuh, observer/mulliken.cuh)
- [x] Retrofit AsyncRecordBase to an internal device_ring_[K]+host_ring_[K] with drop_oldest slot-busy backpressure (include/atomiverse/cuda/async/record_base.cuh)
- [x] Un-gate the benchmark-independent async tests so the tearing + zero-sync gates build without benchmarks (tests/CMakeLists.txt, tests/cuda/async/CMakeLists.txt)
- [x] Add Doxygen with units to the AsyncBridge ring API, the subclass ring buffers, and the AsyncRecordBase ring
- [x] Verify tearing + no-overwrite + zero-sync on the GH200 node (jobs 9586/9590; ATV_BUILD_CUDA=ON)
- [x] Run full CUDA ctest on the GH200 node: normal config 102/102 green (job 9609) + RHF=ON build green / retrofits compile (job 9613). Benchmarks-ON path deferred — bench_scf/bench_rhf have a pre-existing unrelated cpu::Frame::energy break (→ /mol:fix)

## Testing strategy

- **Happy path**：每步捕获、解析轨迹 `x_i(step)=f(tag_i,step)`，慢速主机 worker；断言每个交付的设备环快照均为单步且 `step` 标号正确。
- **Edge cases**：(a) back-to-back 每步捕获时设备环 `device_ring[K]` 不覆盖 in-flight 快照；(b) thermo `scalars->energy` 在次步 `ZeroEnergy` 清零前已被冻结进环（专项覆盖该竞争）；(c) recorder 坐标路径与 thermo 标量路径各自独立断言。
- **Hot-path 不变量**：`DevicePublish`/设备环路径在 step loop 上 `cuda::stream_synchronize` 计数为 0（沿用 `sync_guard.cuh` 计数器）。
- **结构性检查**：`issue_d2h` 不再出现对 `ctx.frame.atoms.*` / `frame.scalars->energy` / `pairs->*` / `workspace().P` 等活跃指针的读取（grep 断言）。
- **Domain validation**：快照一致性不变量在 GH200 上经 CUDA ctest 验证；CUDA 构建与测试为 GH200 节点上的 operator-action（`ATV_BUILD_CUDA=ON`，见 CLAUDE.md GPU 节点说明），登录节点无法执行。

## Out of scope

- link 01（简化版 `AnalysisExecutor`、静态 mask、共享主机 `FrameSnapshot` 环）与 link 03（最终统一/收尾）的内容；本链路假定 link 01 已落地。
- `AsyncBridge` 基类本身的废弃/删除（终态收敛由 link 03 处理）；本链路只在基类内加设备环冻结并改造五个子类。
- 任何 `block` / `grow` 背压策略或 per-component 背压——本链路仅 drop_oldest + 全局 K（已与备选方案对比后排除，理由是与 CPU 路径对称且 hot-path 最廉价）。
- CPU 侧消费函数（`cpu::thermo_print` / `trajectory_append` / `rdf_accumulate` / `mulliken_charge`）的任何改动——它们已共享且保持不变。
