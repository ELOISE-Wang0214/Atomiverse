---
title: MSU PE CUDA components — link 07 — async correctness for CUDA recorders
status: approved
created: 2026-05-23
parent: msu-pe-validation.md (Phase H.1–H.3)
chain: msu-pe-cuda-components-* (link 7/8)
---

# MSU PE CUDA components — link 07 — async correctness for CUDA recorders

## Summary

落地 Phase H 异步契约的 **CUDA 端**:把链节 06 引入的 `CudaCsvRecorder`(以及任何未来的 CUDA-side `Record` / `Thermo` / `Observer` 组件)固化在 `cuda::AsyncRecordBase` 模板基类上,该基类强制 **GPU 必须在 `defer_stream` 上做一次 non-blocking D2H 拷贝**(`cudaMemcpyAsync` + cudaEvent),host 回调通过 CPU 端 `AnalysisExecutor` 异步执行 —— hot-path 永远不阻塞在 `cudaStreamSynchronize` 或 D2H 拷贝上。本链节交付:(a) `cuda::AsyncRecordBase` 模板与 host 回调 dispatcher,(b) 改造链节 06 的 `CudaCsvRecorder` 继承新基类,(c) 一组验证测试:hot-path 上 0 次 `cudaStreamSynchronize` 调用、recorder 频率(每 25 步 vs 每 1 步)对每步 wall-clock 的影响 < 5%、最终 CSV 与同步 baseline 字段一致(行可能丢失但已写入的字节相同)。该链节是父规范 `msu-pe-validation` Phase H.1–H.3 的硬正确性 + 性能门。

## Domain basis

CUDA 异步契约源自 `core/async/bridge.h` 已有的 `AsyncBridge` + `defer_stream` 模式(MD 引擎暴露一个专用流 `cuda::md::MD::async_stream()`,所有 non-hot-path D2H 拷贝必须用这条流而不是默认 compute stream)。`cudaEvent_t` 用于 host 端轮询某个 D2H 是否完成 —— host poll 在 CPU `AnalysisExecutor` 的 worker thread 上跑,主线程从不阻塞。这是 NVIDIA 标准 single-GPU async pattern(《CUDA C++ Programming Guide》§ 3.2.5.6 Concurrent kernel execution and data transfers;CUDA Toolkit 12 docs)。

CPU 端先例:`cpu::AsyncComponentBase` (`include/atomiverse/core/async/component_base.h`,已实现并由 cpu::CsvRecorder 继承)定义了 `AsyncComponentBase::submit(work_fn)` → `AnalysisExecutor` 的契约。CUDA 端 `AsyncRecordBase` 与 CPU 同形,但 `work_fn` 需要先等 `cudaEvent` 完成才能读 host-side snapshot —— `AnalysisExecutor` 在 worker thread 上 `cudaEventSynchronize` 是合法的(它在 worker thread,不在 hot-path)。

性能基准:MSU PE 链节 06 的 `--scale 0.01` smoke run(stage 5 约 1718 步)在 GH200 上 wall-clock ~ 8 秒;按每步 ~ 4.7 ms 算,recorder 频率从 25 步一次提升到 1 步一次必须保持 wall-clock 差 < 5%(否则证明 hot-path 有阻塞)。

单位:wall-clock [ms / step];`cudaStreamSynchronize` 调用计数 [counts];snapshot 字段与 CSV 字段一一对应,无单位转换。

## Design

**总体策略:把链节 06 的 `CudaCsvRecorder` ad-hoc 异步实现抽出来作为可复用基类。** 链节 06 已经在 design 里承诺"通过 `AsyncBridge` 在 `defer_stream` 上 D2H 一份 `SystemScalars`,host 回调里写入 CSV",但当时实现细节落在 `polyethylene_tension.cpp` 里;本链节把它**升格为框架契约**,放在 `include/atomiverse/cuda/async/record_base.cuh`,并加入测试守护。

**新增 `cuda::AsyncRecordBase` 模板基类:**

```cpp
// include/atomiverse/cuda/async/record_base.cuh
namespace ATV::cuda {
template <typename Snapshot>
class AsyncRecordBase : public cuda::ComponentBase {
public:
    /// Owned device-side staging buffer (sized at build()).
    /// Concrete derived classes fill this on the compute stream
    /// before calling submit_async(...).
    Snapshot* device_snapshot_ = nullptr;
    Snapshot  host_snapshot_;
    cudaEvent_t copy_done_event_ = nullptr;
    cudaStream_t defer_stream_ = nullptr;

    /// Hot-path entry: enqueue D2H on defer_stream and submit
    /// host-side work to AnalysisExecutor. Returns immediately;
    /// never calls cudaStreamSynchronize.
    void submit_async(int64_t step,
                      std::function<void(const Snapshot&, int64_t)> work_fn);

    /// Lifecycle.
    void build(BuildContext& ctx) override;
    void destroy(BuildContext& ctx) override;
};
}  // namespace ATV::cuda
```

**`submit_async` 内部步骤(零阻塞):**

1. `cudaMemcpyAsync(&host_snapshot_, device_snapshot_, sizeof(Snapshot), cudaMemcpyDeviceToHost, defer_stream_);`
2. `cudaEventRecord(copy_done_event_, defer_stream_);`
3. 把闭包 `[this, step, work_fn] { cudaEventSynchronize(copy_done_event_); work_fn(host_snapshot_, step); }` 提交给 `AnalysisExecutor::instance().submit(...)`。
4. 立刻返回。`cudaEventSynchronize` 在 worker thread 跑,主线程 (hot-path) 不阻塞。

**`device_snapshot_` 在 `build()` 阶段 `cudaMalloc`,`destroy()` 阶段 `cudaFree` + `cudaEventDestroy`。`defer_stream_` 从 `ctx.md.async_stream()` 拿到(read-only,non-owning)。

**链节 06 的 `CudaCsvRecorder` 改造:**

- 继承 `cuda::AsyncRecordBase<SystemScalarsSnapshot>`,其中 `SystemScalarsSnapshot` 是个 plain struct,字段对应 CSV 列(`step`, `stage`, `T_target`, `T`, `P_bar`, `Pxx`, `Pyy`, `Pzz`, `E_pair`, `E_bond`, `E_angle`, `E_dihedral`, `E_pot`, `E_kin`, `box_x`, `box_y`, `box_z`)。
- `launch(ctx)` 内:在 compute stream 上跑一个轻量 `pack_snapshot<<<1, 32>>>(ctx.frame.scalars, this->device_snapshot_)` kernel,把分散在 device `SystemScalars` 字段里的标量打包到一个连续的 `SystemScalarsSnapshot` 里 → 然后调 `this->submit_async(step, [this](const auto& s, int64_t step) { write_csv_row(s, step); });`。
- `write_csv_row(...)` 是纯 host 函数,跑在 AnalysisExecutor worker thread,无 CUDA 调用。

**Hot-path 守护 — `cudaStreamSynchronize` 0 调用契约:**

新增内部 debug 宏 `ATV_CUDA_TRACK_SYNCS`(只在 `-DATV_BUILD_TESTS=ON` 时启用,默认 off,对 release build 零开销):

- 提供 wrapper `atv::cuda::stream_synchronize(stream)`,内部 atomic 计数器 `n_stream_syncs_during_md_step_`。
- `cuda::md::MD::step()` 在进入 hot-path 前清零计数,在 hot-path 结束(`final_integrate` 之后、`end_of_step` 之前)断言计数 == 0(在 `ATV_BUILD_TESTS=ON` 时 EXPECT,否则 release-build 无开销)。
- 任何在 hot-path 调用 `cudaStreamSynchronize` 的代码(误用旧式同步 D2H)立即让测试 fail。

**性能守护测试:**

新增 `tests/cuda/async/test_async_record_perf.cu`:

- 跑 `polyethylene_tension --backend cuda --mode stage:2 --scale 0.005 --print-every 25`(对照),记录 wall-clock。
- 跑 `polyethylene_tension --backend cuda --mode stage:2 --scale 0.005 --print-every 1`(高频),记录 wall-clock。
- 断言 (`wall_high / wall_low − 1) < 0.05`(5% 上限)。
- 用 `std::chrono::steady_clock` host-side 测量,跑 5 遍取中位数以稳定测量。

**正确性守护测试:**

新增 `tests/cuda/async/test_async_record_correctness.cu`:

- 跑两次 `polyethylene_tension --backend cuda --mode stage:2 --scale 0.001` 同一 seed:
  - Run A: 当前 async path(`AsyncRecordBase`)。
  - Run B: 临时切换到 sync path(测试-only 编译 flag `ATV_TEST_RECORD_SYNC=1`,通过 `cudaMemcpy` 而不是 `cudaMemcpyAsync` + event,等价于强制同步)。
- 比较两份 `thermo.csv`:数值字段逐行 `EXPECT_NEAR(±1e-6)`(单精度 round-trip 容差)。Async path 可能 *丢失* CSV 最末若干行(进程在 destroy 前 join,但万一 SIGINT 中断时);测试只断言"已写入的行内容一致",不断言行数一致(行数容差 ≤ 1)。

**API 表面:**

| 符号 | 文件 | 注 |
|------|------|------|
| `cuda::AsyncRecordBase<Snapshot>` | `include/atomiverse/cuda/async/record_base.cuh` | 新增模板基类 |
| `cuda::AsyncRecordBase::submit_async(step, work_fn)` | 同上 | 唯一公开入口 |
| `cuda::AsyncRecordBase::device_snapshot_` | 同上 | protected,子类 launch() 直接写 |
| `cuda::stream_synchronize(stream)` | `include/atomiverse/cuda/async/sync_guard.cuh`(new) | hot-path sync 计数 wrapper |
| `ATV_CUDA_TRACK_SYNCS` macro | 同上 | 测试-only 计数开关 |
| `CudaCsvRecorder`(链节 06) | `benchmarks/md/polyethylene_tension.cpp` | 改成继承 `AsyncRecordBase` |

**显式不做的事:**

- 不动 CPU 端 `cpu::AsyncComponentBase`(已经存在,链节 07 只是 CUDA 端对位实现)。
- 不实现 `cuda::AsyncThermoBase` / `cuda::AsyncObserverBase`(本链只交付 `AsyncRecordBase`;若链 08 或 Phase H 后续需要 Thermo/Observer,通过本基类的 template + duck typing 复用即可,新基类是 trivial fork)。
- 不动 `AnalysisExecutor` 实现(已存在于 `core/async/executor.h`;本链只是消费它)。
- 不实现多流(multi-stream)调度优化 —— `defer_stream` 复用 `cuda::md::MD::async_stream()` 已有的单条流。

## Files to create or modify

- `/cfs/klemming/home/j/jichenli/work/Atomiverse/include/atomiverse/cuda/async/record_base.cuh` (new) — `AsyncRecordBase<Snapshot>` template + protected staging buffer + submit_async.
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/src/cuda/async/record_base.cu` (new) — explicit template instantiation for `SystemScalarsSnapshot`,与 build/destroy 默认实现。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/include/atomiverse/cuda/async/sync_guard.cuh` (new) — `stream_synchronize` wrapper + `ATV_CUDA_TRACK_SYNCS` 计数。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/src/cuda/md/md.cu` — 在 `step()` 进入/结束处加 `ATV_CUDA_TRACK_SYNCS` reset + assert(仅 `-DATV_BUILD_TESTS=ON` 时插入,否则 macro 展为空)。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/benchmarks/md/polyethylene_tension.cpp` — `CudaCsvRecorder` 改继承 `cuda::AsyncRecordBase<SystemScalarsSnapshot>`,删除链节 06 留下的 ad-hoc 异步代码。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tests/cuda/async/test_async_record_correctness.cu` (new) — async vs sync recorder 字段比对测试。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tests/cuda/async/test_async_record_perf.cu` (new) — 25-step vs 1-step 频率的 wall-clock 守护测试。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tests/cuda/async/test_hot_path_no_sync.cu` (new) — 跑一个简短 cuda::md::MD 循环,断言 hot-path `cudaStreamSynchronize` 调用计数 == 0。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tests/cuda/async/CMakeLists.txt` (new) — 注册三个新测试,labels `cuda;async;gpu_required`。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tests/CMakeLists.txt` — `if(ATV_BUILD_CUDA)` 块下 `add_subdirectory(cuda/async)`。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/.claude/specs/INDEX.md` — 加 `msu-pe-cuda-components-07-async-correctness`。

## Tasks

- [ ] Write failing test `tests/cuda/async/test_hot_path_no_sync.cu` asserting hot-path `cudaStreamSynchronize` count == 0 over a 100-step `cuda::md::MD` loop with `CudaCsvRecorder` bound
- [ ] Write failing test `tests/cuda/async/test_async_record_correctness.cu` comparing async vs sync recorder CSV field-by-field (`EXPECT_NEAR ±1e-6`, row-count tolerance ≤ 1)
- [ ] Write failing test `tests/cuda/async/test_async_record_perf.cu` asserting wall-clock(`--print-every 1`) / wall-clock(`--print-every 25`) − 1 < 0.05 on stage 2 of `polyethylene_tension --backend cuda --scale 0.005`
- [ ] Implement `cuda::AsyncRecordBase<Snapshot>` template in `include/atomiverse/cuda/async/record_base.cuh` + explicit `SystemScalarsSnapshot` instantiation in `src/cuda/async/record_base.cu`
- [ ] Implement `cuda::stream_synchronize` wrapper + `ATV_CUDA_TRACK_SYNCS` counter in `include/atomiverse/cuda/async/sync_guard.cuh`
- [ ] Wire hot-path sync counter reset + assertion into `src/cuda/md/md.cu::step()` under `ATV_BUILD_TESTS=ON` only
- [ ] Refactor `CudaCsvRecorder` in `benchmarks/md/polyethylene_tension.cpp` to inherit `cuda::AsyncRecordBase<SystemScalarsSnapshot>`, deleting the ad-hoc async code from chain link 06
- [ ] Wire `add_subdirectory(cuda/async)` into `tests/CMakeLists.txt` under `if(ATV_BUILD_CUDA)`
- [ ] Register the three new tests in `tests/cuda/async/CMakeLists.txt` with `LABELS "cuda;async;gpu_required"`
- [ ] Add `msu-pe-cuda-components-07-async-correctness` entry to `.claude/specs/INDEX.md`
- [ ] Add doxygen on `AsyncRecordBase` documenting the zero-block contract, the defer-stream + cudaEvent pattern, and the `AnalysisExecutor` host-thread semantics
- [ ] Run full check + test suite (`cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_CUDA=ON -DATV_BUILD_BENCHMARKS=ON && cmake --build build -j && ctest --test-dir build --output-on-failure -L cuda`)

## Testing strategy

- **Happy path — hot-path sync count == 0:** `test_hot_path_no_sync` 跑 100 步 `cuda::md::MD` + `CudaCsvRecorder`(`print_every = 1`),断言 `n_stream_syncs_during_md_step_ == 0` 在每步结束时。
- **Happy path — async vs sync field parity:** `test_async_record_correctness` 跑两次相同 seed,async path vs 强制 sync path(`-DATV_TEST_RECORD_SYNC=1` 编译开关)。CSV 字段逐行 `EXPECT_NEAR ±1e-6`,允许行数差 ≤ 1。
- **Performance gate — 25× recorder frequency:** `test_async_record_perf` 测 `--print-every 25` vs `--print-every 1` 的 wall-clock 比,断言 < 1.05;若超出说明 hot-path 仍被阻塞。
- **Edge — lifecycle:** `build → launch → destroy → build → launch → destroy` 循环 3 遍,`compute-sanitizer --tool memcheck` 报 0 leak;`device_snapshot_` / `copy_done_event_` 在 destroy 后为 `nullptr`。
- **Edge — SIGINT/early exit:** 在 stage 2 跑到一半发送 SIGINT,验证 CSV 文件不损坏(行尾完整,字段无截断)—— 这要求 `AnalysisExecutor` 在析构时 join 所有 worker thread,本链节继承现有契约,只验证不实现。
- **Domain validation:** 不适用 —— 本链节是基础设施 + 性能门,无新增物理。物理由链节 04(NoseHooverNPT)、链节 06(benchmark integration)、链节 08(cross-engine validation)守护。

## Out of scope

- `cuda::AsyncThermoBase` / `cuda::AsyncObserverBase` —— 父规范 Phase H 的扩展空间,但本链节只交付 `AsyncRecordBase`(MSU PE benchmark 当前只用 Record,不需要 Thermo/Observer 在 CUDA 端;CPU 端 Thermo/Observer 用 `cpu::AsyncComponentBase`)。
- CPU 端 `cpu::AsyncComponentBase` 改动 —— 已存在并工作。
- 多 GPU async patterns(`cudaStreamWaitEvent` 跨 device sync)—— 父规范 fixed scope decision #1:CUDA 单卡。
- `AnalysisExecutor` 实现层改动 —— 已存在于 `core/async/executor.h`,本链节只消费。
- 异步轨迹写入(`TrajectoryRecorder` 的 CUDA 版本)—— 本链节只覆盖 `CudaCsvRecorder`;trajectory 输出由后续链节(若需要)按相同模式扩展。
- Multi-stream graph capture / `cudaGraph` 优化 —— 性能优化,不是正确性门;本链节的 5% wall-clock 容差已经足够覆盖 MSU PE 验证需求。
