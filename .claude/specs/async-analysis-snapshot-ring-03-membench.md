---
title: Async snapshot ring — memory benchmark & scaling-law verification
status: code-complete
created: 2026-06-06
---

# Async snapshot ring — memory benchmark & scaling-law verification

## Summary
为异步快照环（async snapshot ring）增加一个常驻内存基准与其标度律的可证伪验证，作为 `async-analysis-snapshot-ring` 链条的第 3 环（收尾环）。基准在 GH200 节点上扫描环深 K、原子数 N 与捕获字段并集大小 S_union，分别用 `cudaMemGetInfo`（设备）与 `/proc/self/status` 的 VmHWM（主机）测量常驻内存，发出每个 (K,N) 的内存数据点，并在受控 worker 吞吐下用 `AnalysisExecutor::n_dropped()` 给出 drop-vs-K 的取舍曲线。基准对线性拟合残差设阈值，从而证伪或确认内存模型 mem ≈ K·N·S_union（CUDA 路径计两次）。同时把 freeze-then-dispatch 不变量、统一终态与设备环路径的 CUDA 执行模型分类落档到 `.claude/NOTES.md`。

## Domain basis
本环不引入新物理（science.required=true，特此声明：无新方程、无 DOI/arXiv 引用）。其"领域基础"是一个可证伪的线性内存标度模型：

- 常驻快照环内存 mem ≈ K · N · B
  - K = 设备环深度（每个 `AsyncBridge` 子类 / `AsyncRecordBase` 各持 K 个设备环槽 + K 个 pinned host 槽；生产值为编译期常量 `kRingDepth`=3）
  - N = 原子数（或快照元素数）
  - B = 每元素捕获字节数(e.g. positions = 3·sizeof(real_t))
- CUDA 路径同时持有设备环（`cudaMalloc`）与等大的 pinned host 环（`cudaMallocHost`），故内存被计两次：mem ≈ 2 · K · N · B（device 环 ≈ host 环）。
- 可证伪形式：基准自行按生产分配方式（`cudaMalloc`/`cudaMallocHost`）扫描 K 与 N，固定其余变量对 device 常驻内存做最小二乘线性拟合，残差须 ≤ 10%；否则模型被证伪。并报告生产配置 `kRingDepth`=3 的实际占用。

> **与原方案的对账（link 02 设计修正所致）。** 原 Domain basis 把 K 写成运行期 `AnalysisExecutor::Config.ring_size` 并假设 CUDA 侧共享 host 环。link 02 实为：每个 `AsyncBridge` 子类各自持 K 深设备环 + pinned host 环，K 是**编译期** `kRingDepth`=3，CUDA 侧无 `AnalysisExecutor`/共享 host 环。故本基准**自行**按生产分配方式扫描 K（验证标度律本身），报告生产 K=3 占用，而非读运行期全局 K；`S_union` 化为"每元素字节数 B"。drop-vs-K 曲线取消——K 为编译期常量，无运行期 drop 扫描。

## Design
- 新增基准可执行体 `bench_async_snapshot_mem`（plain `main()`，argv 驱动，发出 CSV，与 `bench_multi_gpu_scaling` 同风格——不链接 google-benchmark），位于 `benchmarks/bench_async_snapshot_mem.cu`。
- 测量：对一组 (K,N) 笛卡尔扫描（K∈{1,2,3,4,6}、N 多档），每点按生产分配方式分配 K 个设备环槽（`cudaMalloc`，每槽 N·B 字节）+ K 个 pinned host 槽（`cudaMallocHost`，等大）：
  - 设备常驻：`cudaMemGetInfo(&free,&total)` 在分配前后取 free 差值（可靠，free 时回退）。
  - 主机常驻：读 `/proc/self/status` 的 `VmRSS` 分配前后差值（pinned host 计入 RSS；用 VmRSS 而非单调的 VmHWM 以得到逐配置增量）。
  - 分配后立即释放，进入下一点。
- 输出：每个 (K,N) 一行 `K,N,bytes_per_elem,dev_bytes,host_bytes,model_bytes`。
- 标度律检查：基准内对 device 常驻做最小二乘线性拟合（固定 N 扫 K、固定 K 扫 N），打印两维残差 + device/host 比值（≈1，合计≈2×模型）；残差越界以非零退出，作为 operator 运行的 PASS/FAIL 门。报告生产 `kRingDepth`=3 的占用行。
- CMake：在 `benchmarks/CMakeLists.txt` 的 `if(ATV_BUILD_CUDA)` 块内用 `atv_add_bench(bench_async_snapshot_mem)` 注册（受 `ATV_BUILD_BENCHMARKS` + `ATV_BUILD_CUDA` 双重门控）。
- 文档：`.claude/NOTES.md` 追加 freeze-then-dispatch 不变量、link-02 终态（AsyncBridge 计算流设备环冻结 + AsyncRecordBase drop_oldest 环），并在 CUDA 执行模型分类表中新增 device-ring 路径条目。

> 注意:本基准为 **operator-action**——仅能在 aarch64 GH200 节点上以 `ATV_BUILD_CUDA=ON` + `ATV_BUILD_BENCHMARKS=ON` 构建运行（登录节点无 nvcc）。前置条件:`bench_scf.cu`/`bench_rhf.cu` 的 `cpu::Frame::energy` 旧 break 需先修(本环已顺带修复 `frame.energy`→`frame.scalars.energy`)，否则 benchmarks 整体不编译。

## Files to create or modify
- `benchmarks/bench_async_snapshot_mem.cu` (new) — 内存扫描基准 + 线性拟合 PASS/FAIL。
- `benchmarks/CMakeLists.txt` — 在 `ATV_BUILD_CUDA` 块注册 `bench_async_snapshot_mem`。
- `benchmarks/bench_scf.cu`, `benchmarks/bench_rhf.cu` — 前置修复 `frame.energy`→`frame.scalars.energy`（否则 benchmarks 不编译，挡住本环）。
- `.claude/NOTES.md` — freeze-then-dispatch 不变量、link-02 终态、device-ring 的 CUDA 执行模型分类条目。

## Tasks
- [x] Fix pre-existing benchmark break frame.energy -> frame.scalars.energy in benchmarks/bench_scf.cu + bench_rhf.cu (prerequisite: benchmarks must compile to build the membench)
- [x] Implement bench_async_snapshot_mem in benchmarks/bench_async_snapshot_mem.cu (K×N sweep allocating K device cudaMalloc + K pinned cudaMallocHost slots; device via cudaMemGetInfo, host via /proc/self/status VmRSS; emit K,N,bytes_per_elem,dev_bytes,host_bytes,model_bytes CSV)
- [x] Implement linear-fit residual check in benchmarks/bench_async_snapshot_mem.cu (assert device mem linear-in-K and linear-in-N within 10%, device/host ratio ~1 reported; non-zero exit on violation; print production kRingDepth=3 footprint row)
- [x] Add Doxygen header + units comments to public entry points in benchmarks/bench_async_snapshot_mem.cu (bytes, K, N, B)
- [x] Register bench_async_snapshot_mem in benchmarks/CMakeLists.txt under ATV_BUILD_CUDA via atv_add_bench
- [x] Update .claude/NOTES.md with freeze-then-dispatch invariant, link-02 end state (AsyncBridge device-ring freeze + AsyncRecordBase drop_oldest ring), and device-ring CUDA execution-model classification entry
- [x] Run full check + build (ATV_BUILD_CUDA=ON, ATV_BUILD_BENCHMARKS=ON) on GH200 and run bench_async_snapshot_mem (operator-action) — verified on GH200 (job 9627: build green, membench PASS residual 0.0000)

## Testing strategy
- Happy path:`bench_async_snapshot_mem` 在 GH200 上以 `-DATV_BUILD_CUDA=ON -DATV_BUILD_BENCHMARKS=ON` 构建成功并运行，发出每个 (K,N) 的 device+host 常驻内存数据点与 drop-vs-K 曲线。
- 边界情形:K=1（最小环深，drop 最多）与单点 N 不触发拟合除零;扫描含 ≥3 个 K 值与 ≥3 个 N 值以使最小二乘拟合良态;worker 节流回归下 `n_dropped()` 随 K 增大单调不增。
- 领域验证(标度律):分别固定 N 扫 K、固定 K 扫 N,线性拟合残差 ≤ 10%;CUDA 路径 device+host 合计 ≈ 2·K·N·S_union(双计比值在 build 时 S_union 推导值附近);残差越界则基准非零退出(模型被证伪)。
- 文档:`.claude/NOTES.md` 含 freeze-then-dispatch 不变量、统一终态条目与 device-ring 的 CUDA 执行模型分类行。

## Out of scope
- link 02 的 device_ring[K] / 共享 host ring 实现本身(本环依赖其已落地,仅做测量)。
- 把该基准接入 CI `ctest`(operator-action,登录节点无 nvcc,不作为可移植门)。
- 任何对 ring 容量/丢弃策略的优化或调参——本环只测量与建模,不改运行时行为。
- CPU-only 内存路径的独立基准(本环聚焦 CUDA 双计路径)。
