---
title: Async snapshot ring — memory benchmark & scaling-law verification
status: approved
created: 2026-06-06
---

# Async snapshot ring — memory benchmark & scaling-law verification

## Summary
为异步快照环（async snapshot ring）增加一个常驻内存基准与其标度律的可证伪验证，作为 `async-analysis-snapshot-ring` 链条的第 3 环（收尾环）。基准在 GH200 节点上扫描环深 K、原子数 N 与捕获字段并集大小 S_union，分别用 `cudaMemGetInfo`（设备）与 `/proc/self/status` 的 VmHWM（主机）测量常驻内存，发出每个 (K,N) 的内存数据点，并在受控 worker 吞吐下用 `AnalysisExecutor::n_dropped()` 给出 drop-vs-K 的取舍曲线。基准对线性拟合残差设阈值，从而证伪或确认内存模型 mem ≈ K·N·S_union（CUDA 路径计两次）。同时把 freeze-then-dispatch 不变量、统一终态与设备环路径的 CUDA 执行模型分类落档到 `.claude/NOTES.md`。

## Domain basis
本环不引入新物理（science.required=true，特此声明：无新方程、无 DOI/arXiv 引用）。其"领域基础"是一个可证伪的线性内存标度模型：

- 常驻快照内存 mem ≈ K · N · S_union
  - K = 全局环深（`AnalysisExecutor::Config.ring_size`）
  - N = 原子数
  - S_union = 捕获字段并集的每原子字节数（build 时一次性计算的所有订阅者 CaptureMask 的静态并集）
- CUDA 路径同时持有 device_ring[K] 与共享 host ring，故该内存被计两次：mem_cuda ≈ 2 · K · N · S_union。
- 可证伪形式：固定其余变量分别扫描 K 与 N，对测得内存做线性拟合，残差须 ≤ 10%；否则模型被证伪。device/host 两路的比值须在 build 时声明的 S_union 比例附近（CUDA 双计为约 2×）。

## Design
- 新增基准可执行体 `bench_async_snapshot_mem`（plain `main()`，argv 驱动，发出 CSV 单行，与 `bench_multi_gpu_scaling` 同风格——不链接 google-benchmark），位于 `benchmarks/bench_async_snapshot_mem.cu`。
- 测量入口：构造 `cpu::AnalysisExecutor`（经 link 02 后带 device_ring[K]），用真实订阅者集合驱动 build，使 S_union 由 CaptureMask 并集确定;对一组 (K,N) 笛卡尔扫描，每点：
  - 设备常驻：`cudaMemGetInfo(&free,&total)` 在分配前后取 free 差值。
  - 主机常驻：读 `/proc/self/status` 的 `VmHWM`（高水位）。
  - drop 计数：在固定的"慢 worker"节流回归下推送已知数量的快照，读 `AnalysisExecutor::n_dropped()`。
- 输出：每个 (K,N) 一行 `K,N,S_union,dev_bytes,host_bytes,n_dropped`；并发出 drop-vs-K 曲线（固定 N，扫 K）。
- 标度律检查：基准内做线性拟合（最小二乘）并打印 K 维、N 维残差与 CUDA/host 双计比值；残差越界时以非零退出码失败，使其可作为 operator 运行的 PASS/FAIL 门。
- CMake：在 `benchmarks/CMakeLists.txt` 的 `if(ATV_BUILD_CUDA)` 块内用 `atv_add_bench(bench_async_snapshot_mem)` 注册（受 `ATV_BUILD_BENCHMARKS` + `ATV_BUILD_CUDA` 双重门控）。
- 文档：`.claude/NOTES.md` 追加 freeze-then-dispatch 不变量、统一终态（union-freeze-into-ring + drop_oldest + 全局 K + 共享 worker 池），并在 CUDA 执行模型分类表中新增 device-ring 路径条目。

> 注意:本基准为 **operator-action**——仅能在 aarch64 GH200 节点上以 `ATV_BUILD_CUDA=ON` + `ATV_BUILD_BENCHMARKS=ON` 构建运行（登录节点无 nvcc，见 CLAUDE.md GPU-node 指南）。

## Files to create or modify
- `benchmarks/bench_async_snapshot_mem.cu` (new) — 内存扫描基准 + 线性拟合 PASS/FAIL。
- `benchmarks/CMakeLists.txt` — 在 `ATV_BUILD_CUDA` 块注册 `bench_async_snapshot_mem`。
- `.claude/NOTES.md` — freeze-then-dispatch 不变量、统一终态、device-ring 的 CUDA 执行模型分类条目。

## Tasks
- [ ] Implement bench_async_snapshot_mem in benchmarks/bench_async_snapshot_mem.cu (K×N sweep; device via cudaMemGetInfo, host via /proc/self/status VmHWM; emit K,N,S_union,dev_bytes,host_bytes,n_dropped CSV + drop-vs-K curve)
- [ ] Implement linear-fit residual check in benchmarks/bench_async_snapshot_mem.cu (assert mem linear-in-K and linear-in-N within 10%, CUDA dev/host double-count ratio reported; non-zero exit on violation)
- [ ] Add Doxygen header + units comments to public entry points in benchmarks/bench_async_snapshot_mem.cu (bytes, K, N, S_union)
- [ ] Register bench_async_snapshot_mem in benchmarks/CMakeLists.txt under ATV_BUILD_CUDA via atv_add_bench
- [ ] Update .claude/NOTES.md with freeze-then-dispatch invariant, unified end state, and device-ring CUDA execution-model classification entry
- [ ] Run full check + build (ATV_BUILD_CUDA=ON, ATV_BUILD_BENCHMARKS=ON) on GH200 and run bench_async_snapshot_mem (operator-action)

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
