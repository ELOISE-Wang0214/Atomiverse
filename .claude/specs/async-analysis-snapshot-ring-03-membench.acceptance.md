---
slug: async-analysis-snapshot-ring-03-membench
criteria:
  - id: ac-001
    summary: bench_async_snapshot_mem builds and is registered under CUDA gate
    type: code
    pass_when: |
      benchmarks/CMakeLists.txt registers bench_async_snapshot_mem via
      atv_add_bench inside the if(ATV_BUILD_CUDA) block, and
      benchmarks/bench_async_snapshot_mem.cu exists with a main() entry point.
    status: pending
  - id: ac-002
    summary: Benchmark runs on GH200 and emits per-(K,N) memory + drop-vs-K data
    type: performance
    evaluator_hint: operator-action on aarch64 GH200 node (ATV_BUILD_CUDA=ON + ATV_BUILD_BENCHMARKS=ON)
    pass_when: |
      On the GH200 node, configured with -DATV_BUILD_CUDA=ON
      -DATV_BUILD_BENCHMARKS=ON, the built bench_async_snapshot_mem emits one
      CSV row per (K,N) with columns K,N,S_union,dev_bytes,host_bytes,n_dropped
      and a drop-vs-K curve over >=3 K values.
    status: pending
  - id: ac-003
    summary: Measured memory linear in K and N within 10%, CUDA double-counted
    type: scientific
    pass_when: |
      Least-squares linear fit of measured resident memory vs K (fixed N) and
      vs N (fixed K) has residual <= 10%; combined CUDA device+host total
      approximates 2*K*N*S_union (dev/host double-count ratio near the
      build-time S_union value); benchmark exits non-zero if any residual
      exceeds tolerance.
    status: pending
  - id: ac-004
    summary: Drop count is monotonic non-increasing in K under fixed worker throttle
    type: scientific
    pass_when: |
      Under the controlled slow-worker regime, n_dropped() reported per K is
      monotonically non-increasing as K grows across the swept K values.
    status: pending
  - id: ac-005
    summary: Public benchmark entry points carry Doxygen with units
    type: docs
    pass_when: |
      Public functions in benchmarks/bench_async_snapshot_mem.cu have @brief
      and @param/@return documenting units (bytes for memory, dimensionless
      K/N/S_union).
    status: pending
  - id: ac-006
    summary: NOTES.md records invariant, unified end state, device-ring classification
    type: docs
    pass_when: |
      .claude/NOTES.md contains the freeze-then-dispatch invariant, the unified
      end state (union-freeze-into-ring + drop_oldest + global K + shared worker
      pool), and a device-ring entry in the CUDA execution-model classification.
    status: pending
---

# Acceptance criteria

- **ac-001** — 基准源文件存在且在 `ATV_BUILD_CUDA` 门下经 `atv_add_bench` 注册;纯静态检查(code)。
- **ac-002** — operator 在 GH200 上构建并运行,发出 (K,N) 内存点与 drop-vs-K 曲线;performance 类型,因须真实设备测量。
- **ac-003** — 核心可证伪标度律:K、N 线性拟合残差 ≤ 10%,CUDA 双计合计 ≈ 2·K·N·S_union,越界则非零退出。
- **ac-004** — drop 行为正确性:固定节流下 `n_dropped()` 随 K 单调不增。
- **ac-005** — 公共入口的 Doxygen + 单位注释。
- **ac-006** — NOTES 落档:不变量 + 统一终态 + device-ring 分类条目。
