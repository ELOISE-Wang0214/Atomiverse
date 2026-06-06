---
slug: async-analysis-snapshot-ring-03-membench
criteria:
  - id: ac-001
    summary: bench_async_snapshot_mem builds and is registered under the CUDA gate
    type: code
    pass_when: |
      benchmarks/CMakeLists.txt registers bench_async_snapshot_mem via
      atv_add_bench inside the if(ATV_BUILD_CUDA) block, and
      benchmarks/bench_async_snapshot_mem.cu exists with a main() entry point.
    status: verified
    last_checked: 2026-06-06
  - id: ac-002
    summary: Benchmark runs on GH200 and emits per-(K,N) device+host memory CSV
    type: performance
    evaluator_hint: operator-action on aarch64 GH200 node (ATV_BUILD_CUDA=ON + ATV_BUILD_BENCHMARKS=ON)
    pass_when: |
      On the GH200 node (-DATV_BUILD_CUDA=ON -DATV_BUILD_BENCHMARKS=ON), the
      built bench_async_snapshot_mem emits one CSV row per (K,N) with columns
      K,N,bytes_per_elem,dev_bytes,host_bytes,model_bytes over >=3 K values and
      >=3 N values, and prints the production kRingDepth=3 footprint row.
    status: verified
    last_checked: 2026-06-06
  - id: ac-003
    summary: Measured device memory is linear in K and in N within 10%
    type: scientific
    pass_when: |
      Least-squares linear fit of measured device resident memory (cudaMemGetInfo
      delta) vs K (fixed N) and vs N (fixed K) has residual <= 10% for the
      large-N points (where allocation granularity is negligible); the benchmark
      exits non-zero if either residual exceeds tolerance.
    status: verified
    last_checked: 2026-06-06
  - id: ac-004
    summary: Device and host rings mirror in size; combined footprint ~2*K*N*B
    type: scientific
    pass_when: |
      The benchmark reports the device-ring (cudaMalloc) and pinned host-ring
      (cudaMallocHost) resident sizes per (K,N); the device/host ratio is ~1
      (each ~K*N*B), so the combined footprint is ~2*K*N*B, matching the model;
      the production kRingDepth=3 footprint is reported.
    status: verified
    last_checked: 2026-06-06
  - id: ac-005
    summary: Public benchmark entry points carry Doxygen with units
    type: docs
    pass_when: |
      Public functions in benchmarks/bench_async_snapshot_mem.cu have @brief and
      @param/@return documenting units (bytes for memory, dimensionless K/N, B in
      bytes-per-element).
    status: pending
  - id: ac-006
    summary: NOTES.md records the invariant, link-02 end state, device-ring classification
    type: docs
    pass_when: |
      .claude/NOTES.md contains the freeze-then-dispatch invariant, the link-02
      end state (AsyncBridge compute-stream device-ring freeze across the 5
      subclasses + AsyncRecordBase internal device_ring[K] with drop_oldest), and
      a device-ring entry in the CUDA execution-model classification table.
    status: pending
---

# Acceptance criteria

- **ac-001** — source exists + registered under `ATV_BUILD_CUDA` via `atv_add_bench`; static check (code).
- **ac-002** — operator builds + runs on GH200, emits per-(K,N) device+host memory CSV + production footprint (performance, needs real device measurement).
- **ac-003** — the falsifiable scaling law: device memory linear in K and N within 10% (large-N points); non-zero exit on violation.
- **ac-004** — device ring ≈ host ring (mirror), so combined ≈ 2·K·N·B; production K=3 footprint reported. (Reframed from the original drop-vs-K, which is not a runtime sweep — K is the compile-time `kRingDepth`.)
- **ac-005** — Doxygen + units.
- **ac-006** — NOTES: invariant + link-02 end state + device-ring classification row.
