---
slug: msu-pe-cuda-components-07-async-correctness
criteria:
  - id: ac-001
    summary: cuda::AsyncRecordBase template header exists with submit_async + lifecycle
    type: code
    pass_when: |
      `include/atomiverse/cuda/async/record_base.cuh` declares
      `template <typename Snapshot> class cuda::AsyncRecordBase :
      public cuda::ComponentBase` with public members
      `void submit_async(int64_t step, std::function<void(const Snapshot&,
      int64_t)> work_fn)`, `void build(BuildContext&)`, `void destroy(BuildContext&)`,
      protected fields `Snapshot* device_snapshot_`, `Snapshot host_snapshot_`,
      `cudaEvent_t copy_done_event_`, `cudaStream_t defer_stream_`.
      Compiles standalone via `cmake --build build --target atomiverse_cuda`.
    status: verified
    last_checked: 2026-05-23
  - id: ac-002
    summary: cuda::stream_synchronize wrapper and ATV_CUDA_TRACK_SYNCS counter declared
    type: code
    pass_when: |
      `include/atomiverse/cuda/async/sync_guard.cuh` declares a function
      `cuda::stream_synchronize(cudaStream_t)` that increments an atomic
      counter when `ATV_BUILD_TESTS` is defined and the counter has read/reset
      accessors `cuda::detail::sync_count()` and `cuda::detail::reset_sync_count()`.
      Release builds (without `ATV_BUILD_TESTS`) compile the wrapper down to a
      direct `cudaStreamSynchronize` call with zero counter overhead (verified
      by reading the preprocessed source).
    status: verified
    last_checked: 2026-05-23
  - id: ac-003
    summary: CudaCsvRecorder inherits AsyncRecordBase; ad-hoc async code removed
    type: code
    pass_when: |
      In `benchmarks/md/polyethylene_tension.cpp`, the `CudaCsvRecorder`
      definition includes `: public cuda::AsyncRecordBase<SystemScalarsSnapshot>`
      and no longer contains direct calls to `cudaMemcpyAsync` or
      `cudaEventRecord` — all D2H staging goes through
      `AsyncRecordBase::submit_async`. `git grep cudaMemcpyAsync
      benchmarks/md/polyethylene_tension.cpp` returns no matches.
    status: verified
    last_checked: 2026-05-23
  - id: ac-004
    summary: hot-path sync count remains zero across 100 cuda::md::MD steps
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      `ctest --test-dir build -R test_hot_path_no_sync --output-on-failure -L cuda`
      passes on a CUDA-visible host. The test runs
      `cuda::md::MD::step()` 100 times with a bound `CudaCsvRecorder`
      (`print_every=1`); after every step,
      `cuda::detail::sync_count() == 0` (gtest EXPECT_EQ). If
      `ATV_BUILD_TESTS=OFF`, the test is not registered and ctest does not
      list it.
    status: verified
    last_checked: 2026-05-23
  - id: ac-005
    summary: async vs sync recorder produce field-equivalent CSVs
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      `ctest --test-dir build -R test_async_record_correctness --output-on-failure
      -L cuda` passes. The test runs two `polyethylene_tension --backend cuda
      --mode stage:2 --scale 0.001` runs with the same seed (one async via
      `AsyncRecordBase`, one with `-DATV_TEST_RECORD_SYNC=1` forcing
      `cudaMemcpy` sync path). For every row index present in BOTH CSVs, every
      numeric column satisfies `|x_async - x_sync| <= 1e-6 * max(|x_sync|, 1.0)`.
      The row-count difference between the two CSVs is `<= 1`.
    status: verified
    last_checked: 2026-05-23
  - id: ac-006
    summary: 25x recorder-frequency overhead stays below 5% wall-clock
    type: performance
    evaluator_hint: gpu_required
    pass_when: |
      `ctest --test-dir build -R test_async_record_perf --output-on-failure
      -L cuda` passes. The test runs
      `polyethylene_tension --backend cuda --mode stage:2 --scale 0.005
      --print-every 25` (control) and `--print-every 1` (high-frequency)
      five times each, takes the median wall-clock per configuration, and
      asserts `wall_high / wall_low - 1 < 0.05`.
    status: verified
    last_checked: 2026-05-23
  - id: ac-007
    summary: lifecycle build/destroy round-trip leaks no device memory
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      Running the three new tests under
      `compute-sanitizer --tool memcheck` reports zero leaks attributable to
      `device_snapshot_` and `copy_done_event_` after three back-to-back
      `build → launch → destroy` cycles. After the final destroy,
      `device_snapshot_ == nullptr` and `copy_done_event_ == nullptr` (gtest
      EXPECT_EQ).
    status: pending
    last_checked: 2026-05-23
    notes: |
      The `LifecycleRoundTrip` test in tests/cuda/async/test_hot_path_no_sync.cu
      runs the 3× build→launch→destroy cycle and passes; AsyncRecordBase::destroy
      manifestly cudaFrees device_snapshot_ and cudaEventDestroys both events
      (see include/atomiverse/cuda/async/record_base.cuh). `compute-sanitizer` is
      not installed in ~/atv-build-env/cuda so the formal memcheck gate is
      deferred to whoever has the tool; the code-review path covers the
      logical lifecycle invariant.
  - id: ac-008
    summary: doxygen documents the zero-block contract, defer-stream, and worker semantics
    type: docs
    pass_when: |
      `include/atomiverse/cuda/async/record_base.cuh` carries a class-level
      doxygen block that states (a) `submit_async` is non-blocking on the
      caller's thread, (b) the D2H copy uses
      `defer_stream` = `cuda::md::MD::async_stream()` not the compute stream,
      (c) `cudaEventSynchronize` runs on the `AnalysisExecutor` worker
      thread, not the MD hot-path. Every public member has `@brief` and any
      pointer/handle parameter is documented as owning/non-owning.
    status: verified
    last_checked: 2026-05-23
  - id: ac-009
    summary: spec index lists the new chain link
    type: docs
    pass_when: |
      `.claude/specs/INDEX.md` contains a line referencing
      `msu-pe-cuda-components-07-async-correctness` with a one-line
      summary and an `[approved]` (or later) lifecycle tag.
    status: verified
    last_checked: 2026-05-23
  - id: ac-010
    summary: full check + test suite green after the change
    type: runtime
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_BUILD_CUDA=ON
      -DATV_BUILD_BENCHMARKS=ON && cmake --build build -j && ctest --test-dir
      build --output-on-failure` exits 0 (CUDA-gated tests may be `Skipped`
      on non-GPU hosts; no test reports `Failed`).
    status: verified
    last_checked: 2026-05-23
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003** are static (code-shape) gates: the base class exists, the sync-guard wrapper exists with correct release-build cost, and the chain-link-06 recorder is refactored onto the base class with no leftover ad-hoc async code.
- **ac-004** is the headline correctness gate of this link: hot-path `cudaStreamSynchronize` count must be zero. This is the assertion that protects every future CUDA recorder added to the project.
- **ac-005** verifies field-level equivalence between the async and sync code paths on the same seed — proves the async path doesn't drop data or corrupt fields (only possibly drops the trailing row on early exit).
- **ac-006** is the headline performance gate: a 25× increase in recorder frequency must cost < 5% wall-clock. Failure indicates a missed sync somewhere on the hot path.
- **ac-007** ensures lifecycle cleanliness — the only resources that grow with build/destroy cycles are device snapshot buffer and the cudaEvent, both must be freed.
- **ac-008 / ac-009** are documentation gates. ac-010 is the full-suite green gate the chain ends on.

Cross-engine numerical comparison (async vs LAMMPS) is deliberately not an ac here — that lives in chain link 08.
