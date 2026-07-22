# tests/cuda/async/CMakeLists.txt — AsyncRecordBase / AsyncBridge guards.
#
# Unit-only guards (always registered under ctest):
#   * test_hot_path_no_sync      — 0 cuda::stream_synchronize on the hot path.
#   * test_async_bridge_tearing  — AsyncBridge device-ring freeze delivers only
#                                  single-step, correctly-stamped snapshots
#                                  (async-analysis-snapshot-ring-02-device-ring).
#   * test_async_rhf_headers     — header-only async recorder compile surface.
#
# E2E deliberately NOT registered here (tester contract: no e2e under tests/):
#   * test_async_record_correctness.cu shells out via system() to the
#     polyethylene_tension binary — full MD + CSV parity. Keep the source as a
#     reference harness; re-home under benchmarks/ (or regressions/) when that
#     path is rebuilt. The 25× wall-clock gate already lives in
#     benchmarks/md/pe/performance/bench_async_record_overhead.cpp.
#
# All carry the `gpu_required` ctest label so CI can skip them on GPU-less
# runners; they also carry `cuda;async` for partial selection.

atv_add_test(test_hot_path_no_sync        test_hot_path_no_sync.cu)
atv_add_test(test_async_bridge_tearing    test_async_bridge_tearing.cu)
atv_add_test(test_async_rhf_headers       test_async_rhf_headers.cu)

set_tests_properties(
    test_hot_path_no_sync
    test_async_bridge_tearing
    test_async_rhf_headers
    PROPERTIES LABELS "cuda;async;gpu_required"
)
# The tearing test races a lagging host worker against the MD loop on the
# single GPU; serialize it so parallel ctest runs don't oversubscribe.
set_tests_properties(test_async_bridge_tearing PROPERTIES RUN_SERIAL TRUE)
