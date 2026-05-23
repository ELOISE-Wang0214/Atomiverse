# tests/cuda/async/CMakeLists.txt — `msu-pe-cuda-components-07` test suite.
#
# Registers the three guard tests for `cuda::AsyncRecordBase`:
#   * test_hot_path_no_sync  — 0 cuda::stream_synchronize on the hot path.
#   * test_async_record_correctness — async vs sync recorder field parity.
#   * test_async_record_perf — 25× recorder frequency wall-clock < 5%.
#
# All three carry the `gpu_required` ctest label so CI can skip them on
# GPU-less runners; they also carry `cuda;async` for partial selection.

atv_add_test(test_hot_path_no_sync         test_hot_path_no_sync.cu)
atv_add_test(test_async_record_correctness test_async_record_correctness.cu)
atv_add_test(test_async_record_perf        test_async_record_perf.cu)

set_tests_properties(
    test_hot_path_no_sync
    test_async_record_correctness
    test_async_record_perf
    PROPERTIES LABELS "cuda;async;gpu_required"
)

# Correctness + perf tests shell out to the polyethylene_tension benchmark
# and contend on the single GPU; serialize them so parallel ctest runs
# don't oversubscribe the device or race on the shared /tmp output dirs.
set_tests_properties(
    test_async_record_correctness
    test_async_record_perf
    PROPERTIES RUN_SERIAL TRUE
)

# Correctness + perf tests shell out to the polyethylene_tension benchmark.
# Surface its absolute path so the tests do not depend on cwd.
target_compile_definitions(test_async_record_correctness PRIVATE
    ATV_PE_TENSION_BIN="$<TARGET_FILE:polyethylene_tension>")
target_compile_definitions(test_async_record_perf PRIVATE
    ATV_PE_TENSION_BIN="$<TARGET_FILE:polyethylene_tension>")

# The async test files cross-link to the benchmark target's output path
# via add_dependencies so the binary exists before the tests run.
add_dependencies(test_async_record_correctness polyethylene_tension)
add_dependencies(test_async_record_perf        polyethylene_tension)
