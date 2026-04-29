include(FetchContent)

FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.3
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googlebenchmark)

function(atv_add_bench name)
    add_executable(${name} ${name}.cu)
    target_link_libraries(${name} PRIVATE atomiverse benchmark::benchmark)
    set_target_properties(${name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endfunction()

# C++ benchmark: compiled as .cpp, no CUDA required.
function(atv_add_cpp_bench name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE atomiverse benchmark::benchmark)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17)
endfunction()

if(ATV_BUILD_CUDA)
    atv_add_bench(bench_md_step)
    atv_add_bench(bench_scf)
    atv_add_bench(bench_rhf)
endif()

if(ATV_BUILD_RHF)
    atv_add_cpp_bench(bench_libcint)
    atv_add_cpp_bench(bench_native)

    # Phase 1.2 QC suite — doesn't link google-benchmark (plain main())
    add_executable(bench_qc_suite qc/bench_qc_suite.cpp)
    target_link_libraries(bench_qc_suite PRIVATE atomiverse)
    set_target_properties(bench_qc_suite PROPERTIES CXX_STANDARD 17)
    # Absolute path to benchmarks/qc so the binary can locate the
    # checked-in xyz geometries and JSON references regardless of cwd.
    target_compile_definitions(bench_qc_suite PRIVATE
        ATV_QC_DATA_DIR="${CMAKE_SOURCE_DIR}/benchmarks/qc")
endif()

if(ATV_BUILD_TORCH)
    atv_add_bench(bench_torch)
endif()

# ── CPU MD parallel-scaling benchmarks (run under mpirun) ──────────
# These are plain executables (not google-benchmark) so they can be
# launched at varying rank counts and emit CSV one-liners that an
# external driver (tools/run_msu_*.sh, scaling-report scripts) collects.
function(atv_add_cpp_md_bench name)
    add_executable(${name} md/${name}.cpp)
    target_link_libraries(${name} PRIVATE atomiverse)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17)
    if(ATV_USE_MPI AND TARGET MPI::MPI_CXX)
        target_link_libraries(${name} PRIVATE MPI::MPI_CXX)
    endif()
endfunction()

atv_add_cpp_md_bench(scaling)
atv_add_cpp_md_bench(polyethylene_tension)
