# ── Benchmark layout ────────────────────────────────────────────────
# Every benchmark system has exactly two parts:
#   <system>/validation/   scientific correctness vs an external reference
#                          (LAMMPS, PySCF, …) — plain executables, non-zero
#                          exit on any failing case; run manually, not by CI.
#   <system>/performance/  speed / scaling — plain executables or
#                          google-benchmark, emitting CSV for the plotters.
# Systems: md/lj (argon Lennard-Jones), md/pe (polyethylene Dreiding-UA),
# qc (Gaussian-basis quantum chemistry). Shared MD infra (PASS/FAIL
# harness, sbatch scan drivers, plotter) lives in md/common/.

include(FetchContent)

FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.3
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googlebenchmark)

# ── helpers ─────────────────────────────────────────────────────────

# CUDA google-benchmark executable (source path relative to this dir).
function(atv_add_cuda_bench name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE atomiverse benchmark::benchmark)
    set_target_properties(${name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endfunction()

# C++ google-benchmark executable.
function(atv_add_cpp_bench name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE atomiverse benchmark::benchmark)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17)
endfunction()

# Plain C++ executable (no google-benchmark): CSV-emitting performance
# drivers and PASS/FAIL validation benches.
function(atv_add_plain_bench name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE atomiverse)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17)
endfunction()

# ════════════════════════════════════════════════════════════════════
# md/lj — argon Lennard-Jones
# ════════════════════════════════════════════════════════════════════

# validation: Atomiverse MD components vs LAMMPS reference behaviour and
# reference ensemble averages (bench_check.h PASS/FAIL harness).
foreach(name
    bench_setup_lammps_align
    bench_deform_lammps_align
    bench_nve_limit_lammps_align
    bench_thermostat_lammps_align
    bench_pressure_lammps_align
    bench_lammps_langevin_parity
    bench_lammps_npt_parity)
    atv_add_plain_bench(${name} md/lj/validation/${name}.cpp)
    target_include_directories(${name} PRIVATE md/common)
endforeach()

# performance: CPU MPI weak/strong scaling (run by
# md/common/bench_scaling_cpu.sbatch at 1..N ranks; emits CSV).
# `scaling` = MD (NVE), `minimize_scaling` = minimizer (CG/LBFGS, fixed
# iteration budget); both share the sc-LJ system in sc_lj_system.h.
foreach(name scaling minimize_scaling)
    atv_add_plain_bench(${name} md/lj/performance/${name}.cpp)
    if(ATV_USE_MPI AND TARGET MPI::MPI_CXX)
        target_link_libraries(${name} PRIVATE MPI::MPI_CXX)
    endif()
endforeach()

if(ATV_BUILD_CUDA)
    # performance: single-GPU MD step microbench (google-benchmark).
    atv_add_cuda_bench(bench_md_step md/lj/performance/bench_md_step.cu)
    # performance: single-process multi-GPU strong scaling (plain main(),
    # argv-driven GPU count / system size; run by
    # md/common/bench_scaling.sbatch).
    add_executable(bench_multi_gpu_scaling
                   md/lj/performance/bench_multi_gpu_scaling.cu)
    target_link_libraries(bench_multi_gpu_scaling PRIVATE atomiverse)
    set_target_properties(bench_multi_gpu_scaling PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
    # performance: single-GPU CUDA minimizer (CG / LBFGS) vs LAMMPS
    # `minimize` (deck in.lj_min_kk); also emits the LAMMPS data file so
    # both engines relax the IDENTICAL configuration.
    add_executable(bench_minimize md/lj/performance/bench_minimize.cu)
    target_link_libraries(bench_minimize PRIVATE atomiverse)
    set_target_properties(bench_minimize PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
endif()

if(ATV_BUILD_TORCH)
    # performance: LibTorch backend LJ forward/backward + transfer.
    atv_add_cuda_bench(bench_torch md/lj/performance/bench_torch.cu)
endif()

# ════════════════════════════════════════════════════════════════════
# md/pe — polyethylene melt (Dreiding-UA)
# ════════════════════════════════════════════════════════════════════

# validation: C20H42 bonded-energy cross-check vs LAMMPS (1e-6 eV).
atv_add_plain_bench(bench_pe_bonded_lammps_xcheck
                    md/pe/validation/bench_pe_bonded_lammps_xcheck.cpp)
target_include_directories(bench_pe_bonded_lammps_xcheck PRIVATE md/common)

# validation: timer-instrumentation capture gate (per-category rollup >= 95%
# of step-loop wall). A wall-clock-fraction PERF assertion — lives here, not
# in the unit-test suite (UTs assert code correctness only).
atv_add_plain_bench(bench_timing_capture
                    md/pe/validation/bench_timing_capture.cpp)

# performance: PE melt MSU-protocol stages (CPU MPI + optional CUDA
# backend; run by the md/common scan drivers; emits stage timing CSV).
atv_add_plain_bench(polyethylene_tension
                    md/pe/performance/polyethylene_tension.cpp)
if(ATV_USE_MPI AND TARGET MPI::MPI_CXX)
    target_link_libraries(polyethylene_tension PRIVATE MPI::MPI_CXX)
endif()
# The #ifdef-gated `--backend cuda` path reaches the cuda::md::*
# components compiled into the atomiverse library.
if(ATV_BUILD_CUDA)
    target_compile_definitions(polyethylene_tension PRIVATE ATV_BUILD_CUDA)

    # performance: 25× recorder-frequency wall-clock overhead gate (< 5%).
    # Shells out to polyethylene_tension --backend cuda; non-zero exit on
    # violation. Moved from the unit-test suite (perf gates live here).
    atv_add_plain_bench(bench_async_record_overhead
                        md/pe/performance/bench_async_record_overhead.cpp)
    target_compile_definitions(bench_async_record_overhead PRIVATE
        ATV_PE_TENSION_BIN="$<TARGET_FILE:polyethylene_tension>")
    add_dependencies(bench_async_record_overhead polyethylene_tension)
endif()

# ════════════════════════════════════════════════════════════════════
# qc — Gaussian-basis quantum chemistry
# ════════════════════════════════════════════════════════════════════

if(ATV_BUILD_RHF)
    # validation: Minimize + SCFEvaluator vs checked-in reference JSONs.
    add_executable(bench_qc_suite qc/validation/bench_qc_suite.cpp)
    target_link_libraries(bench_qc_suite PRIVATE atomiverse)
    set_target_properties(bench_qc_suite PROPERTIES CXX_STANDARD 17)
    # Absolute path so the binary locates the checked-in xyz geometries
    # and JSON references regardless of cwd.
    target_compile_definitions(bench_qc_suite PRIVATE
        ATV_QC_DATA_DIR="${CMAKE_SOURCE_DIR}/benchmarks/qc/validation")

    # performance: native-kernel vs libcint ERI/integral speed.
    atv_add_cpp_bench(bench_native qc/performance/bench_native.cpp)
    atv_add_cpp_bench(bench_libcint qc/performance/bench_libcint.cpp)
endif()
