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
atv_add_plain_bench(scaling md/lj/performance/scaling.cpp)
if(ATV_USE_MPI AND TARGET MPI::MPI_CXX)
    target_link_libraries(scaling PRIVATE MPI::MPI_CXX)
endif()

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

# ════════════════════════════════════════════════════════════════════
# pw — plane-wave DFT (bulk Si vs CP2K)
# ════════════════════════════════════════════════════════════════════

if(ATV_BUILD_PW)
    # validation: Si fcc equation of state (a0, B0) vs CP2K reference.
    atv_add_plain_bench(bench_pw_si_cp2k_eos
                        pw/validation/bench_pw_si_cp2k_eos.cpp)
    # performance: SCF time-to-solution CSV (cpu/cuda rows; the cp2k
    # rows come from pw/performance/cp2k/run_perf.sh).
    atv_add_plain_bench(bench_pw_scf_perf
                        pw/performance/bench_pw_scf_perf.cpp)
    foreach(name bench_pw_si_cp2k_eos bench_pw_scf_perf)
        # Absolute path so the binaries locate the checked-in UPF and
        # CP2K reference CSV regardless of cwd.
        target_compile_definitions(${name} PRIVATE
            ATV_PW_DATA_DIR="${CMAKE_SOURCE_DIR}/benchmarks/pw/validation/reference")
        # The #ifdef-gated `--backend cuda` path reaches cuda::pw::pw_scf
        # compiled into the atomiverse library.
        if(ATV_BUILD_CUDA)
            target_compile_definitions(${name} PRIVATE ATV_BUILD_CUDA)
        endif()
    endforeach()
endif()
