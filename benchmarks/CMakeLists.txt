# ── Benchmark layout (see README.md) ────────────────────────────────
# <domain>/<case>/{validation,performance} — every case has exactly two
# parts:
#   validation/   scientific correctness vs an external reference
#                 (LAMMPS, PySCF, CP2K, …) — plain executables, non-zero
#                 exit on any failing case; run manually, not by CI.
#   performance/  speed / scaling — plain executables or
#                 google-benchmark, emitting CSV for the plotters.
# Cases: md/lj (argon LJ), md/pe (polyethylene Dreiding-UA),
# qc/rhf (Gaussian-basis small-molecule suite), pw/si (bulk Si LDA).
# Shared per-domain infra (PASS/FAIL harness, sbatch scan drivers,
# plotter) lives in <domain>/common/.

include(FetchContent)

FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.3
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googlebenchmark)

# GoogleTest for the relocated end-to-end / parity validation binaries (moved
# out of the fast unit-test gate — see the fast-gate note in tests/CMakeLists.txt).
# Idempotent with the tests/ fetch when both are built.
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

# ── helpers ─────────────────────────────────────────────────────────

# A gtest-based validation binary: scientific correctness vs an external
# reference (PySCF, libcint, …), registered as a ctest test. These are the
# relocated end-to-end checks — they read as small, orthogonal, learnable
# examples of one capability each, run only under -DATV_BUILD_BENCHMARKS=ON,
# never in the fast unit gate. `src` is relative to this dir; the source file's
# own directory is on the include path so a co-located `reference/` header
# resolves. (validation/ convention is plain executables; these gtest binaries
# are the exception for the relocated suites.)
function(atv_add_validation name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE atomiverse GTest::gtest_main)
    get_filename_component(_dir "${CMAKE_CURRENT_SOURCE_DIR}/${src}" DIRECTORY)
    target_include_directories(${name} PRIVATE "${_dir}")
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

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

# backend-consistency: CPU (all-double reference) writes an energy/force golden;
# CUDA (mixed precision) replays the IDENTICAL config and asserts per-quantity
# tolerances sized to the float-storage floor (NOT bitwise). See
# md/lj/consistency/README.md. CUDA half is nightly/GPU-only.
atv_add_plain_bench(bench_lj_consistency_cpu md/lj/consistency/cpu_reference.cpp)
if(ATV_BUILD_CUDA)
    add_executable(bench_lj_consistency_cuda md/lj/consistency/cuda_compare.cu)
    target_link_libraries(bench_lj_consistency_cuda PRIVATE atomiverse)
    set_target_properties(bench_lj_consistency_cuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
    target_compile_definitions(bench_lj_consistency_cuda PRIVATE ATV_BUILD_CUDA)
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
    add_executable(bench_qc_suite qc/rhf/validation/bench_qc_suite.cpp)
    target_link_libraries(bench_qc_suite PRIVATE atomiverse)
    set_target_properties(bench_qc_suite PROPERTIES CXX_STANDARD 17)
    # Absolute path so the binary locates the checked-in xyz geometries
    # and JSON references regardless of cwd.
    target_compile_definitions(bench_qc_suite PRIVATE
        ATV_QC_DATA_DIR="${CMAKE_SOURCE_DIR}/benchmarks/qc/rhf/validation")

    # performance: native-kernel vs libcint ERI/integral speed.
    atv_add_cpp_bench(bench_native qc/rhf/performance/bench_native.cpp)
    atv_add_cpp_bench(bench_libcint qc/rhf/performance/bench_libcint.cpp)

    # convergence: RHF energy vs basis set (variational monotone decrease).
    atv_add_plain_bench(bench_qc_basis_convergence
                        qc/rhf/convergence/bench_basis_convergence.cpp)

    # validation: vibrational analysis (finite-difference Hessian) — relocated
    # from the unit-test gate (too slow: each is many displaced SCF solves).
    # Each is an orthogonal learnable check of one vib capability: the FD-Hessian
    # Response driver, mass-weighting, translation/rotation projection, harmonic
    # frequencies vs PySCF, and thermochemistry vs PySCF.
    atv_add_validation(val_vib_numerical_hessian  qc/vib/validation/test_numerical_hessian.cpp)
    atv_add_validation(val_vib_mass_weighting     qc/vib/validation/test_mass_weighting.cpp)
    atv_add_validation(val_vib_tr_projection      qc/vib/validation/test_tr_projection.cpp)
    atv_add_validation(val_vib_frequencies_vs_pyscf qc/vib/validation/test_frequencies_vs_pyscf.cpp)
    atv_add_validation(val_vib_thermochem_vs_pyscf  qc/vib/validation/test_thermochem_vs_pyscf.cpp)
    atv_add_validation(val_vib_output_components   qc/vib/validation/test_vib_output_components.cpp)

    # validation: Gaussian SCF end-to-end / parity — relocated from the unit gate
    # (each runs a full SCF solve / analytic-gradient FD sweep / f-shell ERI /
    # grid build to convergence: seconds to ~14 min). Orthogonal learnable checks
    # of one capability each: hybrid single-point, analytic gradients (RKS d-shell,
    # KS grid-response, hybrid, def2-SVP, unified), gradient-vs-libcint, XC grid,
    # f-shell ERI, d-shell basis eval, energy-vs-PySCF, AIMD energy conservation.
    foreach(v
        hybrid_singlepoint ks_dshell_gradient ks_grid_response_gradient
        hybrid_gradient def2svp_gradient gradient gradient_vs_libcint
        grid rys_eri_fshell basis_eval_dshell basis_vs_pyscf aimd_conservation)
        atv_add_validation(val_scf_${v} qc/rhf/validation/test_${v}.cpp)
    endforeach()
    # hybrid_singlepoint locates its molecule battery relative to the source tree.
    target_compile_definitions(val_scf_hybrid_singlepoint PRIVATE
        ATV_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
endif()

# ════════════════════════════════════════════════════════════════════
# pw — plane-wave DFT (bulk Si vs CP2K)
# ════════════════════════════════════════════════════════════════════

if(ATV_BUILD_PW)
    # validation: Si fcc equation of state (a0, B0) vs CP2K reference.
    atv_add_plain_bench(bench_pw_si_cp2k_eos
                        pw/si/validation/bench_pw_si_cp2k_eos.cpp)
    # performance: SCF time-to-solution CSV (cpu/cuda rows; the cp2k
    # rows come from pw/si/performance/cp2k/run_perf.sh).
    atv_add_plain_bench(bench_pw_scf_perf
                        pw/si/performance/bench_pw_scf_perf.cpp)
    foreach(name bench_pw_si_cp2k_eos bench_pw_scf_perf)
        # Absolute path so the binaries locate the checked-in UPF and
        # CP2K reference CSV regardless of cwd.
        target_compile_definitions(${name} PRIVATE
            ATV_PW_DATA_DIR="${CMAKE_SOURCE_DIR}/benchmarks/pw/si/validation/reference")
        # The #ifdef-gated `--backend cuda` path reaches cuda::pw::pw_scf
        # compiled into the atomiverse library.
        if(ATV_BUILD_CUDA)
            target_compile_definitions(${name} PRIVATE ATV_BUILD_CUDA)
        endif()
    endforeach()

    # convergence: Si energy vs ecut & k-mesh (variational / shrinking-delta).
    atv_add_plain_bench(bench_pw_convergence
                        pw/si/convergence/bench_ecut_kmesh_convergence.cpp)
    target_compile_definitions(bench_pw_convergence PRIVATE
        ATV_PW_DATA_DIR="${CMAKE_SOURCE_DIR}/benchmarks/pw/si/validation/reference")
    if(ATV_BUILD_CUDA)
        target_compile_definitions(bench_pw_convergence PRIVATE ATV_BUILD_CUDA)
    endif()
endif()
