cmake_minimum_required(VERSION 3.18)
option(ATV_DOCS_ONLY "Configure documentation targets without the C++/CUDA toolchain" OFF)

if(ATV_DOCS_ONLY)
    project(Atomiverse LANGUAGES NONE)
else()
    project(Atomiverse LANGUAGES CXX)
endif()

find_package(Doxygen QUIET)

if(NOT ATV_DOCS_ONLY)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    option(ATV_BUILD_CUDA "Build CUDA backend (requires CUDA toolkit)" ON)
    if(ATV_BUILD_CUDA)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES native)
        endif()
    endif()

    set(ATV_CUDA_REAL "float" CACHE STRING "CUDA per-atom storage precision (float or double)")
    set_property(CACHE ATV_CUDA_REAL PROPERTY STRINGS float double)

    set(ATV_CUDA_WORK "double" CACHE STRING "CUDA kernel computation precision (float or double)")
    set_property(CACHE ATV_CUDA_WORK PROPERTY STRINGS float double)

    set(ATV_CUDA_PREC "double" CACHE STRING "CUDA accumulator/SCF precision (float or double)")
    set_property(CACHE ATV_CUDA_PREC PROPERTY STRINGS float double)

    option(ATV_BUILD_RHF "Build RHF electronic structure code" OFF)
    option(ATV_BUILD_TORCH "Build LibTorch backend (requires LibTorch)" OFF)
    option(ATV_BUILD_LIBCINT "Fetch libcint for integral validation tests" OFF)
    option(ATV_USE_LIBXC "Use libxc for XC ground-truth validation" OFF)
    set(ATV_MOLRS_PATH "${CMAKE_SOURCE_DIR}/../molrs" CACHE PATH "Path to molrs source tree")
    option(ATV_BUILD_PYTHON "Build Python bindings (nanobind)" OFF)
    option(ATV_BUILD_TESTS "Build tests" OFF)
    option(ATV_BUILD_BENCHMARKS "Build benchmarks" OFF)
    option(ATV_BUILD_EXAMPLES "Build examples" OFF)
    option(ATV_USE_MATPLOT "Enable matplotplusplus for plotting in examples (requires gnuplot)" OFF)

    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/deps")

    # ── Dependencies ──

    include(FetchContent)
    FetchContent_Declare(
        xtl
        GIT_REPOSITORY https://github.com/xtensor-stack/xtl.git
        GIT_TAG 0.8.1
    )
    FetchContent_Declare(
        xtensor
        GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
        GIT_TAG 0.27.1
    )
    FetchContent_Declare(
        xtensor-blas
        GIT_REPOSITORY https://github.com/xtensor-stack/xtensor-blas.git
        GIT_TAG 0.21.0
    )
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.15.0
    )
    FetchContent_Declare(
        xsimd
        GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
        GIT_TAG 13.2.0
    )
    FetchContent_MakeAvailable(xtl xtensor xtensor-blas spdlog xsimd)

    # ── spdlog bundled fmt: Apple Clang 21+ consteval workaround ───────
    # spdlog 1.15's bundled fmt base.h defines FMT_USE_CONSTEVAL=1 for
    # Apple Clang >= 14, but Apple Clang 21 has a stricter consteval
    # enforcement that rejects valid fmt code. Patch the version boundary.
    if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        execute_process(COMMAND sed -i ""
            "s/__apple_build_version__ < 14000029L/__apple_build_version__ < 14000029L || __apple_build_version__ >= 21000000L/"
            "${spdlog_SOURCE_DIR}/include/spdlog/fmt/bundled/base.h"
        )
    endif()

    # ── xtensor 0.27.1 compatibility shim for xtensor-blas 0.21.0 ──────
    # xtensor 0.27.1 moved flat headers into subdirectories. xtensor-blas
    # still uses old paths (e.g. <xtensor/xtensor.hpp>). Generate forwarding
    # headers and prepend the compat dir to atomiverse's BUILD_INTERFACE so
    # that xtensor-blas transitive includes resolve correctly at compile time.
    set(_xt_compat "${CMAKE_BINARY_DIR}/xtensor-compat")
    foreach(_old_new
        "xtensor.hpp:containers/xtensor.hpp"
        "xarray.hpp:containers/xarray.hpp"
        "xstorage.hpp:containers/xstorage.hpp"
        "xview.hpp:views/xview.hpp"
        "xstrided_view.hpp:views/xstrided_view.hpp"
        "xnorm.hpp:reducers/xnorm.hpp"
        "xio.hpp:io/xio.hpp"
        "xeval.hpp:core/xeval.hpp"
        "xnoalias.hpp:core/xnoalias.hpp"
        "xcomplex.hpp:misc/xcomplex.hpp"
        "xmanipulation.hpp:misc/xmanipulation.hpp"
        "xutils.hpp:utils/xutils.hpp"
        "xmath.hpp:core/xmath.hpp"
    )
        string(REPLACE ":" ";" _pair "${_old_new}")
        list(GET _pair 0 _old)
        list(GET _pair 1 _new)
        file(WRITE "${_xt_compat}/xtensor/${_old}"
             "#pragma once\n#include <xtensor/${_new}>\n")
    endforeach()

    # ── Library ──

    # ── OpenMP ──
    # Apple Clang doesn't ship OpenMP; use Homebrew libomp if available.
    if(APPLE AND NOT OpenMP_CXX_FOUND)
        find_library(LIBOMP_LIB omp PATHS /opt/homebrew/opt/libomp/lib NO_DEFAULT_PATH)
        if(LIBOMP_LIB)
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "" FORCE)
            set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
            set(OpenMP_omp_LIBRARY "${LIBOMP_LIB}" CACHE STRING "" FORCE)
            include_directories(/opt/homebrew/opt/libomp/include)
        endif()
    endif()
    find_package(OpenMP)

    # ── MPI (optional) ──
    find_package(MPI QUIET)
    option(ATV_USE_MPI "Enable MPI domain decomposition for CPU MD" ${MPI_CXX_FOUND})
    if(ATV_USE_MPI AND NOT MPI_CXX_FOUND)
        message(FATAL_ERROR
            "ATV_USE_MPI=ON but no MPI installation was found. "
            "Install an MPI implementation (OpenMPI, MPICH, MVAPICH) or "
            "configure with -DATV_USE_MPI=OFF.")
    endif()

    set(ATOMIVERSE_SRCS
        src/cpu/scf/solve.cpp
        src/cpu/molrec.cpp
        src/cpu/observer/rdf.cpp
        src/cpu/observer/density.cpp
        src/cpu/recorder/trajectory.cpp
        src/cpu/recorder/thermo.cpp
        src/cpu/recorder/xyz.cpp
        # CPU backend infrastructure
        src/cpu/backend.cpp
        src/cpu/driver.cpp
        src/cpu/analysis_executor.cpp
        src/cpu/async_component_base.cpp
        src/cpu/snapshot_publisher.cpp
        # CPU MD components
        src/cpu/md/cell_list.cpp
        src/cpu/md/init.cpp
        src/cpu/md/velocity_verlet.cpp
        src/cpu/md/lj.cpp
        src/cpu/md/nose_hoover.cpp
        # CPU MD MPI infrastructure (single-rank stubs in Phase 1)
        src/cpu/md/mpi/domain.cpp
        src/cpu/md/mpi/halo_exchange.cpp
        src/cpu/md/mpi/atom_migration.cpp
        src/cpu/md/mpi/reduce.cpp
        src/cpu/observer/thermo_component.cpp
        # CPU bonded components
        src/cpu/bonded/harmonic_bond.cpp
        src/cpu/bonded/harmonic_angle.cpp
        src/cpu/bonded/multi_harmonic_dihedral.cpp
    )

    if(ATV_BUILD_CUDA)
        list(APPEND ATOMIVERSE_SRCS
            src/cuda/log_callback.cpp
            src/cuda/async_bridge.cu
            src/cuda/transfer.cu
            src/cuda/driver.cu
            src/cuda/md/velocity_verlet.cu
            src/cuda/md/init.cu
            src/cuda/md/kinetic_energy.cu
            src/cuda/neighbor_list.cu
            src/cuda/md/neighborlist/sort_and_cluster.cu
            src/cuda/md/neighborlist/grid.cu
            src/cuda/md/neighborlist/pair_search.cu
            src/cuda/md/neighborlist/fused_search.cu
            src/cuda/md/neighborlist/setup.cu
            src/cuda/md/neighborlist/export_pairs.cu
            src/cuda/minimizer/sd.cu
            src/cuda/potentials/lj.cu
            src/cuda/potentials/harmonic_bond.cu
            src/cuda/potentials/harmonic_angle.cu
            src/cuda/potentials/multi_harmonic_dihedral.cu
            src/cuda/md/nose_hoover.cu
            src/cuda/md/deform.cu
            src/cuda/md/barostat.cu
            src/cuda/observer/stress.cu
            src/cuda/scf/hamiltonian.cu
            src/cuda/scf/eigen_solver.cu
            src/cuda/scf/density_builder.cu
            src/cuda/scf/convergence_check.cu
            src/cuda/scf/driver.cu
        )
    endif()

    add_library(atomiverse ${ATOMIVERSE_SRCS})

    if(ATV_BUILD_RHF)
        target_sources(atomiverse PRIVATE
            src/cpu/scf/basis_set.cpp
            src/cpu/scf/one_electron.cpp
            src/cpu/scf/two_electron.cpp
            src/cpu/scf/integrals.cpp
            src/cpu/scf/native_integrals.cpp
            src/cpu/scf/rys_roots.cpp
            src/cpu/scf/rys_eri.cpp
            src/cpu/scf/shell_pair.cpp
            src/cpu/scf/operators.cpp
            src/cpu/scf/diis.cpp
            src/cpu/scf/minao_data.cpp
            src/cpu/scf/minao.cpp
            src/cpu/scf/scf_solve.cpp
            src/cpu/scf/rhf_solve.cpp
            src/cpu/scf/uhf_solve.cpp
            src/cpu/scf/xc_functional.cpp
            src/cpu/scf/grid.cpp
            src/cpu/scf/uks_solve.cpp
            src/cpu/scf/integral_derivatives.cpp
            src/cpu/scf/gradient.cpp
            src/cpu/scf/libcint_bridge.cpp
            src/cpu/scf/aspc_predictor.cpp
            src/cpu/scf/scf_evaluator.cpp
            src/cpu/scf/orbital_info_observer.cpp
            src/cpu/minimizer/init.cpp
            src/cpu/minimizer/lbfgs.cpp
            src/cpu/minimizer/trajectory_writer.cpp
            src/cpu/minimizer/opt_log.cpp
            src/cpu/minimizer/final_geometry.cpp
            src/cpu/recorder/cube.cpp
            src/cpu/recorder/density_cube_writer.cpp
            src/cpu/recorder/orbital_cube_writer.cpp
            src/cpu/observer/mulliken.cpp
        )
        if(ATV_BUILD_CUDA)
            target_sources(atomiverse PRIVATE
                src/cuda/scf/rhf_fock.cu
                src/cuda/scf/rhf_density.cu
                src/cuda/scf/device_integrals.cu
            )
        endif()
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_RHF)
    endif()

    if(ATV_BUILD_TORCH)
        find_package(Torch REQUIRED)
        target_sources(atomiverse PRIVATE
            src/torch/transfer.cu
            src/torch/stream_guard.cu
            src/torch/potentials/lj.cu
            src/torch/neighbor_list.cu
            src/torch/potentials/harmonic_bond.cu
        )
        target_link_libraries(atomiverse PUBLIC ${TORCH_LIBRARIES})
        target_include_directories(atomiverse PUBLIC ${TORCH_INCLUDE_DIRS})
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_TORCH)
    endif()

    # Optional: fetch libcint as a test-only ground-truth reference
    if(ATV_BUILD_LIBCINT)
        FetchContent_Declare(
            libcint
            GIT_REPOSITORY https://github.com/sunqm/libcint.git
            GIT_TAG v6.1.2
        )
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(libcint)
        target_include_directories(cint PUBLIC
            ${libcint_BINARY_DIR}/include
            ${libcint_SOURCE_DIR}/include
        )
        target_link_libraries(atomiverse PUBLIC cint)
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_LIBCINT)
    endif()

    # Optional: fetch libxc for XC ground-truth validation
    if(ATV_USE_LIBXC)
        find_package(Libxc 6 CONFIG QUIET)
        if(Libxc_FOUND)
            message(STATUS "Found system libxc: ${Libxc_DIR}")
            target_link_libraries(atomiverse PUBLIC Libxc::xc)
        else()
            message(STATUS "System libxc not found; fetching from source (tag 6.2.2)")
            set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
            set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
            set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
            set(ENABLE_PYTHON OFF CACHE BOOL "" FORCE)
            set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)
            FetchContent_Declare(
                libxc
                GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
                GIT_TAG 6.2.2
            )
            FetchContent_MakeAvailable(libxc)
            unset(CMAKE_POLICY_VERSION_MINIMUM CACHE)
            target_link_libraries(atomiverse PUBLIC xc)
            # libxc only sets INSTALL_INTERFACE includes; add BUILD_INTERFACE manually
            target_include_directories(atomiverse PUBLIC
                "${libxc_SOURCE_DIR}/src"
                "${libxc_BINARY_DIR}/src"
                "${libxc_BINARY_DIR}"
            )
        endif()
        target_sources(atomiverse PRIVATE
            src/cpu/scf/xc_libxc.cpp
        )
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_LIBXC)
    endif()

    if(ATV_BUILD_CUDA)
        find_package(CUDAToolkit REQUIRED)
    endif()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    if(ATV_BUILD_CUDA)
        target_include_directories(atomiverse PRIVATE src/cuda/md/neighborlist)
    endif()
    target_include_directories(atomiverse PUBLIC include)
    # xtensor-blas 0.21.0 uses old-style flat headers; prepend compat shim.
    # BUILD_INTERFACE: only applies at build time, not to install paths.
    target_include_directories(atomiverse PUBLIC
        $<BUILD_INTERFACE:${_xt_compat}>)
    target_link_libraries(atomiverse PUBLIC xtensor xtensor-blas xsimd spdlog::spdlog ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    # xsimd is linked for direct use in native_integrals.cpp ERI hot loops.
    # XTENSOR_USE_XSIMD is NOT set globally: enabling it triggers SIMD-backed
    # xtensor storage which breaks reshape_view in xc_grid_engine.h.
    # Native CPU tuning: unlock NEON on ARM64, AVX2/AVX-512 on x86-64.
    # Language-guarded to prevent NVCC from seeing -march=native on CUDA files.
    target_compile_options(atomiverse PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
    # Apple Clang 21+ enforces consteval strictly; spdlog 1.15 fmt bundled lib
    # triggers a false consteval-in-non-constexpr-context error. Disable it.
    if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        target_compile_definitions(atomiverse PUBLIC FMT_USE_CONSTEVAL=0)
    endif()
    if(ATV_BUILD_CUDA)
        target_link_libraries(atomiverse PUBLIC CUDA::cuda_driver CUDA::cublas CUDA::cusolver)
        set_target_properties(atomiverse PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_CUDA)
    endif()
    target_compile_definitions(atomiverse PUBLIC ATV_CUDA_REAL_TYPE=${ATV_CUDA_REAL} ATV_CUDA_WORK_TYPE=${ATV_CUDA_WORK} ATV_CUDA_PREC_TYPE=${ATV_CUDA_PREC})

    if(OpenMP_CXX_FOUND)
        target_link_libraries(atomiverse PUBLIC OpenMP::OpenMP_CXX)
    endif()

    if(ATV_USE_MPI)
        target_link_libraries(atomiverse PUBLIC MPI::MPI_CXX)
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_MPI)
    endif()

    # ── molrs (hard dependency via CXX FFI) ──
    FetchContent_Declare(Corrosion
        GIT_REPOSITORY https://github.com/corrosion-rs/corrosion
        GIT_TAG v0.6.1
    )
    FetchContent_MakeAvailable(Corrosion)
    corrosion_import_crate(
        MANIFEST_PATH "${ATV_MOLRS_PATH}/Cargo.toml"
        CRATES molcrafts-molrs-cxxapi
    )
    # molrs-core hardcodes F = f64.  Bridge always uses double.
    # CPU backend uses double; zero-copy across the FFI boundary.
    # Package is `molcrafts-molrs-cxxapi` in Cargo.toml, but `[lib] name` is
    # `molrs_cxxapi` -- corrosion exposes the target under the lib name.
    corrosion_add_cxxbridge(molrs_cxx
        CRATE molrs_cxxapi
        FILES bridge.rs
    )
    target_link_libraries(atomiverse PUBLIC molrs_cxx)

    if(ATV_USE_MATPLOT)
        FetchContent_Declare(
            matplotplusplus
            GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus.git
            GIT_TAG v1.2.1
        )
        FetchContent_MakeAvailable(matplotplusplus)
    endif()

    # ── Python bindings (nanobind) ──
    # When built via scikit-build-core (pip install), nanobind is installed into
    # the build venv before CMake runs, so find_package(nanobind) resolves via
    # `python -m nanobind --cmake_dir`.  For standalone cmake builds, nanobind
    # must be installed in the active Python environment (pip install nanobind).
    if(ATV_BUILD_PYTHON)
        find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)

        # Locate nanobind from the active Python environment
        execute_process(
            COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE
            OUTPUT_VARIABLE NB_CMAKE_DIR
            RESULT_VARIABLE NB_RESULT
        )
        if(NOT NB_RESULT EQUAL 0)
            message(FATAL_ERROR "nanobind not found in Python environment. "
                "Run: pip install nanobind>=2.4")
        endif()
        list(APPEND CMAKE_PREFIX_PATH "${NB_CMAKE_DIR}")
        find_package(nanobind CONFIG REQUIRED)

        nanobind_add_module(_atomiverse
            python/src/bindings/bind.cpp
            python/src/bindings/bind_frame.cpp
            python/src/bindings/bind_md.cpp
            python/src/bindings/bind_scf.cpp
            python/src/bindings/bind_optimizer.cpp
        )
        target_link_libraries(_atomiverse PRIVATE atomiverse)
        target_include_directories(_atomiverse PRIVATE include python/src/bindings)
        install(TARGETS _atomiverse LIBRARY DESTINATION atomiverse)
    endif()

    if(ATV_BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()

    if(ATV_BUILD_BENCHMARKS)
        add_subdirectory(benchmarks)
    endif()

    if(ATV_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()

if(DOXYGEN_FOUND)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating Doxygen documentation"
        VERBATIM
    )

    add_custom_target(docs-open
        COMMAND ${CMAKE_COMMAND}
                -DINPUT_FILE=${CMAKE_SOURCE_DIR}/docs/html/index.html
                -P ${CMAKE_SOURCE_DIR}/cmake/utils/open_generated_file.cmake
        DEPENDS docs
        COMMENT "Opening generated documentation"
        VERBATIM
    )
else()
    message(STATUS "Doxygen not found; docs targets are disabled.")
endif()
