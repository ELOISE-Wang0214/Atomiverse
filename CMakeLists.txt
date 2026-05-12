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

    set(ATV_REAL "float" CACHE STRING "Floating-point precision (float or double)")
    set_property(CACHE ATV_REAL PROPERTY STRINGS float double)

    set(ATV_INT "int" CACHE STRING "Integer precision (int or int64_t)")
    set_property(CACHE ATV_INT PROPERTY STRINGS int int64_t)

    option(ATV_BUILD_RHF "Build RHF electronic structure code" OFF)
    option(ATV_BUILD_PW "Build plane-wave DFT module (requires FFTW3)" OFF)
    option(ATV_BUILD_TORCH "Build LibTorch backend (requires LibTorch)" OFF)
    option(ATV_BUILD_LIBCINT "Fetch libcint for integral validation tests" OFF)
    set(ATV_MOLRS_PATH "${CMAKE_SOURCE_DIR}/../molrs" CACHE PATH "Path to molrs source tree")
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
        GIT_TAG 0.7.7
    )
    FetchContent_Declare(
        xtensor
        GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
        GIT_TAG 0.25.0
    )
    FetchContent_Declare(
        xtensor-blas
        GIT_REPOSITORY https://github.com/xtensor-stack/xtensor-blas.git
        GIT_TAG 0.21.0
    )
    FetchContent_MakeAvailable(xtl xtensor xtensor-blas)

    # ── Library ──

    # ── OpenMP ──
    find_package(OpenMP)

    # ── MPI (optional) ──
    find_package(MPI QUIET)

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
        # CPU MD components
        src/cpu/md/cell_list.cpp
        src/cpu/md/init.cpp
        src/cpu/md/velocity_verlet.cpp
        src/cpu/md/lj.cpp
        src/cpu/md/nose_hoover.cpp
        src/cpu/observer/thermo_component.cpp
        # CPU bonded components
        src/cpu/bonded/harmonic_bond.cpp
        src/cpu/bonded/harmonic_angle.cpp
        src/cpu/bonded/multi_harmonic_dihedral.cpp
    )

    if(ATV_BUILD_CUDA)
        list(APPEND ATOMIVERSE_SRCS
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
            src/cpu/scf/integrals.cpp
            src/cpu/scf/native_integrals.cpp
            src/cpu/scf/operators.cpp
            src/cpu/scf/diis.cpp
            src/cpu/scf/scf_solve.cpp
            src/cpu/scf/rhf_solve.cpp
            src/cpu/scf/uhf_solve.cpp
            src/cpu/scf/xc_functional.cpp
            src/cpu/scf/grid.cpp
            src/cpu/scf/uks_solve.cpp
            src/cpu/scf/integral_derivatives.cpp
            src/cpu/scf/gradient.cpp
            src/cpu/scf/aspc_predictor.cpp
            src/cpu/scf/rhf_force.cpp
            src/cpu/scf/uhf_force.cpp
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

    if(ATV_BUILD_PW)
        find_library(FFTW3_LIB fftw3 HINTS /opt/homebrew/opt/fftw/lib /usr/local/lib)
        find_path(FFTW3_INC fftw3.h HINTS /opt/homebrew/opt/fftw/include /usr/local/include)
        if(NOT FFTW3_LIB OR NOT FFTW3_INC)
            message(FATAL_ERROR "FFTW3 required for ATV_BUILD_PW but not found. "
                               "Install via: brew install fftw")
        endif()
        add_library(fftw3_stub INTERFACE)
        target_link_libraries(fftw3_stub INTERFACE ${FFTW3_LIB})
        target_include_directories(fftw3_stub INTERFACE ${FFTW3_INC})

        target_sources(atomiverse PRIVATE
            src/cpu/pw/internal/g_vector.cpp
            src/cpu/pw/internal/fft.cpp
            src/cpu/pw/internal/upf_reader.cpp
            src/cpu/pw/internal/kpoint.cpp
            src/cpu/pw/internal/crystal.cpp
            src/cpu/pw/internal/ewald.cpp
            src/cpu/pw/pw_poisson.cpp
            src/cpu/pw/pw_lda.cpp
            src/cpu/pw/pw_hamiltonian.cpp
            src/cpu/pw/pw_veff.cpp
            src/cpu/pw/pw_davidson.cpp
            src/cpu/pw/pw_density.cpp
            src/cpu/pw/pw_convergence.cpp
        )
        target_link_libraries(atomiverse PUBLIC fftw3_stub)
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_PW)
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

    if(ATV_BUILD_CUDA)
        find_package(CUDAToolkit REQUIRED)
    endif()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    if(ATV_BUILD_CUDA)
        target_include_directories(atomiverse PRIVATE src/cuda/md/neighborlist)
    endif()
    target_include_directories(atomiverse PUBLIC include)
    target_link_libraries(atomiverse PUBLIC xtensor xtensor-blas ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    if(ATV_BUILD_CUDA)
        target_link_libraries(atomiverse PUBLIC CUDA::cublas CUDA::cusolver)
        set_target_properties(atomiverse PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_CUDA)
    endif()
    target_compile_definitions(atomiverse PUBLIC ATV_REAL_TYPE=${ATV_REAL} ATV_INT_TYPE=${ATV_INT})

    if(OpenMP_CXX_FOUND)
        target_link_libraries(atomiverse PUBLIC OpenMP::OpenMP_CXX)
    endif()

    if(MPI_CXX_FOUND)
        target_link_libraries(atomiverse PUBLIC MPI::MPI_CXX)
        target_compile_definitions(atomiverse PUBLIC ATV_HAS_MPI)
    endif()

    # ── molrs (hard dependency via CXX FFI) ──
    FetchContent_Declare(Corrosion
        GIT_REPOSITORY https://github.com/corrosion-rs/corrosion
        GIT_TAG v0.5.1
    )
    FetchContent_MakeAvailable(Corrosion)
    corrosion_import_crate(
        MANIFEST_PATH "${ATV_MOLRS_PATH}/Cargo.toml"
        CRATES molcrafts-molrs-cxxapi
    )
    # molrs-core hardcodes F = f64.  Bridge always uses double.
    # Zero-copy when ATV_REAL=double; promotion when ATV_REAL=float.
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
