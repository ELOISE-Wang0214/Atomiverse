include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

# CUDA test: compiled as .cu with separable compilation.
function(atv_add_test name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE atomiverse GTest::gtest_main)
    set_target_properties(${name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# C++ test: compiled as .cpp, no CUDA required.
function(atv_add_cpp_test name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE atomiverse GTest::gtest_main)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# ── CPU backend tests (always available, no CUDA required) ──────────

atv_add_cpp_test(test_cpu_md cpu/md/test_md.cpp)

if(ATV_BUILD_RHF)
    atv_add_cpp_test(test_cpu_scf_uhf_solve    cpu/scf/test_uhf_solve.cpp)
    atv_add_cpp_test(test_cpu_scf_uks_solve    cpu/scf/test_uks_solve.cpp)
    atv_add_cpp_test(test_cpu_scf_xc_functional cpu/scf/test_xc_functional.cpp)
    atv_add_cpp_test(test_cpu_scf_grid          cpu/scf/test_grid.cpp)
    atv_add_cpp_test(test_cpu_scf_gradient      cpu/scf/test_gradient.cpp)
    atv_add_cpp_test(test_cpu_scf_evaluator     cpu/scf/test_scf_evaluator.cpp)
    atv_add_cpp_test(test_cpu_scf_integrals     cpu/scf/test_integrals.cpp)
    atv_add_cpp_test(test_cpu_scf_rys_roots     cpu/scf/test_rys_roots.cpp)
    atv_add_cpp_test(test_cpu_scf_rys_eri       cpu/scf/test_rys_eri.cpp)
    atv_add_cpp_test(test_cpu_scf_primitives    cpu/scf/test_primitives.cpp)
    atv_add_cpp_test(test_cpu_scf_accelerators  cpu/scf/test_accelerators.cpp)
    atv_add_cpp_test(test_cpu_scf_minao         cpu/scf/test_minao.cpp)
    atv_add_cpp_test(test_cpu_scf_gradient_vs_libcint cpu/scf/test_gradient_vs_libcint.cpp)
    atv_add_cpp_test(test_cpu_scf_rys_eri_deriv cpu/scf/test_rys_eri_deriv.cpp)
    atv_add_cpp_test(test_cpu_minimizer_lbfgs   cpu/minimizer/test_lbfgs.cpp)
    atv_add_cpp_test(test_xc_reference         test_xc_reference.cpp)
    if(ATV_USE_LIBXC)
        atv_add_cpp_test(test_cpu_scf_xc_vs_libxc cpu/scf/test_xc_vs_libxc.cpp)
    endif()
endif()

# ── CUDA backend tests ──────────────────────────────────────────────

if(ATV_BUILD_CUDA)
    atv_add_test(test_cuda_frame      cuda/test_frame.cu)
    atv_add_test(test_cuda_backend    cuda/test_backend.cu)
    atv_add_test(test_facade_plan     facade/test_plan.cu)
    atv_add_test(test_cuda_md         cuda/md/test_md.cu)
    atv_add_test(test_cuda_md_thermostat cuda/md/test_thermostat.cu)
    atv_add_test(test_cuda_potentials_bonded cuda/potentials/test_bonded.cu)
    atv_add_test(test_cuda_minimizer_sd cuda/minimizer/test_sd.cu)
    atv_add_test(test_cuda_scf        cuda/scf/test_scf.cu)
    atv_add_test(test_cuda_scf_rhf    cuda/scf/test_rhf.cu)
    atv_add_test(test_cuda_scf_device_integrals cuda/scf/test_device_integrals.cu)
endif()

# ── Torch backend tests ─────────────────────────────────────────────

if(ATV_BUILD_TORCH)
    atv_add_test(test_torch_frame          torch/test_frame.cu)
    atv_add_test(test_torch_transfer       torch/test_transfer.cu)
    atv_add_test(test_torch_lj             torch/test_lj.cu)
    atv_add_test(test_torch_harmonic_bond  torch/test_harmonic_bond.cu)
    atv_add_test(test_torch_neighbor_list  torch/test_neighbor_list.cu)
endif()
