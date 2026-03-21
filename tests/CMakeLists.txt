include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

# CUDA test: compiled as .cu with separable compilation.
function(atv_add_test name)
    add_executable(${name} ${name}.cu)
    target_link_libraries(${name} PRIVATE atomiverse GTest::gtest_main)
    set_target_properties(${name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# C++ test: compiled as .cpp, no CUDA required.
function(atv_add_cpp_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE atomiverse GTest::gtest_main)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

if(ATV_BUILD_CUDA)
    atv_add_test(test_frame)
    atv_add_test(test_plan)
    atv_add_test(test_backend)
    atv_add_test(test_md)
    atv_add_test(test_minimizer)
    atv_add_test(test_bonded)
    atv_add_test(test_thermostat)
    atv_add_test(test_scf)
    atv_add_test(test_rhf)
    atv_add_test(test_native_integrals)
    atv_add_test(test_device_integrals)
endif()

# CPU backend tests (always available, no CUDA required)
atv_add_cpp_test(test_cpu_md)

if(ATV_BUILD_RHF)
    atv_add_cpp_test(test_uhf)
    atv_add_cpp_test(test_xc_functional)
    atv_add_cpp_test(test_grid)
    atv_add_cpp_test(test_uks)
    atv_add_cpp_test(test_gradient)
    atv_add_cpp_test(test_aimd)
endif()

if(ATV_BUILD_TORCH)
    atv_add_test(test_torch_frame)
    atv_add_test(test_torch_transfer)
    atv_add_test(test_torch_lj)
    atv_add_test(test_torch_harmonic_bond)
    atv_add_test(test_torch_neighbor_list)
endif()
