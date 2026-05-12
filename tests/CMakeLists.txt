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

# Core type tests (always available, no CUDA required)
add_executable(test_mat3 core/test_mat3.cpp)
target_link_libraries(test_mat3 PRIVATE atomiverse GTest::gtest_main)
add_test(NAME test_mat3 COMMAND test_mat3)

add_executable(test_vec3 core/test_vec3.cpp)
target_link_libraries(test_vec3 PRIVATE atomiverse GTest::gtest_main)
add_test(NAME test_vec3 COMMAND test_vec3)

add_executable(test_box core/test_box.cpp)
target_link_libraries(test_box PRIVATE atomiverse GTest::gtest_main)
add_test(NAME test_box COMMAND test_box)

add_executable(test_units core/test_units.cpp)
target_link_libraries(test_units PRIVATE atomiverse GTest::gtest_main)
add_test(NAME test_units COMMAND test_units)

# CPU backend tests (always available, no CUDA required)
atv_add_cpp_test(test_cpu_md)

# PW infrastructure tests (require ATV_BUILD_PW)
if(ATV_BUILD_PW)
    add_executable(test_g_vector cpu/pw/test_g_vector.cpp)
    target_link_libraries(test_g_vector PRIVATE atomiverse GTest::gtest_main)
    add_test(NAME test_g_vector COMMAND test_g_vector)

    add_executable(test_fft cpu/pw/test_fft.cpp)
    target_link_libraries(test_fft PRIVATE atomiverse GTest::gtest_main fftw3_stub)
    add_test(NAME test_fft COMMAND test_fft)

    add_executable(test_upf cpu/pw/test_upf.cpp)
    target_link_libraries(test_upf PRIVATE atomiverse GTest::gtest_main)
    add_test(NAME test_upf COMMAND test_upf)

    add_executable(test_kpoint cpu/pw/test_kpoint.cpp)
    target_link_libraries(test_kpoint PRIVATE atomiverse GTest::gtest_main)
    add_test(NAME test_kpoint COMMAND test_kpoint)

    add_executable(test_crystal cpu/pw/test_crystal.cpp)
    target_link_libraries(test_crystal PRIVATE atomiverse GTest::gtest_main)
    add_test(NAME test_crystal COMMAND test_crystal)

    add_executable(test_ewald cpu/pw/test_ewald.cpp)
    target_link_libraries(test_ewald PRIVATE atomiverse GTest::gtest_main)
    add_test(NAME test_ewald COMMAND test_ewald)
endif()

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
