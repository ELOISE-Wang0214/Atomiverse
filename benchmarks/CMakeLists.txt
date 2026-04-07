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

atv_add_bench(bench_md_step)
atv_add_bench(bench_scf)
atv_add_bench(bench_rhf)

if(ATV_BUILD_RHF)
    atv_add_cpp_bench(bench_native)
endif()

if(ATV_BUILD_TORCH)
    atv_add_bench(bench_torch)
endif()
