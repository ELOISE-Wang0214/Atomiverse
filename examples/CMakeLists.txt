# ── Example layout ──────────────────────────────────────────────────
# Flat, numbered tutorial sequence (see README.md):
#   NN_topic.cpp   CPU canonical example — one concept per file
#   NN_topic.cu    optional CUDA twin of the same number (target NN_topic_cuda)
# Examples demonstrate API usage only; speed / science comparisons
# live in benchmarks/.

# CUDA example (source = <name>.cu, may carry an explicit target name).
function(atv_add_cuda_example name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE atomiverse)
    set_target_properties(${name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    if(ATV_USE_MATPLOT)
        target_link_libraries(${name} PRIVATE Matplot++::matplot)
        target_compile_definitions(${name} PRIVATE ATV_HAS_MATPLOT)
    endif()
endfunction()

# C++ example (no CUDA required).
function(atv_add_cpp_example name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE atomiverse)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17)
    if(ATV_USE_MATPLOT)
        target_link_libraries(${name} PRIVATE Matplot++::matplot)
        target_compile_definitions(${name} PRIVATE ATV_HAS_MATPLOT)
    endif()
endfunction()

# 01 — bare SCF loop on xtensor only (intentionally does not link atomiverse).
add_executable(01_scf 01_scf.cpp)
target_link_libraries(01_scf PRIVATE xtensor xtensor-blas ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
# xtensor-blas still uses old flat-header paths; route through the compat shim
# that atomiverse also uses.
target_include_directories(01_scf PRIVATE $<BUILD_INTERFACE:${_xt_compat}>)
set_target_properties(01_scf PROPERTIES CXX_STANDARD 17)

if(ATV_BUILD_RHF)
    atv_add_cpp_example(02_rhf)
    atv_add_cpp_example(03_uhf)
    atv_add_cpp_example(04_uks)
    atv_add_cpp_example(05_geom_opt)
    atv_add_cpp_example(06_aimd)
    if(ATV_BUILD_CUDA)
        atv_add_cuda_example(02_rhf_cuda 02_rhf.cu)
    endif()
endif()

if(ATV_BUILD_CUDA)
    atv_add_cuda_example(07_md_lj 07_md_lj.cu)
endif()
