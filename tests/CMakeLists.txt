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

# MPI test executable: one build target, can be launched at multiple rank
# counts via `atv_add_mpi_run`. The test source provides its own `main`
# that wraps MPI init/finalize.
function(atv_add_mpi_exec target_name)
    add_executable(${target_name} ${ARGN})
    target_link_libraries(${target_name} PRIVATE atomiverse GTest::gtest MPI::MPI_CXX)
    set_target_properties(${target_name} PROPERTIES CXX_STANDARD 17)
    target_include_directories(${target_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

# Register one CTest invocation at a given rank count for an MPI executable.
# Tagged with the "mpi" CTest label so CI can skip them on MPI-less environments.
function(atv_add_mpi_run test_name target_name n)
    add_test(
        NAME ${test_name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${n}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${target_name}> ${MPIEXEC_POSTFLAGS}
    )
    # RUN_SERIAL keeps two MPI tests from running at the same time.
    # RESOURCE_LOCK "mpi_singleton" further blocks concurrent execution with
    # any other MPI test under ctest -jN. Together these prevent OpenMPI
    # session-directory collisions on macOS where two mpirun invocations in
    # the same user account briefly contend for the same prte session id.
    # PASS_REGULAR_EXPRESSION: OpenMPI 5.0.x on macOS occasionally segfaults
    # inside its BTL teardown during MPI_Finalize even when all user tests
    # have completed cleanly; prterun then reports a non-zero exit code. We
    # parse the GoogleTest summary line instead, which appears in output
    # only when every test passed.
    set_tests_properties(${test_name} PROPERTIES
        LABELS "mpi"
        RUN_SERIAL TRUE
        RESOURCE_LOCK "mpi_singleton"
        PASS_REGULAR_EXPRESSION "PASSED[^a-zA-Z]"
        FAIL_REGULAR_EXPRESSION "FAILED[^a-zA-Z]")
endfunction()

# ── CPU backend tests (always available, no CUDA required) ──────────

atv_add_cpp_test(test_cpu_md cpu/md/test_md.cpp)
atv_add_cpp_test(test_cpu_md_pe_single_chain cpu/md/test_pe_single_chain.cpp)
atv_add_cpp_test(test_cpu_md_pe_nve               cpu/md/test_pe_nve.cpp)
atv_add_cpp_test(test_cpu_md_pe_fd_force          cpu/md/test_pe_fd_force.cpp)
atv_add_cpp_test(test_cpu_md_pe_multi_chain       cpu/md/test_pe_multi_chain.cpp)
atv_add_cpp_test(test_cpu_md_lj_exclusions        cpu/md/test_lj_exclusions.cpp)
atv_add_cpp_test(test_cpu_md_pe_nvt               cpu/md/test_pe_nvt.cpp)
atv_add_cpp_test(test_cpu_md_frame_tag_image      cpu/md/test_frame_tag_image.cpp)
atv_add_cpp_test(test_cpu_md_domain_decomposition cpu/md/test_domain_decomposition.cpp)
atv_add_cpp_test(test_cpu_md_mpi_stubs            cpu/md/test_mpi_stubs.cpp)
atv_add_cpp_test(test_cpu_analysis_executor       cpu/test_analysis_executor.cpp)
atv_add_cpp_test(test_cpu_snapshot_publisher      cpu/test_snapshot_publisher.cpp)
atv_add_cpp_test(test_cpu_async_component_base    cpu/test_async_component_base.cpp)
atv_add_cpp_test(test_cpu_observer_rdf_async      cpu/observer/test_rdf_observer.cpp)

# ── MPI tests (launched under `mpirun`) ─────────────────────────────────

if(ATV_USE_MPI AND MPIEXEC_EXECUTABLE)
    atv_add_mpi_exec(test_mpi_topology mpi/test_topology.cpp)
    atv_add_mpi_run(test_mpi_topology_n1 test_mpi_topology 1)
    atv_add_mpi_run(test_mpi_topology_n2 test_mpi_topology 2)
    atv_add_mpi_run(test_mpi_topology_n4 test_mpi_topology 4)
    atv_add_mpi_run(test_mpi_topology_n8 test_mpi_topology 8)

    atv_add_mpi_exec(test_mpi_halo_exchange mpi/test_halo_exchange.cpp)
    atv_add_mpi_run(test_mpi_halo_exchange_n1 test_mpi_halo_exchange 1)
    atv_add_mpi_run(test_mpi_halo_exchange_n2 test_mpi_halo_exchange 2)
    atv_add_mpi_run(test_mpi_halo_exchange_n8 test_mpi_halo_exchange 8)

    atv_add_mpi_exec(test_mpi_lj_nve mpi/test_lj_nve.cpp)
    atv_add_mpi_run(test_mpi_lj_nve_n1 test_mpi_lj_nve 1)
    atv_add_mpi_run(test_mpi_lj_nve_n8 test_mpi_lj_nve 8)

    atv_add_mpi_exec(test_mpi_atom_migration mpi/test_atom_migration.cpp)
    atv_add_mpi_run(test_mpi_atom_migration_n1 test_mpi_atom_migration 1)
    atv_add_mpi_run(test_mpi_atom_migration_n2 test_mpi_atom_migration 2)

    atv_add_mpi_exec(test_mpi_pe_chains_multirank mpi/test_pe_chains_multirank.cpp)
    atv_add_mpi_run(test_mpi_pe_chains_multirank_n2 test_mpi_pe_chains_multirank 2)

    atv_add_mpi_exec(test_mpi_pe_chains_8rank mpi/test_pe_chains_8rank.cpp)
    atv_add_mpi_run(test_mpi_pe_chains_8rank_n8 test_mpi_pe_chains_8rank 8)

    atv_add_mpi_exec(test_mpi_pe_nvt_multirank mpi/test_pe_nvt_multirank.cpp)
    atv_add_mpi_run(test_mpi_pe_nvt_multirank_n8 test_mpi_pe_nvt_multirank 8)
endif()

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
    atv_add_cpp_test(test_cpu_scf_basis_vs_pyscf cpu/scf/test_basis_vs_pyscf.cpp)
    atv_add_cpp_test(test_cpu_scf_def2svp_gradient cpu/scf/test_def2svp_gradient.cpp)
    atv_add_cpp_test(test_cpu_scf_basis_eval_dshell cpu/scf/test_basis_eval_dshell.cpp)
    atv_add_cpp_test(test_cpu_minimizer_lbfgs   cpu/minimizer/test_lbfgs.cpp)
    atv_add_cpp_test(test_cpu_minimizer_recorders cpu/minimizer/test_trajectory_writer.cpp)
    atv_add_cpp_test(test_cpu_scf_orbital_info_observer cpu/scf/test_orbital_info_observer.cpp)
    atv_add_cpp_test(test_cpu_recorder_density_cube cpu/recorder/test_density_cube.cpp)
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
