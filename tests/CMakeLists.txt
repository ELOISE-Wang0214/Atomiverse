include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
)
FetchContent_MakeAvailable(googletest)

# Unit tests assert CODE CORRECTNESS on tiny systems (tens to hundreds of
# atoms); performance and LAMMPS-reference validation live in benchmarks/.
# On many-core nodes the OpenMP default (= all hardware threads, e.g. 288 on
# a GH200 node) makes those tiny systems WALL-DOMINATED by fork/join barriers
# — 20-atom 5000-step runs took minutes, and `ctest -j` multiplied the
# oversubscription. Cap every test's default; a test that genuinely needs a
# specific thread count can still call omp_set_num_threads itself.
set(ATV_TEST_OMP_THREADS 8 CACHE STRING
    "OMP_NUM_THREADS for single-process unit tests")
# MPI tests exercise multi-RANK behaviour (domain decomposition / halo), not
# OpenMP — so 1 thread/rank is both sufficient and correct. An n8 run at 4
# threads/rank spawns 8×4=32 contending threads, which on a shared login node
# oversubscribes badly: test_mpi_npt_multirank_n8 took ~103 s AND its
# thread-reduction-order-sensitive pressure check flaked. At 1 thread/rank it
# passes deterministically in ~20 s.
set(ATV_TEST_MPI_OMP_THREADS 1 CACHE STRING
    "OMP_NUM_THREADS per rank for MPI unit tests")

# Fast-gate discipline: this suite is unit tests only — each exercises one
# function or method at millisecond scale. End-to-end runs (a full SCF/MD/vib
# solve to convergence, a multi-sweep MC equilibration, an external-engine
# parity check) do not belong here at any granularity; their home is
# benchmarks/<domain>/validation, run at acceptance time. There is deliberately
# no opt-in switch to run e2e from ctest.

# CUDA test: compiled as .cu with separable compilation.
# A .cu source registered while the CUDA language is disabled is silently
# dropped by CMake — the binary links, runs zero tests, and reports PASS.
# Refuse the configuration instead of shipping a green that tested nothing.
function(atv_guard_cuda_sources name)
    foreach(src ${ARGN})
        if(src MATCHES "\\.cu$" AND NOT ATV_BUILD_CUDA)
            message(FATAL_ERROR
                "test '${name}' has a CUDA source (${src}) but ATV_BUILD_CUDA "
                "is OFF: CMake would silently drop it and the test binary "
                "would report PASS with zero tests. Enable CUDA or do not "
                "register this test.")
        endif()
    endforeach()
endfunction()

function(atv_add_test name)
    atv_guard_cuda_sources(${name} ${ARGN})
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE atomiverse GTest::gtest_main)
    set_target_properties(${name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=${ATV_TEST_OMP_THREADS}")
endfunction()

# C++ test: compiled as .cpp, no CUDA required.
function(atv_add_cpp_test name)
    atv_guard_cuda_sources(${name} ${ARGN})
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE atomiverse GTest::gtest_main)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 20)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=${ATV_TEST_OMP_THREADS}")
endfunction()

# MPI test executable: one build target, can be launched at multiple rank
# counts via `atv_add_mpi_run`. The test source provides its own `main`
# that wraps MPI init/finalize.
function(atv_add_mpi_exec target_name)
    add_executable(${target_name} ${ARGN})
    target_link_libraries(${target_name} PRIVATE atomiverse GTest::gtest MPI::MPI_CXX)
    set_target_properties(${target_name} PROPERTIES CXX_STANDARD 20)
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
    set_tests_properties(${test_name} PROPERTIES
        LABELS "mpi"
        RUN_SERIAL TRUE
        RESOURCE_LOCK "mpi_singleton"
        ENVIRONMENT "OMP_NUM_THREADS=${ATV_TEST_MPI_OMP_THREADS}"
        PASS_REGULAR_EXPRESSION "PASSED[^a-zA-Z]"
        FAIL_REGULAR_EXPRESSION "FAILED[^a-zA-Z]")
endfunction()

# ── Core type tests (always available, no CUDA required) ──────────────

atv_add_cpp_test(test_mat3  core/test_mat3.cpp)
atv_add_cpp_test(test_vec3  core/test_vec3.cpp)
atv_add_cpp_test(test_real_harmonics core/test_real_harmonics.cpp)
atv_add_cpp_test(test_units core/test_units.cpp)

add_test(
    NAME test_molrs_source_contract_mismatch
    COMMAND ${CMAKE_COMMAND}
            -DATV_SOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/molrs-contract-mismatch
            -P ${CMAKE_CURRENT_SOURCE_DIR}/architecture/molrs_contract_mismatch.cmake)
set_tests_properties(test_molrs_source_contract_mismatch PROPERTIES
    WILL_FAIL TRUE)

# ── CPU backend tests (always available, no CUDA required) ──────────

# Mirrors include/atomiverse/cpu/frame/box.h (+ box_simd.h).
atv_add_cpp_test(test_cpu_box cpu/frame/test_box.cpp)
# Context service slices (arch-refactor-05-context): flat core fields + nested
# md/mc/opt/analysis pure-aggregate slices, exact field mapping, all-null
# defaults, trivially-copyable / non-polymorphic guards. Layout-only, no
# physics.
atv_add_cpp_test(test_cpu_component_context_slices cpu/component/test_context_slices.cpp)
# esmm-contract-05-accumulate-switch (U8 + U7b): the kEnergyElectronic per-term
# key — ZeroEnergy must clear it in lockstep with kEnergy, and a co-bound
# Evaluate pass must leave the QM share alone in it. Registered UNCONDITIONALLY
# (ZeroEnergy is always-on infrastructure); the two co-binding cases inside are
# #ifdef ATV_HAS_HF. Serial: the RHF pin must not drift with reduction order.
atv_add_cpp_test(test_cpu_zero_energy_electronic cpu/component/test_zero_energy_electronic.cpp)
set_tests_properties(test_cpu_zero_energy_electronic PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
atv_add_cpp_test(test_cpu_mc_scaffold cpu/test_mc_scaffold.cpp)
atv_add_cpp_test(test_cpu_mc_move cpu/montecarlo/test_mc_move.cpp)
atv_add_cpp_test(test_cpu_mc_transaction_rollback cpu/montecarlo/test_transaction_rollback.cpp)
atv_add_cpp_test(test_cpu_mc_ensemble_policies cpu/montecarlo/test_ensemble_policies.cpp)
atv_add_cpp_test(test_cpu_mc_gcmc_npt_movers cpu/montecarlo/test_gcmc_npt_movers.cpp)
atv_add_cpp_test(test_cpu_mc_row_order cpu/montecarlo/test_mc_row_order.cpp)
# The MC physics-validation gates (pair_distribution / acceptance_ratio /
# lj_nist) drive MC::run over many sweeps — e2e, not unit. They now live in
# benchmarks/mc/lj/validation/ and run at acceptance time, never in ctest.
atv_add_cpp_test(test_cpu_md cpu/dynamics/test_md.cpp)
atv_add_cpp_test(test_cpu_md_parrinello_rahman_npt cpu/dynamics/test_parrinello_rahman_npt.cpp)
atv_add_cpp_test(test_cpu_md_nph_integrator        cpu/dynamics/test_nph_integrator.cpp)
atv_add_cpp_test(test_cpu_md_surface_tension_npt   cpu/dynamics/test_surface_tension_npt.cpp)
atv_add_cpp_test(test_cpu_md_pe_thermo_schema cpu/dynamics/test_pe_thermo_schema.cpp)
atv_add_cpp_test(test_cpu_md_pe_single_chain cpu/dynamics/test_pe_single_chain.cpp)
# Multi-step PE bonded NVE conservation is a regression
# (regressions/md/pe_bonded_nve), not a unit test.
# test_cpu_md_csvr / test_cpu_md_nvt relocated to benchmarks/md/lj/validation/
# (val_md_csvr / val_md_nvt) — full thermostat equilibrations, not ms unit tests.
atv_add_cpp_test(test_cpu_md_atom_sort            cpu/dynamics/test_atom_sort.cpp)
atv_add_cpp_test(test_cpu_md_step_timing          cpu/dynamics/test_step_timing.cpp)
atv_add_cpp_test(test_cpu_md_timing_report        cpu/dynamics/test_timing_report.cpp)
atv_add_cpp_test(test_cpu_md_pe_fd_force          cpu/dynamics/test_pe_fd_force.cpp)
# Multi-chain multi-step NVE was e2e; multi-chain unit surface is
# special-bonds / net-force one-launch only in pe_multi_chain.
atv_add_cpp_test(test_cpu_md_pe_multi_chain       cpu/dynamics/test_pe_multi_chain.cpp)
atv_add_cpp_test(test_cpu_md_lj_exclusions        cpu/dynamics/test_lj_exclusions.cpp)
atv_add_cpp_test(test_cpu_md_pair_lj_cut          cpu/dynamics/test_pair_lj_cut.cpp)
atv_add_cpp_test(test_cpu_md_const_pot_bo         cpu/dynamics/test_const_pot_bo.cpp)
atv_add_cpp_test(test_cpu_md_const_charge_electrode cpu/dynamics/test_const_charge_electrode.cpp)
atv_add_cpp_test(test_cpu_md_qeq                  cpu/dynamics/test_qeq.cpp)
atv_add_cpp_test(test_cpu_md_const_pot_cg         cpu/dynamics/test_const_pot_cg.cpp)
atv_add_cpp_test(test_cpu_md_const_pot_per_eta    cpu/dynamics/test_const_pot_per_eta.cpp)
atv_add_cpp_test(test_cpu_md_charge_updater_probe cpu/dynamics/test_charge_updater_probe.cpp)
atv_add_cpp_test(test_cpu_comm_halo_refresh_charges cpu/communication/test_halo_refresh_charges.cpp)
atv_add_cpp_test(test_cpu_comm_tag_index_lifetime cpu/communication/test_tag_index_lifetime.cpp)
atv_add_cpp_test(test_cpu_md_const_pot_ffield     cpu/dynamics/test_const_pot_ffield.cpp)
atv_add_cpp_test(test_cpu_md_const_pot_maze_symplectic cpu/dynamics/test_const_pot_maze_symplectic.cpp)
atv_add_cpp_test(test_cpu_md_const_pot_site_volume cpu/dynamics/test_const_pot_site_volume.cpp)
atv_add_cpp_test(test_cpu_voronoi_volume          cpu/observe/test_voronoi_volume.cpp)
atv_add_cpp_test(test_cpu_md_dipole_electrode_joint cpu/dynamics/test_dipole_electrode_joint.cpp)
atv_add_cpp_test(test_cpu_classical_polarization  cpu/observe/test_classical_polarization.cpp)
atv_add_cpp_test(test_cpu_pair_coul_cut           cpu/potentials/test_pair_coul_cut.cpp)
atv_add_cpp_test(test_cpu_charging_work           cpu/electrodes/test_charging_work.cpp)
atv_add_cpp_test(test_cpu_electrode_charge_aspc   cpu/electrodes/test_charge_aspc.cpp)
atv_add_cpp_test(test_cpu_pair_lj_cut_tip4p_cut   cpu/potentials/test_pair_lj_cut_tip4p_cut.cpp)
atv_add_cpp_test(test_cpu_pair_lj_soft            cpu/potentials/test_pair_lj_soft.cpp)
atv_add_cpp_test(test_cpu_pair_coul_soft          cpu/potentials/test_pair_coul_soft.cpp)
atv_add_cpp_test(test_cpu_lambda_ti               cpu/potentials/test_lambda_ti.cpp)
atv_add_cpp_test(test_cpu_pair_special_scale      cpu/potentials/test_pair_special_scale.cpp)
atv_add_cpp_test(test_cpu_coul_long_special       cpu/potentials/test_coul_long_special.cpp)
atv_add_cpp_test(test_cpu_pair_lj_cut_coul_long   cpu/potentials/test_pair_lj_cut_coul_long.cpp)
atv_add_cpp_test(test_cpu_dihedral_charmm         cpu/potentials/test_dihedral_charmm.cpp)
atv_add_cpp_test(test_cpu_improper_harmonic       cpu/potentials/test_improper_harmonic.cpp)
atv_add_cpp_test(test_cpu_compiled_forcefield     cpu/io/test_compiled_forcefield.cpp)
atv_add_cpp_test(test_cpu_lammps_data_full        cpu/io/test_lammps_data_full.cpp)
atv_add_cpp_test(test_cpu_lammps_real             cpu/io/test_lammps_real.cpp)
atv_add_cpp_test(test_cpu_fix_efield              cpu/dynamics/test_fix_efield.cpp)
atv_add_cpp_test(test_cpu_fix_dfield              cpu/dynamics/test_fix_dfield.cpp)
atv_add_cpp_test(test_cpu_fix_addforce            cpu/dynamics/test_fix_addforce.cpp)
atv_add_cpp_test(test_cpu_fix_oscillating_efield  cpu/dynamics/test_fix_oscillating_efield.cpp)
atv_add_cpp_test(test_cpu_fix_wall_lj1043         cpu/dynamics/test_fix_wall_lj1043.cpp)
atv_add_cpp_test(test_cpu_bias_distance_cv        cpu/bias/test_distance_cv.cpp)
atv_add_cpp_test(test_cpu_bias_angle_cv           cpu/bias/test_angle_cv.cpp)
atv_add_cpp_test(test_cpu_bias_dihedral_cv        cpu/bias/test_dihedral_cv.cpp)
atv_add_cpp_test(test_cpu_bias_metadynamics_1d    cpu/bias/test_metadynamics_1d.cpp)
atv_add_cpp_test(test_cpu_bias_metadynamics_nd    cpu/bias/test_metadynamics_nd.cpp)
atv_add_cpp_test(test_cpu_wham                    cpu/observe/test_wham.cpp)
atv_add_cpp_test(test_cpu_pair_yukawa             cpu/potentials/test_pair_yukawa.cpp)
atv_add_cpp_test(test_cpu_gaussian_charge         cpu/potentials/test_gaussian_charge.cpp)
atv_add_cpp_test(test_cpu_pair_buck               cpu/potentials/test_pair_buck.cpp)
atv_add_cpp_test(test_cpu_pair_born               cpu/potentials/test_pair_born.cpp)
atv_add_cpp_test(test_cpu_pair_virial             cpu/potentials/test_pair_virial.cpp)
atv_add_cpp_test(test_cpu_pim_induced_dipole      cpu/potentials/test_pim_induced_dipole.cpp)
atv_add_cpp_test(test_cpu_pim_xl_dipole           cpu/potentials/test_pim_xl_dipole.cpp)
atv_add_cpp_test(test_cpu_minimize_aniso_cell     cpu/minimize/test_anisotropic_cell.cpp)
atv_add_cpp_test(test_cpu_xtb_gfn0                cpu/semiempirical/test_xtb_gfn0.cpp)
# QM/MM embedding units. src/cpu/multiscale/*.cpp compiles only under
# ATV_BUILD_SEMIEMPIRICAL (root CMakeLists), and the electrostatic /
# polarizable objects additionally reference ATV::cpu::scf::RHF, which only
# links under ATV_BUILD_HF. Registering them unconditionally broke the default
# configuration (HF=OFF) at link time.
if(ATV_BUILD_SEMIEMPIRICAL)
    atv_add_cpp_test(test_cpu_qmmm_mechanical     cpu/multiscale/test_qmmm_mechanical.cpp)
    if(ATV_BUILD_HF)
        atv_add_cpp_test(test_cpu_qmmm_electrostatic  cpu/multiscale/test_qmmm_electrostatic.cpp)
        # test_link_atom includes multiscale/qmmm_electrostatic.h and drives the
        # RHF-embedded link-atom cut — same RHF link dependency.
        atv_add_cpp_test(test_cpu_link_atom           cpu/multiscale/test_link_atom.cpp)
        atv_add_cpp_test(test_cpu_qmmm_polarizable    cpu/multiscale/test_qmmm_polarizable.cpp)
    endif()
endif()
atv_add_cpp_test(test_cpu_minimize_triclinic      cpu/minimize/test_triclinic_cell.cpp)
atv_add_cpp_test(test_cpu_bond_fene               cpu/potentials/test_bond_fene.cpp)
atv_add_cpp_test(test_cpu_pair_wca                cpu/potentials/test_pair_wca.cpp)
atv_add_cpp_test(test_cpu_pair_dpd                cpu/potentials/test_pair_dpd.cpp)
atv_add_cpp_test(test_cpu_path_hremd              cpu/path/test_hremd.cpp)
atv_add_cpp_test(test_cpu_ewald_sum               cpu/longrange/test_ewald_sum.cpp)
atv_add_cpp_test(test_cpu_pme_math                cpu/longrange/test_pme_math.cpp)
atv_add_cpp_test(test_cpu_ewald_coul              cpu/longrange/test_ewald_coul.cpp)
atv_add_cpp_test(test_cpu_slab_ewald              cpu/longrange/test_slab_ewald.cpp)
# Backend-agnostic FMM operator math (absorbed into cpu/longrange/, always on)
atv_add_cpp_test(test_multipole_expansion cpu/longrange/test_multipole_expansion.cpp)
atv_add_cpp_test(test_fmm_operators cpu/longrange/test_fmm_operators.cpp)
atv_add_cpp_test(test_fmm_tree cpu/longrange/test_fmm_tree.cpp)
# hf-capability-gaps-21-cfmm: WJGH Gaussian extent (header-only facade)
atv_add_cpp_test(test_gaussian_extent cpu/longrange/test_gaussian_extent.cpp)
if(ATV_BUILD_KSPACE)
    atv_add_cpp_test(test_cpu_pme_coul            cpu/longrange/test_pme_coul.cpp)
    atv_add_cpp_test(test_cpu_pppm_coul           cpu/longrange/test_pppm_coul.cpp)
endif()
atv_add_cpp_test(test_cpu_nddo_api_surface       cpu/semiempirical/test_nddo_api_surface.cpp)
atv_add_cpp_test(test_cpu_am1_scaffold            cpu/semiempirical/test_am1_scaffold.cpp)
atv_add_cpp_test(test_cpu_rm1                     cpu/semiempirical/test_rm1.cpp)
atv_add_cpp_test(test_cpu_dftb2                   cpu/semiempirical/test_dftb2.cpp)
atv_add_cpp_test(test_cpu_dftb0                   cpu/semiempirical/test_dftb0.cpp)
atv_add_cpp_test(test_cpu_dftb2_sk                cpu/semiempirical/test_dftb2_sk.cpp)
atv_add_cpp_test(test_cpu_dftb2_3ob               cpu/semiempirical/test_dftb2_3ob.cpp)
atv_add_cpp_test(test_cpu_dftb2_dispersion_d3     cpu/semiempirical/test_dftb2_dispersion_d3.cpp)
atv_add_cpp_test(test_cpu_dftb3                   cpu/semiempirical/test_dftb3.cpp)
atv_add_cpp_test(test_cpu_dftb3_scc               cpu/semiempirical/test_dftb3_scc.cpp)
atv_add_cpp_test(test_cpu_scc_mixer               cpu/semiempirical/test_scc_mixer.cpp)

# GPW module tests. All of them run unconditionally on the in-tree GTH + MOLOPT +
# D3 tables (cpu/gpw/basis_data/) — no external data directory is consulted at
# build or test time. (Until 2026-07-11 the 14 tests that touched a real molecule
# only registered when ATV_CP2K_DATA_DIR pointed at a CP2K install, so with the
# gate off GPW had no real-molecule physics and no force validation at all.)
if(ATV_BUILD_GPW)
    # Every GPW test below runs on the in-tree data tables
    # (cpu/gpw/basis_data/gpw_data.h) — no external data directory, no
    # third-party file read. E2E (converged SCF / AIMD / geoopt / anchors) lives
    # in benchmarks/qc/gpw/validation and never runs under ctest.
    atv_add_cpp_test(test_cpu_gpw_gth_integrals    cpu/gpw/test_cpu_gpw_gth_integrals.cpp)
    atv_add_cpp_test(test_cpu_gpw_collocate        cpu/gpw/test_cpu_gpw_collocate.cpp)
    atv_add_cpp_test(test_cpu_gpw_multigrid        cpu/gpw/test_cpu_gpw_multigrid.cpp)
    if(ATV_BUILD_PW)  # FFT transfer (prolong_fft / restrict_fft) needs FFTW3
        target_link_libraries(test_cpu_gpw_multigrid PRIVATE fftw3_stub)
        target_compile_definitions(test_cpu_gpw_multigrid PRIVATE ATV_GPW_HAS_FFT)
    endif()
    atv_add_cpp_test(test_cpu_gpw_xc_grid          cpu/gpw/test_cpu_gpw_xc_grid.cpp)
    atv_add_cpp_test(test_cpu_gpw_one_electron     cpu/gpw/test_cpu_gpw_one_electron.cpp)
    atv_add_cpp_test(test_cpu_xc_pade              cpu/xc/test_cpu_xc_pade.cpp)
    if(ATV_BUILD_PW)  # FFT-backed pieces
        atv_add_cpp_test(test_cpu_gpw_poisson_fft  cpu/gpw/test_cpu_gpw_poisson_fft.cpp)
        atv_add_cpp_test(test_cpu_gpw_hartree_grid cpu/gpw/test_cpu_gpw_hartree_grid.cpp)
        atv_add_cpp_test(test_cpu_gpw_grid_fft     cpu/gpw/test_cpu_gpw_grid_fft.cpp)
        atv_add_cpp_test(test_cpu_gpw_xc_grid_gga  cpu/gpw/test_cpu_gpw_xc_grid_gga.cpp)
        atv_add_cpp_test(test_cpu_gpw_mt_poisson   cpu/gpw/test_cpu_gpw_mt_poisson.cpp)
        atv_add_cpp_test(test_cpu_gpw_method       cpu/gpw/test_gpw_method.cpp)
        atv_add_cpp_test(test_cpu_gpw_solver       cpu/gpw/test_cpu_gpw_solver.cpp)
        atv_add_cpp_test(test_cpu_gpw_kpoints      cpu/gpw/test_cpu_gpw_kpoints.cpp)
        atv_add_cpp_test(test_cpu_gpw_periodic_wannier
                                                   cpu/gpw/test_cpu_gpw_periodic_wannier.cpp)
        atv_add_cpp_test(test_cpu_gpw_smearing     cpu/gpw/test_cpu_gpw_smearing.cpp)
        atv_add_cpp_test(test_cpu_gpw_ot           cpu/gpw/test_cpu_gpw_ot.cpp)
        # Γ-point GPW MP2 kernels (spin-orbital cross-check, MO-on-grid,
        # (ia|jb) vs Poisson Hartree identity, guard throws). No SCF here.
        atv_add_cpp_test(test_cpu_gpw_mp2          cpu/gpw/test_cpu_gpw_mp2.cpp)
        # Γ-point GPW dRPA plasmon kernel + GpwDft::direct_rpa guards.
        atv_add_cpp_test(test_cpu_gpw_rpa          cpu/gpw/test_cpu_gpw_rpa.cpp)
        atv_add_cpp_test(test_cpu_gpw_g0w0         cpu/gpw/test_cpu_gpw_g0w0.cpp)
        # XL-BOMD propagator (Table-I pin, fixed point, bounded response)
        # + GpwDft guess_P plumbing on the inline ms-scale H2.
        atv_add_cpp_test(test_cpu_gpw_xl_bomd     cpu/gpw/test_gpw_xl_bomd.cpp)
        atv_add_cpp_test(test_cpu_gpw_cpmd        cpu/gpw/test_gpw_cpmd.cpp)
        atv_add_cpp_test(test_cpu_gpw_second_gen_cpmd cpu/gpw/test_gpw_second_gen_cpmd.cpp)
        if(ATV_BUILD_MD)
            atv_add_cpp_test(test_cpu_gpw_bomd    cpu/gpw/test_gpw_bomd.cpp)
        endif()
    endif()

    # Data-table unit tests: grammar/ingest of the in-tree GTH + MOLOPT tables and
    # the D3 parameter tables. No SCF is run here — the converged-SCF checks that
    # used to live alongside these are now benchmarks/qc/gpw/validation.
    atv_add_cpp_test(test_cpu_gpw_data_table    cpu/gpw/test_cpu_gpw_data_table.cpp)
    # Minimal model system for the Γ-point lattice sum: S must be positive definite
    # at every cell size (Gram matrix). Pinned against PySCF.
    atv_add_cpp_test(test_cpu_gpw_periodic_overlap cpu/gpw/test_cpu_gpw_periodic_overlap.cpp)
    # Heavy elements: Zn/Na/Cl GTH-PBE +
    # DZVP-MOLOPT-SR-GTH ingest, incl. the Zn d/f channels.
    atv_add_cpp_test(test_cpu_gpw_heavy_elements cpu/gpw/test_cpu_gpw_heavy_elements.cpp)
    # Grimme DFT-D3 (D3): energy + analytic force incl. the CN-derivative
    # term, against a finite-difference check. Reference-free.
    atv_add_cpp_test(test_cpu_gpw_dftd3         cpu/gpw/test_cpu_gpw_dftd3.cpp)
    atv_add_cpp_test(test_cpu_gpw_dftd4         cpu/gpw/test_cpu_gpw_dftd4.cpp)
    atv_add_cpp_test(test_cpu_gpw_ts            cpu/gpw/test_cpu_gpw_ts.cpp)
    atv_add_cpp_test(test_cpu_gpw_xdm           cpu/gpw/test_cpu_gpw_xdm.cpp)
    atv_add_cpp_test(test_cpu_gpw_mbd           cpu/gpw/test_cpu_gpw_mbd.cpp)
    atv_add_cpp_test(test_cpu_gpw_vdw_df        cpu/gpw/test_cpu_gpw_vdw_df.cpp)
    atv_add_cpp_test(test_cpu_gpw_rvv10         cpu/gpw/test_cpu_gpw_rvv10.cpp)
    atv_add_cpp_test(test_cpu_gpw_hirshfeld     cpu/gpw/test_cpu_gpw_hirshfeld.cpp)
endif()
atv_add_cpp_test(test_cpu_md_pair_morse_smoke     cpu/dynamics/test_pair_morse_smoke.cpp)
atv_add_cpp_test(test_cpu_md_pe_ff_constants      cpu/dynamics/test_pe_ff_constants.cpp)
atv_add_cpp_test(test_cpu_md_pe_energy_xcheck     cpu/dynamics/test_pe_energy_xcheck.cpp)
# pe_tension_protocol grepped benchmarks/md/pe/... which was deleted; re-add
# when the PE tension app is restored under regressions/ or examples/.
atv_add_cpp_test(test_cpu_md_neighbor_list_arena  cpu/dynamics/test_neighbor_list_arena.cpp)
atv_add_cpp_test(test_cpu_md_neighbor_list_rebuild  cpu/dynamics/test_neighbor_list_rebuild.cpp)
atv_add_cpp_test(test_cpu_md_full_neighbor_list  cpu/dynamics/test_full_neighbor_list.cpp)
atv_add_cpp_test(test_cpu_md_cell_list_grid_stale  cpu/dynamics/test_cell_list_grid_stale.cpp)
atv_add_cpp_test(test_cpu_atom_variable_n  cpu/test_atom_variable_n.cpp)
atv_add_cpp_test(test_cpu_frame_registry  cpu/frame/test_frame_registry.cpp)
atv_add_cpp_test(test_cpu_eam  cpu/potentials/test_eam.cpp)
atv_add_cpp_test(test_cpu_sw  cpu/potentials/test_sw.cpp)
if(ATV_BUILD_VIB)
    atv_add_cpp_test(test_cpu_vib_phonon_supercell cpu/vibration/test_phonon_supercell.cpp)
    atv_add_cpp_test(test_cpu_vib_phonon_dispersion cpu/vibration/test_phonon_dispersion.cpp)
    # esmm-contract-01-zeroing: the Hessian driver must zero kAtom_f* / kEnergy
    # per displacement pass, or an accumulating source halves the curvature.
    atv_add_cpp_test(test_cpu_vib_hessian_zeroing
        cpu/vibration/test_hessian_zeroing.cpp)
    # GPW force-consumer wiring for Phonon (no NeighborList). Full FD spectrum
    # lives in regressions/gpw/gpw-fd-phonon-si (6N SCF — not unit budget).
    if(ATV_BUILD_GPW)
        atv_add_cpp_test(test_cpu_vib_phonon_gpw_wiring
            cpu/vibration/test_phonon_gpw_wiring.cpp)
    endif()
endif()
atv_add_cpp_test(test_cpu_pair_lj_types  cpu/potentials/test_pair_lj_types.cpp)
atv_add_cpp_test(test_cpu_md_fix_rattle  cpu/dynamics/test_fix_rattle.cpp)
atv_add_cpp_test(test_cpu_md_fix_freeze  cpu/dynamics/test_fix_freeze.cpp)
atv_add_cpp_test(test_cpu_md_piston_electrode cpu/dynamics/test_piston_electrode.cpp)
atv_add_cpp_test(test_cpu_md_const_potential_bo  cpu/dynamics/test_const_potential_bo.cpp)
atv_add_cpp_test(test_cpu_md_fix_settle  cpu/dynamics/test_fix_settle.cpp)
atv_add_cpp_test(test_cpu_md_virtual_sites  cpu/dynamics/test_virtual_sites.cpp)
atv_add_cpp_test(test_cpu_md_brownian     cpu/dynamics/test_brownian.cpp)
atv_add_cpp_test(test_cpu_md_write_sets  cpu/dynamics/test_write_sets.cpp)
atv_add_cpp_test(test_cpu_path_replica    cpu/path/test_replica.cpp)
atv_add_cpp_test(test_cpu_path_remd       cpu/path/test_remd.cpp)
atv_add_cpp_test(test_cpu_path_pimd       cpu/path/test_pimd.cpp)
atv_add_cpp_test(test_cpu_path_dimer      cpu/path/test_dimer.cpp)
atv_add_cpp_test(test_cpu_path_string     cpu/path/test_string_method.cpp)
atv_add_cpp_test(test_cpu_path_irc        cpu/path/test_irc.cpp)
atv_add_cpp_test(test_cpu_path_rpmd_est   cpu/path/test_rpmd_estimator.cpp)
# esmm-contract-01-zeroing: PathDriver must zero kAtom_f* / kEnergy before each
# per-replica force pass, or an accumulating source grows with the iteration.
atv_add_cpp_test(test_cpu_path_zeroing    cpu/path/test_path_zeroing.cpp)
# esmm-contract-02-validation-walks: PathDriver::validate_components must see a
# nested `Driver<scf::Phase> + HF` method node as the band's force source.
atv_add_cpp_test(test_cpu_path_validation cpu/path/test_path_validation.cpp)
atv_add_cpp_test(test_cpu_minimizer_joint_cell cpu/minimize/test_joint_cell.cpp)
# esmm-contract-03-runslot-nested: Minimize's StrongWolfe trials must re-enter a
# nested `Driver<scf::Phase> + HF` force source (run_stage → Runtime::run_slot_).
# Stub physics (harmonic well) — no SCF, milliseconds.
atv_add_cpp_test(test_cpu_minimizer_line_search_nested_node
                 cpu/minimize/test_line_search_nested_node.cpp)
atv_add_cpp_test(test_cpu_criterion       cpu/driver/test_criterion.cpp)
atv_add_cpp_test(test_cpu_driver_nesting  cpu/driver/test_driver_nesting.cpp)
# scf-driver-refactor-01-plumb: scf::Phase / PhaseKind::SCF / PhaseTraits<scf::Phase>
# / the Criterion component kind + driver pull-poll. Pure type-foundation plumbing —
# no consumer yet.
atv_add_cpp_test(test_cpu_scf_phase_plumbing cpu/driver/test_scf_phase_plumbing.cpp)
atv_add_cpp_test(test_cpu_component_ownership cpu/driver/test_component_ownership.cpp)
# esmm-contract-04-forces-requested: `forces_requested` as a compiled Plan-slot
# property — MC's energy phase is energy-only, MD/Minimize/Evaluate still ask for
# forces, and the veto ANDs through a nested driver's own slots. Stubs only, no
# physics, microseconds; the QM-in-MC e2e is regressions/mc/mc_bias_es.
atv_add_cpp_test(test_cpu_forces_requested cpu/driver/test_forces_requested.cpp)
# test_cpu_md_pe_nvt relocated to benchmarks/md/pe/validation/ (val_md_pe_nvt) —
# a full NVT equilibration of the bonded PE melt, not a ms unit test.
atv_add_cpp_test(test_cpu_md_frame_tag_image      cpu/dynamics/test_frame_tag_image.cpp)
atv_add_cpp_test(test_cpu_md_domain_decomposition cpu/dynamics/test_domain_decomposition.cpp)
atv_add_cpp_test(test_cpu_md_mpi_stubs            cpu/dynamics/test_mpi_stubs.cpp)
atv_add_cpp_test(test_cpu_analysis_executor       cpu/test_analysis_executor.cpp)
atv_add_cpp_test(test_cpu_snapshot_publisher      cpu/test_snapshot_publisher.cpp)
atv_add_cpp_test(test_cpu_async_component_base    cpu/test_async_component_base.cpp)
atv_add_cpp_test(test_cpu_observer_rdf_async      cpu/observe/test_rdf_observer.cpp)
atv_add_cpp_test(test_cpu_observer_msd            cpu/observe/test_msd.cpp)
atv_add_cpp_test(test_cpu_observer_vacf           cpu/observe/test_vacf.cpp)
atv_add_cpp_test(test_cpu_observer_widom          cpu/observe/test_widom.cpp)
atv_add_cpp_test(test_cpu_observer_electrode_charge cpu/observe/test_electrode_charge_observer.cpp)
atv_add_cpp_test(test_cpu_observer_capacitance_charge_fluctuation
                 cpu/observe/test_capacitance_charge_fluctuation.cpp)
atv_add_cpp_test(test_cpu_observer_potential_profile
                 cpu/observe/test_potential_profile.cpp)
atv_add_cpp_test(test_cpu_observer_species_density_profile
                 cpu/observe/test_species_density_profile.cpp)

atv_add_cpp_test(test_cpu_md_no_alloc              cpu/test_md_no_alloc.cpp)
atv_add_cpp_test(test_cpu_driver_role_validation cpu/test_cpu_driver_role_validation.cpp)
# Registration-time PhaseKind gate + five-driver plan-shape regression
# (arch-refactor-01-bindings). Uses AM1Evaluator (unconditional, no RHF gate).
atv_add_cpp_test(test_cpu_driver_binding_kind             cpu/test_driver_binding_kind.cpp)
atv_add_cpp_test(test_cpu_evaluate_driver_capabilities    cpu/test_evaluate_driver_capabilities.cpp)
atv_add_cpp_test(test_cpu_optimize_driver_early_stop      cpu/test_optimize_driver_early_stop.cpp)
atv_add_cpp_test(test_cpu_dynamics_full_budget            cpu/test_dynamics_full_budget.cpp)
atv_add_cpp_test(test_cpu_install_infra_hook              cpu/test_install_infra_hook.cpp)
atv_add_cpp_test(test_cpu_md_set_api                      cpu/dynamics/test_md_set_api.cpp)
atv_add_cpp_test(test_cpu_md_plan_equivalence             cpu/dynamics/test_md_plan_equivalence.cpp)
# TimestepProvider/TimestepConsumer opt-in protocol + the MD build-time push
# (esmm-cpmd-02-timestep-protocol): dt is single-sourced from the integrator,
# conflicting/absent/invalid publications throw at build. Local stubs only —
# no GPW dependency, one no-op MD step to fire the component build() pass.
atv_add_cpp_test(test_cpu_md_timestep_protocol            cpu/dynamics/test_timestep_protocol.cpp)
# Car-Parrinello conserved-energy assembler (esmm-cpmd-03-observability):
# kConserved = kEnergy + kKinetic + the two electronic keys, end_of_step @ -90,
# throw when the electronic keys are absent. Pure meta arithmetic on a seeded
# one-atom frame — no GPW, no SCF, no MD step.
atv_add_cpp_test(test_cpu_md_cpmd_conserved               cpu/dynamics/test_cpmd_conserved.cpp)
atv_add_cpp_test(test_cpu_minimize_evaluate_set_api       cpu/test_minimize_evaluate_set_api.cpp)
atv_add_cpp_test(test_cpu_materialize                     cpu/test_materialize.cpp)
atv_add_cpp_test(test_cpu_molrs_contract                  cpu/test_molrs_contract.cpp)

# cpu/task infrastructure cluster (mpi-launcher-embedded-python-02-codec):
# TaskSpec data model + MPI-free byte codec + rank-agnostic materializer.
# Serial, no mpirun; the AM1 single point is the ms-scale H2 solve already
# established by test_cpu_am1_scaffold.
atv_add_cpp_test(test_cpu_task_codec                      cpu/test_task_codec.cpp)
atv_add_cpp_test(test_cpu_task_materializer               cpu/test_task_materializer.cpp)

# LauncherGate (mpi-launcher-embedded-python-04, ac-002): the PID-qualified
# ATV_LAUNCHER marker probe. Serial gtest — MPI is never initialized in this
# binary, which is exactly the state the degradation contract is about.
atv_add_cpp_test(test_cpu_launcher_gate                   cpu/task/test_launcher_gate.cpp)

# LauncherArgs (same link, ac-002/ac-006): header-only class in
# app/atomiverse/, which is NOT part of libatomiverse — the test compiles the
# header directly via this include dir.
atv_add_cpp_test(test_app_launcher_args                   app/test_launcher_args.cpp)
target_include_directories(test_app_launcher_args PRIVATE
    ${CMAKE_SOURCE_DIR}/app/atomiverse)

# InterpreterLocator (mpi-launcher-embedded-python-05, ac-004/ac-005): pure
# path-resolution + version-string logic — Python.h-free, no popen/exec, no
# interpreter ever executed. The TU lives in app/atomiverse/ (NOT part of
# libatomiverse), so this target compiles it directly as a second source.
atv_add_cpp_test(test_app_interpreter_locator
    app/test_interpreter_locator.cpp
    ${CMAKE_SOURCE_DIR}/app/atomiverse/interpreter_locator.cpp)
target_include_directories(test_app_interpreter_locator PRIVATE
    ${CMAKE_SOURCE_DIR}/app/atomiverse)

# ExtXYZ loading through molrs — inline fixtures only, no external file.
atv_add_cpp_test(test_cpu_io_extxyz                       cpu/io/test_cpu_io_extxyz.cpp)

# ── MPI tests (launched under `mpirun`) ─────────────────────────────────

if(ATV_USE_MPI AND MPIEXEC_EXECUTABLE)
    atv_add_mpi_exec(test_mpi_topology mpi/test_topology.cpp)
    atv_add_mpi_run(test_mpi_topology_n1 test_mpi_topology 1)
    atv_add_mpi_run(test_mpi_topology_n2 test_mpi_topology 2)
    atv_add_mpi_run(test_mpi_topology_n4 test_mpi_topology 4)
    atv_add_mpi_run(test_mpi_topology_n8 test_mpi_topology 8)

    atv_add_mpi_exec(test_mpi_halo_dual_image_ghosts mpi/test_halo_dual_image_ghosts.cpp)
    atv_add_mpi_run(test_mpi_halo_dual_image_ghosts_n1 test_mpi_halo_dual_image_ghosts 1)
    atv_add_mpi_run(test_mpi_halo_dual_image_ghosts_n2 test_mpi_halo_dual_image_ghosts 2)

    atv_add_mpi_exec(test_mpi_bcast_frame mpi/test_bcast_frame.cpp)
    atv_add_mpi_run(test_mpi_bcast_frame_n1 test_mpi_bcast_frame 1)
    atv_add_mpi_run(test_mpi_bcast_frame_n2 test_mpi_bcast_frame 2)
    atv_add_mpi_run(test_mpi_bcast_frame_n4 test_mpi_bcast_frame 4)

    atv_add_mpi_exec(test_mpi_halo_exchange mpi/test_halo_exchange.cpp)
    atv_add_mpi_run(test_mpi_halo_exchange_n1 test_mpi_halo_exchange 1)
    atv_add_mpi_run(test_mpi_halo_exchange_n2 test_mpi_halo_exchange 2)
    atv_add_mpi_run(test_mpi_halo_exchange_n8 test_mpi_halo_exchange 8)

    # Ghost-charge refresh after a mid-step Q rewrite (const-pot MPI DD).
    atv_add_mpi_exec(test_mpi_halo_refresh_charges mpi/test_halo_refresh_charges.cpp)
    atv_add_mpi_run(test_mpi_halo_refresh_charges_n1 test_mpi_halo_refresh_charges 1)
    atv_add_mpi_run(test_mpi_halo_refresh_charges_n2 test_mpi_halo_refresh_charges 2)

    atv_add_mpi_exec(test_mpi_lj_nve mpi/test_lj_nve.cpp)
    atv_add_mpi_run(test_mpi_lj_nve_n1 test_mpi_lj_nve 1)
    atv_add_mpi_run(test_mpi_lj_nve_n8 test_mpi_lj_nve 8)

    atv_add_mpi_exec(test_mpi_lj_energy_bruteforce
                     mpi/test_lj_energy_bruteforce.cpp)
    atv_add_mpi_run(test_mpi_lj_energy_bruteforce_n1
                    test_mpi_lj_energy_bruteforce 1)
    atv_add_mpi_run(test_mpi_lj_energy_bruteforce_n8
                    test_mpi_lj_energy_bruteforce 8)

    atv_add_mpi_exec(test_mpi_coul_energy_bruteforce
                     mpi/test_coul_energy_bruteforce.cpp)
    atv_add_mpi_run(test_mpi_coul_energy_bruteforce_n1
                    test_mpi_coul_energy_bruteforce 1)
    atv_add_mpi_run(test_mpi_coul_energy_bruteforce_n8
                    test_mpi_coul_energy_bruteforce 8)

    if(ATV_BUILD_KSPACE)
        atv_add_mpi_exec(test_mpi_pme_nacl_parity
                         mpi/test_pme_nacl_parity.cpp)
        atv_add_mpi_run(test_mpi_pme_nacl_parity_n1
                        test_mpi_pme_nacl_parity 1)
        atv_add_mpi_run(test_mpi_pme_nacl_parity_n2
                        test_mpi_pme_nacl_parity 2)
        atv_add_mpi_exec(test_mpi_pppm_brick_parity
                         mpi/test_pppm_brick_parity.cpp)
        atv_add_mpi_run(test_mpi_pppm_brick_parity_n1
                        test_mpi_pppm_brick_parity 1)
        atv_add_mpi_run(test_mpi_pppm_brick_parity_n2
                        test_mpi_pppm_brick_parity 2)
    endif()

    atv_add_mpi_exec(test_mpi_bonded_energy_bruteforce
                     mpi/test_bonded_energy_bruteforce.cpp)
    atv_add_mpi_run(test_mpi_bonded_energy_bruteforce_n1
                    test_mpi_bonded_energy_bruteforce 1)
    atv_add_mpi_run(test_mpi_bonded_energy_bruteforce_n8
                    test_mpi_bonded_energy_bruteforce 8)

    atv_add_mpi_exec(test_mpi_atom_migration mpi/test_atom_migration.cpp)
    atv_add_mpi_run(test_mpi_atom_migration_n1 test_mpi_atom_migration 1)
    atv_add_mpi_run(test_mpi_atom_migration_n2 test_mpi_atom_migration 2)

    atv_add_mpi_exec(test_mpi_minimize_lj mpi/test_minimize_lj.cpp)
    atv_add_mpi_run(test_mpi_minimize_lj_n1 test_mpi_minimize_lj 1)
    atv_add_mpi_run(test_mpi_minimize_lj_n2 test_mpi_minimize_lj 2)
    atv_add_mpi_run(test_mpi_minimize_lj_n4 test_mpi_minimize_lj 4)
    atv_add_mpi_run(test_mpi_minimize_lj_n8 test_mpi_minimize_lj 8)

    atv_add_mpi_exec(test_mpi_pe_chains_multirank mpi/test_pe_chains_multirank.cpp)
    atv_add_mpi_run(test_mpi_pe_chains_multirank_n2 test_mpi_pe_chains_multirank 2)

    atv_add_mpi_exec(test_mpi_pe_chains_8rank mpi/test_pe_chains_8rank.cpp)
    atv_add_mpi_run(test_mpi_pe_chains_8rank_n8 test_mpi_pe_chains_8rank 8)

    atv_add_mpi_exec(test_mpi_nose_hoover_t_ramp mpi/test_nose_hoover_t_ramp.cpp)
    atv_add_mpi_run(test_mpi_nose_hoover_t_ramp_n2 test_mpi_nose_hoover_t_ramp 2)

    atv_add_mpi_exec(test_mpi_pe_chain_crossing mpi/test_pe_chain_crossing.cpp)
    atv_add_mpi_run(test_mpi_pe_chain_crossing_n2 test_mpi_pe_chain_crossing 2)

    atv_add_mpi_exec(test_mpi_pe_full_multirank mpi/test_pe_full_multirank.cpp)
    atv_add_mpi_run(test_mpi_pe_full_multirank_n2 test_mpi_pe_full_multirank 2)

    atv_add_mpi_exec(test_mpi_pe_chain_drifting mpi/test_pe_chain_drifting.cpp)
    atv_add_mpi_run(test_mpi_pe_chain_drifting_n2 test_mpi_pe_chain_drifting 2)

    atv_add_mpi_exec(test_mpi_npt_multirank mpi/test_npt_multirank.cpp)
    atv_add_mpi_run(test_mpi_npt_multirank_n8 test_mpi_npt_multirank 8)

    atv_add_mpi_exec(test_mpi_nose_hoover_npt_multirank
                     mpi/test_nose_hoover_npt_multirank.cpp)
    atv_add_mpi_run(test_mpi_nose_hoover_npt_multirank_n8
                    test_mpi_nose_hoover_npt_multirank 8)

    atv_add_mpi_exec(test_mpi_deform_multirank mpi/test_deform_multirank.cpp)
    atv_add_mpi_run(test_mpi_deform_multirank_n8 test_mpi_deform_multirank 8)

    atv_add_mpi_exec(test_mpi_eam_multirank mpi/test_eam_multirank.cpp)
    atv_add_mpi_run(test_mpi_eam_multirank_n2 test_mpi_eam_multirank 2)
    atv_add_mpi_run(test_mpi_eam_multirank_n4 test_mpi_eam_multirank 4)

    atv_add_mpi_exec(test_mpi_rattle_multirank mpi/test_rattle_multirank.cpp)
    atv_add_mpi_run(test_mpi_rattle_multirank_n2 test_mpi_rattle_multirank 2)
    atv_add_mpi_run(test_mpi_rattle_multirank_n4 test_mpi_rattle_multirank 4)

    atv_add_mpi_exec(test_mpi_sw_multirank mpi/test_sw_multirank.cpp)
    atv_add_mpi_run(test_mpi_sw_multirank_n2 test_mpi_sw_multirank 2)
    atv_add_mpi_run(test_mpi_sw_multirank_n4 test_mpi_sw_multirank 4)

    atv_add_mpi_exec(test_mpi_pe_msu_pipeline mpi/test_pe_msu_pipeline.cpp)
    atv_add_mpi_run(test_mpi_pe_msu_pipeline_n8 test_mpi_pe_msu_pipeline 8)

    atv_add_mpi_exec(test_mpi_gather mpi/test_gather.cpp)
    atv_add_mpi_run(test_mpi_gather_n1 test_mpi_gather 1)
    atv_add_mpi_run(test_mpi_gather_n2 test_mpi_gather 2)
    atv_add_mpi_run(test_mpi_gather_n8 test_mpi_gather 8)

    atv_add_mpi_exec(test_mpi_rcb mpi/test_rcb.cpp)
    atv_add_mpi_run(test_mpi_rcb_n1 test_mpi_rcb 1)
    atv_add_mpi_run(test_mpi_rcb_n8 test_mpi_rcb 8)

    atv_add_mpi_exec(test_mpi_task_protocol mpi/test_task_protocol.cpp)
    atv_add_mpi_run(test_mpi_task_protocol_n1 test_mpi_task_protocol 1)
    atv_add_mpi_run(test_mpi_task_protocol_n2 test_mpi_task_protocol 2)
    atv_add_mpi_run(test_mpi_task_protocol_n4 test_mpi_task_protocol 4)

    # serial-on-replica seam (mpi-launcher-embedded-python-04): a driver with
    # set_serial_on_replica() must NOT decompose — every rank solves the full
    # replica (AM1 H2 single point) and lands on the identical energy.
    atv_add_mpi_exec(test_mpi_serial_on_replica mpi/test_serial_on_replica.cpp)
    atv_add_mpi_run(test_mpi_serial_on_replica_n1 test_mpi_serial_on_replica 1)
    atv_add_mpi_run(test_mpi_serial_on_replica_n2 test_mpi_serial_on_replica 2)
    atv_add_mpi_run(test_mpi_serial_on_replica_n4 test_mpi_serial_on_replica 4)
endif()

# ── PW infrastructure tests (require ATV_BUILD_PW) ─────────────────────

if(ATV_BUILD_PW)
    # Naming convention: test_{backend}_{domain}_{component}.
    atv_add_cpp_test(test_cpu_pw_g_vector   cpu/planewave/test_g_vector.cpp)
    atv_add_cpp_test(test_cpu_pw_upf        cpu/planewave/test_upf.cpp)
    atv_add_cpp_test(test_cpu_pw_kpoint     cpu/planewave/test_kpoint.cpp)
    atv_add_cpp_test(test_cpu_pw_crystal    cpu/planewave/test_crystal.cpp)
    atv_add_cpp_test(test_cpu_pw_ewald      cpu/planewave/test_ewald.cpp)
    atv_add_cpp_test(test_cpu_pw_components cpu/planewave/test_pw_components.cpp)
    atv_add_cpp_test(test_cpu_pw_si         cpu/planewave/test_pw_si.cpp)
    atv_add_cpp_test(test_cpu_pw_evaluator  cpu/planewave/test_pw_evaluator.cpp)
    atv_add_cpp_test(test_cpu_pw_method     cpu/planewave/test_pw_method.cpp)
    atv_add_cpp_test(test_cpu_pw_spin       cpu/planewave/test_pw_spin.cpp)
    atv_add_cpp_test(test_cpu_pw_forces     cpu/planewave/test_pw_forces.cpp)
    atv_add_cpp_test(test_cpu_pw_bomd       cpu/planewave/test_pw_bomd.cpp)
    atv_add_cpp_test(test_cpu_pw_stress     cpu/planewave/test_pw_stress.cpp)
    atv_add_cpp_test(test_cpu_pw_dfpt       cpu/planewave/test_dfpt.cpp)
        atv_add_cpp_test(test_cpu_gpw_forces       cpu/gpw/test_cpu_gpw_forces.cpp)
    atv_add_cpp_test(test_cpu_gpw_stress       cpu/gpw/test_cpu_gpw_stress.cpp)
    atv_add_cpp_test(test_cpu_gpw_uks          cpu/gpw/test_cpu_gpw_uks.cpp)
    atv_add_cpp_test(test_cpu_gpw_dftu         cpu/gpw/test_cpu_gpw_dftu.cpp)
    atv_add_cpp_test(test_cpu_pw_smearing   cpu/planewave/test_pw_smearing.cpp)
    atv_add_cpp_test(test_cpu_pw_bands_dos  cpu/planewave/test_pw_bands_dos.cpp)
    atv_add_cpp_test(test_cpu_pw_berry      cpu/planewave/test_berry_polarization.cpp)
    atv_add_cpp_test(test_cpu_pw_periodic_wannier cpu/planewave/test_periodic_wannier.cpp)

    # FFT test needs FFTW3 linkage
    add_executable(test_cpu_pw_fft cpu/planewave/test_fft.cpp)
    target_link_libraries(test_cpu_pw_fft PRIVATE atomiverse GTest::gtest_main fftw3_stub)
    set_target_properties(test_cpu_pw_fft PROPERTIES CXX_STANDARD 20)
    target_include_directories(test_cpu_pw_fft PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME test_cpu_pw_fft COMMAND test_cpu_pw_fft)
    set_tests_properties(test_cpu_pw_fft PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=${ATV_TEST_OMP_THREADS}")

    # MPI k-point parallel parity (launched under mpirun)
    if(ATV_USE_MPI AND MPIEXEC_EXECUTABLE)
        atv_add_mpi_exec(test_mpi_pw_kpoint mpi/test_pw_kpoint.cpp)
        atv_add_mpi_run(test_mpi_pw_kpoint_n1 test_mpi_pw_kpoint 1)
        atv_add_mpi_run(test_mpi_pw_kpoint_n2 test_mpi_pw_kpoint 2)
        atv_add_mpi_run(test_mpi_pw_kpoint_n4 test_mpi_pw_kpoint 4)
    endif()
    if(ATV_USE_MPI AND MPIEXEC_EXECUTABLE AND ATV_BUILD_GPW)
        atv_add_mpi_exec(test_mpi_gpw_kpoint mpi/test_gpw_kpoint.cpp)
        atv_add_mpi_run(test_mpi_gpw_kpoint_n1 test_mpi_gpw_kpoint 1)
        atv_add_mpi_run(test_mpi_gpw_kpoint_n2 test_mpi_gpw_kpoint 2)
        atv_add_mpi_run(test_mpi_gpw_kpoint_n4 test_mpi_gpw_kpoint 4)
    endif()
endif()

# ── Shared dense linear algebra (cpu/linalg) ──────────────────────────

# Canonical orthogonalization with S-eigenvalue thresholding
# (hf-robust-02). Header-only cpu/linalg/eigensolve.h — no RHF gate.
atv_add_cpp_test(test_cpu_linalg_canonical_orth
                 cpu/linalg/test_canonical_orthogonalization.cpp)

# Symmetric matrix power + nonsymmetric (dgeev) eigensolver
# (scf-review-01-linalg). Header-only cpu/linalg/eigensolve.h — no RHF gate.
atv_add_cpp_test(test_cpu_linalg_eigensolve
                 cpu/linalg/test_eigensolve.cpp)

# Preconditioned conjugate gradient (scf-review-01-linalg). Header-only
# cpu/linalg/linear_solve.h — no RHF gate.
atv_add_cpp_test(test_cpu_linalg_linear_solve
                 cpu/linalg/test_linear_solve.cpp)

# ── Electronic-state propagator tests ─────────────────────────────────
#
# esmm-cpmd-01-electronic-home: Cpmd/CpmdCplx (Car-Parrinello RATTLE),
# Cpmd/CpmdCplx and ASPCPredictor(Cplx) live in cpu/electronic/ —
# method-agnostic dense algebra over cpu/linalg, consumed by Gaussian SCF
# *and* GPW. Their sources compile under ATV_BUILD_HF (see the root
# CMakeLists HF source block), so these units gate on HF only — deliberately
# NOT on ATV_BUILD_GPW. Every test here is pure propagator algebra; the GPW
# build()/guess_P integration guards stay in the ATV_BUILD_GPW section above.
# XlBomd/XlBomdCplx compile in the always-on set (classical PIM XL dipoles
# consume the recursion in MD-only builds), so its unit is unconditional.

atv_add_cpp_test(test_cpu_electronic_xl_bomd cpu/electronic/test_xl_bomd.cpp)

if(ATV_BUILD_HF)
    atv_add_cpp_test(test_cpu_electronic_cpmd   cpu/electronic/test_cpmd.cpp)
    atv_add_cpp_test(test_cpu_electronic_aspc_predictor
                     cpu/electronic/test_aspc_predictor.cpp)
endif()

# ── RHF/SCF tests ─────────────────────────────────────────────────────

if(ATV_BUILD_HF)
    atv_add_cpp_test(test_cpu_scf_cart2sph     cpu/scf/test_cart2sph.cpp)
    atv_add_cpp_test(test_cpu_scf_multipole_ints cpu/scf/test_multipole_integrals.cpp)
    atv_add_cpp_test(test_cpu_scf_wannier        cpu/scf/test_wannier.cpp)
    atv_add_cpp_test(test_cpu_observer_multipole cpu/observe/test_multipole.cpp)
    # Löwdin populations + MK/CHELPG ESP fit — hand-built shells/determinant
    # and synthetic point-charge potentials only (hf-capability-gaps-04).
    atv_add_cpp_test(test_cpu_observer_lowdin cpu/observe/test_lowdin.cpp)
    atv_add_cpp_test(test_cpu_observer_esp_charges cpu/observe/test_esp_charges.cpp)
    atv_add_cpp_test(test_cpu_scf_spherical    cpu/scf/test_spherical_basis.cpp)
    atv_add_cpp_test(test_cpu_scf_ecp_gradient cpu/scf/test_ecp_gradient.cpp)
    # Static basis-table structure for the third-period + diffuse extension
    # (hf-capability-gaps-05-basis) — table pins vs PySCF, no SCF solve.
    atv_add_cpp_test(test_basis_tables         cpu/scf/test_basis_tables.cpp)
    # Semi-local def2-ECP: static table pins, EcpPotential lookup,
    # EcpIntegrals U-matrix vs external PySCF ECPscalar pins, IntegralPack
    # valence wiring, x2c double-counting gate (hf-capability-gaps-06-ecp).
    atv_add_cpp_test(test_ecp                  cpu/scf/test_ecp.cpp)
    # hf-capability-gaps-18-x2c: spin-free pVp kernel pins, X2cWorkspace
    # nonrel/exact-decoupling identities, ECP↔x2c mutual exclusion both
    # ways, forces() intercept. No heavy-atom SCF in this unit binary.
    atv_add_cpp_test(test_cpu_scf_x2c           cpu/scf/test_x2c.cpp)
    atv_add_cpp_test(test_cpu_scf_uhf_solver   cpu/scf/test_uhf_solver.cpp)
    # hf-capability-gaps-10-rohf: ROHF (Roothaan effective Fock, PySCF
    # canonical / Guest-Saunders convention) vs externally pinned PySCF.
    atv_add_cpp_test(test_cpu_scf_rohf_solver  cpu/scf/test_rohf_solver.cpp)
    # hf-capability-gaps-11-ghf: real two-component GHF — spinor
    # SlaterDeterminant invariants, H2 collapse to UHF, noncollinear win on
    # the frustrated H3 triangle vs externally pinned PySCF, collinear trap.
    atv_add_cpp_test(test_cpu_scf_ghf          cpu/scf/test_ghf.cpp)
    atv_add_cpp_test(test_cpu_scf_uks_solver   cpu/scf/test_uks_solver.cpp)
    atv_add_cpp_test(test_cpu_scf_xc_functional cpu/scf/test_xc_functional.cpp)
    atv_add_cpp_test(test_cpu_scf_evaluator     cpu/scf/test_scf_evaluator.cpp)
    # scf-driver-refactor-04b-wire: method-driver nesting + MeanFieldOutputs.
    atv_add_cpp_test(test_cpu_scf_method_driver_wire
                     cpu/scf/test_method_driver_wire.cpp)
    # esmm-cpmd-04-ao-cpmd: AO Car-Parrinello on the mean-field methods —
    # the refusal lattice, the TimestepConsumer seam, ScfGradient::
    # compute_at_density parity at a converged determinant, and the
    # zero-SCF-iteration production signature. H2/STO-3G only (2 AOs); the
    # E_cons drift physics is regressions/dynamics/cpmd-{rhf-h2,rks-water,
    # uhf-o2}.
    atv_add_cpp_test(test_cpu_scf_cpmd            cpu/scf/test_scf_cpmd.cpp)
    # BOMD+d3bj/d4 conservativeness: `MeanFieldMethod` adds E_disp to the total
    # energy, so `internal::ScfGradient` must add −∂E_disp/∂R to the force. That
    # reduces to whether D3BJ/D4 `forces()` differentiates their own `energy()`
    # — measured here by central difference, no SCF (ms). The end-to-end
    # conservation witness is regressions/scf/bomd-d3bj-ar2-nve.
    atv_add_cpp_test(test_cpu_scf_dispersion_gradient
                     cpu/scf/test_dispersion_gradient.cpp)
    # esmm-cpmd-05-scf-xlbomd: XL-BOMD (Niklasson JCP 130, 214109) as the
    # Gaussian-SCF guess — the xl_bomd/xl_bomd_k member surface, the
    # exclusivity gates (ASPC predictor, cpmd, Table-I order K in [3,9]),
    # auxiliary-history integrity (advance only on converged nested runs), and
    # the guess ladder explicit > xl > ASPC > last-P. H2O/STO-3G on a
    # ballistic nested-MD trajectory, ~1 s; the warm-step iteration margin and
    # the drift envelope are regressions/dynamics/xl-bomd-rhf-water.
    atv_add_cpp_test(test_cpu_scf_xl_bomd         cpu/scf/test_scf_xl_bomd.cpp)
    # esmm-contract-04-forces-requested: the ASPC-predictor-under-MC guard on
    # MeanFieldMethod::run_nested_erased — extrapolated density histories are
    # MD-only (rejected trials break the continuity the extrapolation assumes),
    # while plain last-P reuse under MC stays legal. H2/STO-3G, one sweep.
    atv_add_cpp_test(test_cpu_scf_mc_predictor_guard
                     cpu/scf/test_scf_mc_predictor_guard.cpp)
    # esmm-contract-05-accumulate-switch: the ONE entry that runs two SCF
    # epilogues — the ASPC corrector-miss fallback — must leave exactly one
    # contribution in the driver-owned accumulators (kEnergy /
    # kEnergyElectronic / kAtom_f*) now that the in-plan entry accumulates.
    # Four 6-step H2/STO-3G nested trajectories, ~20 ms. Serial for the same
    # ULP-determinism reason as test_cpu_scf_method_driver_parity: the pin is
    # bit-identity against the method's own energy() readback.
    atv_add_cpp_test(test_cpu_scf_aspc_fallback_accumulate
                     cpu/scf/test_scf_aspc_fallback_accumulate.cpp)
    set_tests_properties(test_cpu_scf_aspc_fallback_accumulate PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    # hf-capability-gaps-08-stability: MO orbital-Hessian A/B kernel +
    # Seeger-Pople stability analysis. NO SCF loop: the converged H2/STO-3G
    # determinant is symmetry-exact in closed form (one Fock build, no
    # iteration); the full RHF.run()+stability() path is
    # regressions/scf/h2_rhf_stability.
    atv_add_cpp_test(test_cpu_scf_orbital_hessian cpu/scf/test_orbital_hessian.cpp)
    atv_add_cpp_test(test_cpu_scf_stability       cpu/scf/test_stability.cpp)
    # hf-capability-gaps-14-tdhf: TDA + TDHF/RPA lowest roots over the same
    # A/B kernel (dense vs full-supermatrix vs Davidson parity, the
    # triplet-instability throw) + the RHF::tda()/tdhf() method surface.
    # Formaldehyde oracle e2e lives in regressions/scf/ch2o_tdhf.
    atv_add_cpp_test(test_cpu_scf_tdhf_pack       cpu/scf/test_tdhf_pack.cpp)
    atv_add_cpp_test(test_cpu_scf_tdhf_method     cpu/scf/test_tdhf_method.cpp)
    atv_add_cpp_test(test_cpu_scf_rks_tda         cpu/scf/test_rks_tda.cpp)
    atv_add_cpp_test(test_cpu_scf_tda_forces      cpu/scf/test_tda_forces.cpp)
    atv_add_cpp_test(test_cpu_scf_analytic_tda_forces cpu/scf/test_analytic_tda_forces.cpp)
    atv_add_cpp_test(test_cpu_scf_tda_nac         cpu/scf/test_tda_nac.cpp)
    atv_add_cpp_test(test_cpu_scf_nac_vectors        cpu/scf/test_nac_vectors.cpp)
    atv_add_cpp_test(test_cpu_scf_mp2_public      cpu/scf/test_mp2_public.cpp)
    atv_add_cpp_test(test_cpu_scf_ri_mp2          cpu/scf/test_ri_mp2.cpp)
    atv_add_cpp_test(test_cpu_scf_alda_diff       cpu/scf/test_alda_diff.cpp)
    atv_add_cpp_test(test_cpu_scf_full_grid_alda  cpu/scf/test_full_grid_alda.cpp)
    atv_add_cpp_test(test_cpu_scf_embedding       cpu/scf/test_embedding.cpp)
    # hf-capability-gaps-16-solvent: C-PCM cavity/Born oracle, f(ε), H2O
    # solvent vs gas, solvent_eps=0 vacuum identity, forces throw.
    atv_add_cpp_test(test_cpu_scf_pcm             cpu/scf/test_pcm.cpp)
    # hf-capability-gaps-17-bsse: ghost AtomList/ShellList/IntegralPack
    # workspace semantics + Boys–Bernardi Counterpoise 5-SCF identities
    # (He–He/STO-3G; locked bsse ≤ 0). Water-dimer e2e lives in
    # regressions/scf/h2o_dimer_cp_hf631gs.
    atv_add_cpp_test(test_ghost_basis             cpu/scf/test_ghost_basis.cpp)
    atv_add_cpp_test(test_counterpoise            cpu/scf/test_counterpoise.cpp)
    # hf-capability-gaps-13-cphf: second/skeleton derivative-integral FD
    # gates, B^x sum rules, CPHF iterative-vs-dense parity + the link-08
    # sigma_plus LHS identity (micro HeH+), and the H2/STO-3G analytic-
    # Hessian-vs-FD-of-analytic-forces first gate. The H2O/6-31G* e2e lives
    # in regressions/scf/h2o_rhf_freq.
    atv_add_cpp_test(test_cpu_scf_cphf            cpu/scf/test_cphf.cpp)
    atv_add_cpp_test(test_cpu_scf_analytic_hessian
                     cpu/scf/test_analytic_hessian.cpp)
    atv_add_cpp_test(test_cpu_scf_uhf_hessian
                     cpu/scf/test_uhf_hessian.cpp)
    # scf-review-04-meanfield-dedup: the FD nuclear-Hessian skeleton hoisted to
    # the protected CRTP template
    # MeanFieldMethod<Derived>::central_difference_nuclear_hessian, pinned in
    # isolation against an analytic quadratic mock force field (no SCF); and the
    # Ha/Bohr->eV/Angstrom gradient glue hoisted onto the tag ctor
    # GradientResult(FromHartreeBohr, ...), pinned against hand-computed
    # eV/Angstrom + eV literals (1e-12).
    atv_add_cpp_test(test_cpu_scf_fd_hessian_template
                     cpu/scf/test_fd_hessian_template.cpp)
    atv_add_cpp_test(test_cpu_scf_gradient_result
                     cpu/scf/test_gradient_result.cpp)
    # hf-capability-gaps-09-soscf: CIAH stepper units — exp(κ) unitarity,
    # gradient vs FD, matrix-free CIAH vs dense augmented eigensolve, trust
    # region, MixStage latch; plus the allowed micro H2/STO-3G default-off
    # invariance pin (milliseconds). The hard-case contrast is
    # regressions/scf/soscf_hard_case.
    atv_add_cpp_test(test_soscf                   cpu/scf/test_soscf.cpp)
    atv_add_cpp_test(test_cpu_scf_mp2           cpu/scf/test_mp2.cpp)
    # scf-review-02-eri-hoists: hosted spin-orbital ERI helpers
    # (SpinOrbital / antisym_eri / mp_denominator) in internal/spin_orbital_eri.h.
    atv_add_cpp_test(test_cpu_scf_spin_orbital_eri cpu/scf/test_spin_orbital_eri.cpp)
    atv_add_cpp_test(test_cpu_scf_fci           cpu/scf/test_fci.cpp)
    # scf-review-03-posthf-dedup: canonical determinant-CI evaluators
    # (choose/enumerate_strings/jw_single_sign/Slater-Condon elements/ci_element)
    # hoisted to internal/ci_engine.h; hand-computed 2-orbital pins (1e-12).
    atv_add_cpp_test(test_cpu_scf_ci_engine     cpu/scf/test_ci_engine.cpp)
    # scf-review-03-posthf-dedup: the single hoisted make_tau in
    # internal/eom_ccsd_common.h; hand-computed tau pins (1e-12, fac=1 and 0.5).
    atv_add_cpp_test(test_cpu_scf_eom_ccsd_common cpu/scf/test_eom_ccsd_common.cpp)
    atv_add_cpp_test(test_cpu_scf_koopmans      cpu/scf/test_koopmans.cpp)
    atv_add_cpp_test(test_cpu_scf_delta_scf     cpu/scf/test_delta_scf.cpp)
    atv_add_cpp_test(test_cpu_scf_integrals     cpu/scf/test_integrals.cpp)
    # scf-review-07-posthf-perf: block-assembled dense ERI (item 17) equivalence
    # guard — IntegralPack::eri_dense vs per-AO native::eri reference (≤1e-12) +
    # hard-coded anchors. Mirrors src/cpu/scf/two_electron.cpp.
    atv_add_cpp_test(test_cpu_scf_two_electron  cpu/scf/test_two_electron.cpp)
    atv_add_cpp_test(test_cpu_scf_rys_roots     cpu/scf/test_rys_roots.cpp)
    atv_add_cpp_test(test_cpu_scf_rys_eri       cpu/scf/test_rys_eri.cpp)
    # scf-review-06-scf-perf: fused one-sweep jk_rhf vs two-pass g_rhf+k
    # (≤1e-12 Ha F gate, d-shell fixture) + coulomb_exchange_fused thread
    # invariance. Mirrors cpu/scf/operators.cpp.
    atv_add_cpp_test(test_cpu_scf_operators     cpu/scf/test_operators.cpp)
    # hf-capability-gaps-12-df: def2-JKFIT table, g-shell 2c/3c integrals,
    # metric whitening, fitted J/K vs dense (single calls, no SCF).
    atv_add_cpp_test(test_density_fitting       cpu/scf/test_density_fitting.cpp)
    # hf-capability-gaps-21-cfmm: shell-pair distributions, near/far tables,
    # single H2O/STO-3G j(P) vs dense direct J (no converged SCF).
    atv_add_cpp_test(test_cfmm_j                cpu/scf/test_cfmm_j.cpp)
    atv_add_cpp_test(test_cpu_scf_primitives    cpu/scf/test_primitives.cpp)
    atv_add_cpp_test(test_cpu_scf_accelerators  cpu/scf/test_accelerators.cpp)
    # scf-review-08-ccsd-spatial-diis: CC amplitude DIIS accelerator
    # (internal/cc_diis.h) over flattened [t1;t2] vectors — Pulay B-matrix on a
    # contractive linear toy + single/near-singular/state edges (pure algebra).
    atv_add_cpp_test(test_cpu_scf_cc_diis       cpu/scf/test_cc_diis.cpp)
    # scf-review-08-ccsd-spatial-diis: closed-shell spatial-orbital CCSD kernel
    # (internal/ccsd_spatial.h) — H2/STO-3G E_corr pin (same PySCF pin as
    # regressions/scf/ccsd-h2-sto3g) + He/STO-3G zero-virtual edge.
    atv_add_cpp_test(test_cpu_scf_ccsd_spatial  cpu/scf/test_ccsd_spatial.cpp)
    # hf-capability-gaps-07-guess: Hückel/SAP guess matrices + EDIIS/ADIIS
    # far-field mixers (single-call math only; iteration-count contrast is
    # regressions/scf/h2o_scf_convergence).
    atv_add_cpp_test(test_cpu_scf_init_guess    cpu/scf/test_init_guess.cpp)
    atv_add_cpp_test(test_cpu_scf_ediis_adiis   cpu/scf/test_ediis_adiis.cpp)
    # Charge/multiplicity resolver + method run() validity guard
    # (hf-robust-01); tutorial ch02/ch03 regression anchor.
    atv_add_cpp_test(test_cpu_scf_spin_resolution cpu/scf/test_spin_resolution.cpp)
    # SlaterDeterminant value type: validating ctors, per-spin density,
    # expect/homo/lumo, with_occupation, validate(S)
    # (scf-slater-determinant-01-type).
    atv_add_cpp_test(test_cpu_scf_slater_determinant cpu/scf/test_slater_determinant.cpp)
    # Public SCF vocabulary types re-homed into dedicated headers: scf_enums.h,
    # gradient_result.h, scf_guess.h — without including scf_result.h
    # (scf-slater-determinant-02-rehome).
    atv_add_cpp_test(test_scf_type_homes cpu/scf/test_scf_type_homes.cpp)
    # SCF linear-dependence handling via canonical orthogonalization
    # (hf-robust-02); tutorial ch02 regression anchor.
    atv_add_cpp_test(test_cpu_scf_lindep         cpu/scf/test_scf_lindep.cpp)
    # Virtual-orbital level shift, PySCF get_fock ordering parity
    # (hf-robust-03); tutorial ch02 regression anchor.
    atv_add_cpp_test(test_cpu_scf_level_shift    cpu/scf/test_level_shift.cpp)
    atv_add_cpp_test(test_cpu_scf_solver        cpu/scf/test_scf_solver.cpp)
    # scf-driver-refactor-02-stages: five scf::Phase stage components wrapping
    # the internal SCF helpers, isolated bit-identical parity (unit tier).
    atv_add_cpp_test(test_cpu_scf_stages cpu/scf/stages/test_scf_stages.cpp)
    # Bit-identical parity compares two independent builds of the same kernels;
    # parallel reduction order (OpenMP / BLAS) is ULP-non-deterministic, so this
    # test pins serial execution. Do NOT loosen the assertions to tolerances —
    # the spec mandates bit-identity.
    set_tests_properties(test_cpu_scf_stages PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    # scf-driver-refactor-03a-convert: RHF/UHF/RKS/UKS as Driver<scf::Phase>,
    # regression-pinned against the 03a-parity energies/eigenvalues. Pins serial
    # execution for the same ULP-determinism reason as test_cpu_scf_stages.
    atv_add_cpp_test(test_cpu_scf_method_driver_parity
                     cpu/scf/test_method_driver_parity.cpp)
    set_tests_properties(test_cpu_scf_method_driver_parity PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    atv_add_cpp_test(test_cpu_scf_method_driver_no_alloc
                     cpu/scf/test_method_driver_no_alloc.cpp)
    atv_add_cpp_test(test_cpu_evaluate_driver        cpu/test_evaluate_driver.cpp)
    # esmm-contract-05-accumulate-switch (U7): two energy sources co-bound into
    # ONE cpu::Evaluate pass must both land in kEnergy, in either binding order.
    # HF-gated because the second source is a real cpu::scf::RHF (H2/STO-3G,
    # literature-pinned) — one SCF on 2 atoms, milliseconds. Serial for the same
    # ULP-determinism reason as test_cpu_scf_method_driver_parity.
    atv_add_cpp_test(test_cpu_evaluate_source_cobinding
                     cpu/evaluate/test_source_cobinding.cpp)
    set_tests_properties(test_cpu_evaluate_source_cobinding PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    atv_add_cpp_test(test_cpu_scf_minao         cpu/scf/test_minao.cpp)
    atv_add_cpp_test(test_cpu_scf_rys_eri_deriv cpu/scf/test_rys_eri_deriv.cpp)
    # test_cpu_minimizer_lbfgs relocated to benchmarks/qc/rhf/validation/
    # (val_min_lbfgs_geoopt) — full SCF-driven H2 geometry optimizations.
    atv_add_cpp_test(test_cpu_minimizer_line_search cpu/minimize/test_line_search.cpp)
    atv_add_cpp_test(test_cpu_minimizer_fire cpu/minimize/test_fire.cpp)
    atv_add_cpp_test(test_cpu_minimizer_cell_relax cpu/minimize/test_cell_relax.cpp)
    atv_add_cpp_test(test_cpu_minimizer_recorders cpu/minimize/test_trajectory_writer.cpp)
    atv_add_cpp_test(test_cpu_minimize_driver_swap_optimizer cpu/test_minimize_driver_swap_optimizer.cpp)
    atv_add_cpp_test(test_cpu_scf_orbital_info_observer cpu/scf/test_orbital_info_observer.cpp)
    # Vibrational analysis (finite-difference Hessian) is an end-to-end /
    # PySCF-parity suite — too slow for the unit gate. RELOCATED to
    # benchmarks/qc/vib/validation/ (run under -DATV_BUILD_BENCHMARKS=ON).
    atv_add_cpp_test(test_cpu_recorder_density_cube cpu/record/test_density_cube.cpp)
    # Mirrors cpu/scf/xc_functional.h — analytical hard-coded goldens (no live libxc).
    atv_add_cpp_test(test_xc_reference         cpu/scf/test_xc_reference.cpp)
    # scf-review-06-scf-perf: geometry-fixed XCGridBasis cache equivalence —
    # cached-basis vs rebuild-every-call evaluate_xc_batched are bitwise-equal
    # for LDA and PBE. Mirrors cpu/scf/xc_grid_engine.h.
    atv_add_cpp_test(test_xc_grid_engine       cpu/scf/test_xc_grid_engine.cpp)
    # libxc test oracle (sub-spec 02 of libxc-xc-fix creates the oracle/
    # subdirectory; this block lands ahead of it so the option flips on
    # cleanly the moment the directory exists).
    if(ATV_BUILD_XC_ORACLE)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cpu/scf/oracle/CMakeLists.txt")
            add_subdirectory(cpu/scf/oracle)
        else()
            message(STATUS
                "ATV_BUILD_XC_ORACLE=ON but tests/cpu/scf/oracle/ does "
                "not exist yet (sub-spec 02 of libxc-xc-fix will create it). "
                "Skipping oracle build.")
        endif()
    endif()
endif()

# ── CUDA backend tests ──────────────────────────────────────────────

if(ATV_BUILD_CUDA)
    atv_add_test(test_cuda_device_manager cuda/test_device_manager.cu)
    atv_add_test(test_cuda_domain_decomposition cuda/test_domain_decomposition.cu)
    atv_add_test(test_cuda_multi_gpu  cuda/dynamics/test_multi_gpu.cu)
    atv_add_test(test_cuda_frame      cuda/test_frame.cu)
    atv_add_test(test_cuda_materialize cuda/test_materialize.cu)
    atv_add_test(test_cuda_backend    cuda/test_backend.cu)
    atv_add_test(test_core_plan       core/test_plan.cu)
    atv_add_test(test_cuda_md         cuda/dynamics/test_md.cu)
    atv_add_test(test_cuda_md_thermostat cuda/dynamics/test_thermostat.cu)
    atv_add_test(test_cuda_md_zero_momentum cuda/dynamics/test_zero_momentum.cu)
    atv_add_test(test_cuda_md_components cuda/dynamics/test_cuda_md_components.cu)
    atv_add_test(test_cuda_neighborlist_skin_gating cuda/dynamics/test_neighborlist_skin_gating.cu)
    atv_add_test(test_cuda_cluster_contract cuda/dynamics/test_cuda_cluster_contract.cu)
    atv_add_test(test_cuda_cluster_mask cuda/dynamics/test_cuda_cluster_mask.cu)
    atv_add_test(test_pair_lj_tile_step0 cuda/potentials/test_pair_lj_tile_step0.cu)
    atv_add_test(test_cuda_atom_cell_sort cuda/dynamics/test_atom_cell_sort.cu)
    # Pulls in the private primitive header src/cuda/neighbor/atom_cell_sort.cuh
    target_include_directories(test_cuda_atom_cell_sort PRIVATE
        ${CMAKE_SOURCE_DIR}/src/cuda/neighbor)
    atv_add_test(test_cuda_pair_lj_no_d2h cuda/test_pair_lj_no_d2h.cu)
    target_compile_definitions(test_cuda_pair_lj_no_d2h PRIVATE
        ATV_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    atv_add_test(test_cuda_driver_no_per_phase_sync cuda/test_driver_no_per_phase_sync.cu)
    target_compile_definitions(test_cuda_driver_no_per_phase_sync PRIVATE
        ATV_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    atv_add_test(test_cuda_md_nve_limit cuda/dynamics/test_nve_limit.cu)
    atv_add_test(test_cuda_md_fix_langevin cuda/dynamics/test_fix_langevin.cu)
    atv_add_test(test_cuda_md_nose_hoover_npt cuda/dynamics/test_cuda_nose_hoover_npt.cu)
    # cuda::RDF on the molrs streaming accumulator (arch-refactor-07-molrs):
    # static-trajectory g(r) parity vs a host-fed cpu::compute::RDFAccumulator
    # + fail-fast ctor edges. Needs a GPU node (runs a cuda::md::MD).
    atv_add_test(test_cuda_observer_rdf_molrs cuda/observe/test_rdf_molrs.cu)
    atv_add_test(test_cuda_pair_lj_special_bonds cuda/potentials/test_pair_lj_special_bonds.cu)
    atv_add_test(test_cuda_potentials_bonded cuda/potentials/test_bonded.cu)
    # pair_kernel is cuda-owned (backend-mirror-02-absorb); its scalar-math test
    # moved here from the former cpu/potentials/test_pair_kernel.cpp.
    atv_add_test(test_cuda_pair_kernel cuda/potentials/test_pair_kernel.cu)
    atv_add_test(test_cuda_ewald_coul cuda/longrange/test_ewald_coul.cu)
    if(ATV_BUILD_KSPACE)
        atv_add_test(test_cuda_pme_coul cuda/longrange/test_pme_coul.cu)
    endif()
    # Manual A/B (not ctest): OpenMM cluster vs LAMMPS Verlet nonbonded path.
    # Run: ATV_CUDA_NONBONDED=openmm|lammps ./cuda_nonbonded_style_bench [n_side] [steps]
    add_executable(cuda_nonbonded_style_bench
                   ${CMAKE_SOURCE_DIR}/tools/cuda_nonbonded_style_bench.cu)
    target_link_libraries(cuda_nonbonded_style_bench PRIVATE atomiverse)
    set_target_properties(cuda_nonbonded_style_bench PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    if(ATV_BUILD_PW)
        atv_add_test(test_cuda_pw_fft cuda/planewave/test_fft.cu)
        target_link_libraries(test_cuda_pw_fft PRIVATE CUDA::cufft)
        atv_add_test(test_cuda_pw_poisson cuda/planewave/test_poisson.cu)
        target_link_libraries(test_cuda_pw_poisson PRIVATE CUDA::cufft)
        atv_add_test(test_cuda_pw_grid_map cuda/planewave/test_grid_map.cu)
        target_link_libraries(test_cuda_pw_grid_map PRIVATE CUDA::cufft)
        atv_add_test(test_cuda_pw_xc_lda cuda/planewave/test_xc_lda.cu)
        target_link_libraries(test_cuda_pw_xc_lda PRIVATE CUDA::cufft)
        atv_add_test(test_cuda_pw_veff cuda/planewave/test_veff.cu)
        target_link_libraries(test_cuda_pw_veff PRIVATE CUDA::cufft)
    endif()
    atv_add_test(test_cuda_minimizer_sd cuda/minimize/test_sd.cu)
    atv_add_test(test_cuda_minimizer_cg cuda/minimize/test_cg.cu)
    atv_add_test(test_cuda_minimizer_lbfgs cuda/minimize/test_lbfgs.cu)
    atv_add_test(test_cuda_scf        cuda/scf/test_scf.cu)
    atv_add_test(test_cuda_scf_rhf    cuda/scf/test_rhf.cu)
    atv_add_test(test_cuda_scf_device_integrals cuda/scf/test_device_integrals.cu)
    if(ATV_BUILD_PW)
        atv_add_test(test_cuda_pw_si  cuda/planewave/test_pw_si.cu)
        # Parity needs a gapped system: uses the checked-in ONCV LDA Si
        # pseudopotential from the pw benchmark reference data.
        target_compile_definitions(test_cuda_pw_si PRIVATE
            ATV_PW_UPF_PATH="${CMAKE_SOURCE_DIR}/tests/data/pw/Si_lda_nc.upf")
    endif()
    atv_add_test(test_cuda_driver_role_validation cuda/test_cuda_driver_role_validation.cu)
    # Registration-time PhaseKind gate, CUDA mirror
    # (arch-refactor-01-bindings). add()-only — needs no CUDA device.
    atv_add_test(test_cuda_driver_binding_kind cuda/test_driver_binding_kind.cu)

    # AsyncRecordBase / AsyncBridge contract guards. The subdir registers
    # benchmark-independent tests unconditionally and gates the
    # polyethylene_tension-dependent ones on ATV_BUILD_BENCHMARKS internally.
    add_subdirectory(cuda/asynchronous)

endif()

# Grep gate: the frame-block registry is THE single per-atom field
# enumeration (frame-block-01, ac-001). set_capacity / copy_row /
# pack_exchange / unpack_exchange / max_exchange_doubles iterate the
# registry; a hand-written member name inside any of them fails the build.
add_test(
    NAME test_frame_field_list_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/cpu/frame/field_list_scrub.sh
            ${CMAKE_SOURCE_DIR})

# Grep gate: the direct `atoms.<field>[` / `atoms.<field>(` per-atom column
# access shape (frame-block-02, ac-001) must not remain in src/cpu/ or
# include/atomiverse/cpu/ -- every consumer resolves a typed
# `Frame::span(kAtom_*)` once at the launch boundary instead. Pure grep —
# runs in every build configuration.
add_test(
    NAME test_frame_direct_access_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/cpu/frame/direct_access_scrub.sh
            ${CMAKE_SOURCE_DIR})

# Grep gate: no hand-written static_cast<int> on a Phase enum. The seven workflow
# Phase enums all start at 0, so an int phase cannot be told apart across
# workflows; ATV::slot_of is the one sanctioned conversion and the typed entry
# points (slot_of / slot_count / ctx.at / run_stage / make_binding) are the only
# way in. Pure grep — runs in every build configuration.
add_test(
    NAME test_stage_typing_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/core/stage_typing_scrub.sh
            ${CMAKE_SOURCE_DIR})

# backend-mirror-01-core gate: the facade layer is gone; contracts live under
# core/ and nothing includes the old facade include path. Pure grep — every config.
add_test(
    NAME test_no_facade_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/core/no_facade_scrub.sh
            ${CMAKE_SOURCE_DIR})

# backend-mirror-02-absorb gate: no top-level shared-math module; reciprocal /
# pair kernels are owned per backend, and cuda uses the full word longrange
# (never kspace). Pure grep — every build configuration.
add_test(
    NAME test_no_toplevel_shared_math_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/no_toplevel_shared_math_scrub.sh
            ${CMAKE_SOURCE_DIR})

# backend-mirror-03-paths gate: backend subdirectories are full English words,
# mirrored across cpu/ and cuda/ — no abbreviated legacy dir names or include
# paths. Pure grep — every build configuration.
add_test(
    NAME test_path_words_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/path_words_scrub.sh
            ${CMAKE_SOURCE_DIR})

# backend-mirror-04-driver gate: the loop node is Driver (no backend-prefixed
# driver class), erased base is DriverNode, CUDA termination API is Criterion.
add_test(
    NAME test_driver_names_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/driver_names_scrub.sh
            ${CMAKE_SOURCE_DIR})

# backend-mirror-05-dissolve gate: no CPU Backend god class; the loop runtime is
# composed from communication / device / asynchronous sub-runtimes.
add_test(
    NAME test_no_backend_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/no_backend_scrub.sh
            ${CMAKE_SOURCE_DIR})

# backend-mirror-06-gates: cross-backend layout invariants not covered by the
# per-phase scrubs — no runtime/loop/shared dirs, cuda mirrors cpu, tests mirror
# source. Pure grep/find — every build configuration.
add_test(
    NAME test_layout_mirror_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/layout_mirror_scrub.sh
            ${CMAKE_SOURCE_DIR})

# backend-mirror-07-docs gate: CLAUDE.md / AGENTS.md Architecture prose is free
# of the pre-restructure vocabulary (legacy class + path names). Pure grep.
add_test(
    NAME test_docs_vocab_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/docs_vocab_scrub.sh
            ${CMAKE_SOURCE_DIR})

# regression-harness gate: a case must resolve its data through the compiled-in
# absolute ATV_REG_DATA_DIR (never by probing the CWD or reading the
# environment), must not silently swap in a different Hamiltonian when an input
# is missing, and must declare its OpenMP width at registration. Pure grep —
# every build configuration.
add_test(
    NAME test_regression_harness_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/regression_harness_scrub.sh
            ${CMAKE_SOURCE_DIR})

# cuda-resurrect gate: the CUDA backend must not reference the deleted host CPU
# API (catch-all Real / core/real.h, renamed unit constants, stale host
# cpu::Frame accessor shapes). Pure grep over the cuda trees — runs on the login
# node without CUDA. Device cuda::Frame .atoms/.box/.scalars pointers stay legal.
add_test(
    NAME test_cuda_no_stale_api_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/cuda_no_stale_api_scrub.sh
            ${CMAKE_SOURCE_DIR})
add_test(
    NAME test_cuda_openmm_md_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/cuda_openmm_md_scrub.sh
            ${CMAKE_SOURCE_DIR})

# Grep gate: domain objects are constructed via class constructors — never
# free-function factories (build_integrals / *_config / make_atom_list / …).
# Pure grep — runs in every build configuration. See AGENTS.md OOP rule.
# Product-module DAG + source gating (python-on-demand-install-01)
add_test(
    NAME test_module_dag_auto_enable
    COMMAND ${CMAKE_COMMAND}
            -DATV_SOURCE_DIR=${CMAKE_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/architecture/module_dag_auto_enable.cmake)
set_tests_properties(test_module_dag_auto_enable PROPERTIES
    LABELS "architecture;cmake"
    TIMEOUT 300)

add_test(
    NAME test_module_off_excludes_sources
    COMMAND ${CMAKE_COMMAND}
            -DATV_SOURCE_DIR=${CMAKE_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/architecture/module_off_excludes_sources.cmake)
set_tests_properties(test_module_off_excludes_sources PROPERTIES
    LABELS "architecture;cmake"
    TIMEOUT 300)

add_test(
    NAME test_no_factory_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/no_factory_scrub.sh
            ${CMAKE_SOURCE_DIR})

# Grep gate (scf-review-02-eri-hoists): the AO→MO transforms live in the
# private internal/ao2mo.h, not public mp2.h, and the spin-orbital
# antisymmetrization helpers (struct SO / so_eri / spat_ / spin_ / phys_ /
# antisym_) exist in exactly one hosted place. Pure grep — every build config.
add_test(
    NAME test_no_ao2mo_public_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/no_ao2mo_public_scrub.sh
            ${CMAKE_SOURCE_DIR})

# Grep gate (scf-review-03-posthf-dedup): the EOM-CCSD intermediates
# (make_tau / struct EomImds; no free build_imds) and the determinant-CI engine
# (choose / enumerate_strings / occ_list / jw_single_sign / Slater-Condon
# elements / ci_element / dense_ci_solve / make_rdm12 / struct CiSolution) are
# each DEFINED in exactly one internal engine file. Pure grep — every config.
add_test(
    NAME test_posthf_dedup_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/posthf_dedup_scrub.sh
            ${CMAKE_SOURCE_DIR})

# Grep gate (scf-review-04-meanfield-dedup): the FD nuclear-Hessian skeleton
# lives only in the protected CRTP template
# MeanFieldMethod<Derived>::central_difference_nuclear_hessian (no FD
# displacement-loop sentinels hess_to_ev_ang2 / f_to_ha_bohr / (2.0 * h_bohr)
# in rks/uks/rohf.cpp), and the factory-named gradient glue make_result /
# make_force_result appears nowhere under src/cpu/scf (it is the tag ctor
# GradientResult(FromHartreeBohr, ...)). Pure grep — every build config.
add_test(
    NAME test_scf_dedup_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/architecture/scf_dedup_scrub.sh
            ${CMAKE_SOURCE_DIR})

# Grep gate: compatibility-only CPU aliases/overloads stay deleted and the
# TaskSpec wire format is strict v3 with an explicit ComponentId mapping.
# Deliberately excludes cpu/frame/ and CUDA/Torch sources.
add_test(
    NAME test_cpu_compatibility_surface_scrub
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/cpu/compatibility_surface_scrub.sh
            ${CMAKE_SOURCE_DIR})

# Grep gate: cuda::RDF must be migrated off the deleted legacy host histogram
# helpers (arch-refactor-07-molrs) — rdf.cuh must consume
# cpu::compute::RDFAccumulator via trajectory_metrics.h; only rdf_write_curve
# (the one g(r) emitter) may remain. The forbidden-token list lives in the
# gate script. Pure grep — runs in every build configuration, no CUDA needed.
add_test(
    NAME test_cuda_rdf_molrs_migration
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/cuda/observe/rdf_molrs_migration.sh
            ${CMAKE_SOURCE_DIR})

# ── Torch backend tests ─────────────────────────────────────────────

if(ATV_BUILD_TORCH)
    atv_add_test(test_torch_frame          torch/test_frame.cu)
    atv_add_test(test_torch_materialize    torch/test_materialize.cu)
    atv_add_test(test_torch_transfer       torch/test_transfer.cu)
    atv_add_test(test_torch_lj             torch/test_lj.cu)
    atv_add_test(test_torch_harmonic_bond  torch/test_harmonic_bond.cu)
    atv_add_test(test_torch_neighbor_list  torch/test_neighbor_list.cu)
endif()

# Diagnostic helper (not registered as a test): fixed-V_eff band parity.
if(ATV_BUILD_CUDA AND ATV_BUILD_PW)
    add_executable(diag_band_parity cuda/planewave/diag_band_parity.cu)
    target_link_libraries(diag_band_parity PRIVATE atomiverse)
    set_target_properties(diag_band_parity PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

# Every end-to-end / parity check has been physically relocated to
# benchmarks/{qc,md,mc,pw}/**/validation/ (val_scf_*, val_vib_*, val_md_*,
# val_min_lbfgs_geoopt, val_md_pe_nvt_multirank, the MC physics gates), built
# with -DATV_BUILD_BENCHMARKS=ON and run at acceptance time. Nothing is
# disabled-but-compiled here any more: if it is registered above, it runs.
