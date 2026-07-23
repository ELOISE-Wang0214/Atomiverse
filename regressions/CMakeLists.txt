# ── regressions/ — out-of-library textbook end-to-end cases ──────────────
#
# See README.md for the rules. In short: each case is ONE flat main() that
# consumes the PUBLIC library only, checks a cited LITERATURE value, runs in
# <= 60 s, and exits non-zero on failure.
#
# These are deliberately NOT registered with add_test(): CLAUDE.md's rule is
# that tests/ is unit-only and "if it is registered, it runs". Regressions are
# end-to-end, so they get their own target instead:
#
#   cmake -B build -DATV_BUILD_REGRESSIONS=ON
#   cmake --build build --target regressions      # builds AND runs every case

add_custom_target(regressions
    COMMENT "Running all end-to-end textbook regression cases")

# atv_add_regression(<domain> <case>)
#   expects regressions/<domain>/<case>/reg_<case>.cpp
function(atv_add_regression domain case)
    set(target reg_${case})
    add_executable(${target} ${domain}/${case}/${target}.cpp)
    target_link_libraries(${target} PRIVATE atomiverse)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 17)

    add_custom_target(run_${target}
        COMMAND ${target}
        DEPENDS ${target}
        COMMENT "regression [${domain}/${case}]")
    add_dependencies(regressions run_${target})
endfunction()

# Layout guarantee (always on): out-of-library consumer of the core contract
# headers — backend-mirror-01-core moved these from facade/ into core/.
atv_add_regression(layout core_contract_includes)
atv_add_regression(layout backend_owned_longrange)

if(ATV_BUILD_MD)
    # Layout: full-word backend include paths (backend-mirror-03-paths).
    atv_add_regression(layout path_words)
    atv_add_regression(layout driver_names)
    atv_add_regression(layout no_backend)
    atv_add_regression(layout architecture_gates)
    atv_add_regression(md lj_argon)
    atv_add_regression(md pr_npt_lj)
    atv_add_regression(md nph_lj)
    atv_add_regression(md pe_bonded_nve)
    atv_add_regression(md eam_cu)
    atv_add_regression(md sw_si)
    atv_add_regression(md rigid_water)
    atv_add_regression(md dihedral_charmm)
    atv_add_regression(md improper_harmonic)
    atv_add_regression(md gaff_protocol_smoke)
endif()
if(ATV_BUILD_MC)
    atv_add_regression(mc gcmc_lj)
    atv_add_regression(mc npt_lj)
    atv_add_regression(mc widom_hs)
endif()
if(ATV_BUILD_PW)
    atv_add_regression(pw h_force_fd)
    atv_add_regression(pw h_force_analytic)
    atv_add_regression(pw h_stress_cubic)
    if(ATV_BUILD_MINIMIZER)
        atv_add_regression(pw h_cell_isotropic)
    endif()
endif()
if(ATV_BUILD_HF)
    atv_add_regression(scf h2o_rhf_sto3g)
    atv_add_regression(scf h2o_rhf_def2svp_sph)
    atv_add_regression(scf h2o_rhf_def2svp_grad)
    atv_add_regression(scf h2o_rhf_dipole)
    atv_add_regression(scf h2o_embedding_field)
    atv_add_regression(scf embed_h2o_grad)
    atv_add_regression(scf field_h2o_grad)
    atv_add_regression(scf h2o_rhf_631gs_cpcm)
    atv_add_regression(scf h2o_charges_pyscf)
    atv_add_regression(scf n2_ionization)
    atv_add_regression(scf diffuse_anion)
    atv_add_regression(scf ecp_heavy_atom)
    atv_add_regression(scf ecp_i2_grad)
    atv_add_regression(scf h2o_scf_convergence)
    atv_add_regression(scf h2_rhf_stability)
    atv_add_regression(scf soscf_hard_case)
    atv_add_regression(scf o2_rohf_sto3g)
    atv_add_regression(scf o2_rohf_grad)
    atv_add_regression(scf h2o_ri_jk)
    atv_add_regression(scf c10h22_cfmm_j)
    atv_add_regression(scf h3_ghf_sto3g)
    atv_add_regression(scf h2o_rhf_freq)
    atv_add_regression(scf h2o_rhf_ir)
    atv_add_regression(scf ch2o_tdhf)
    atv_add_regression(scf h2co_rks_tda)
    atv_add_regression(scf h2_rks_alda_tda)
    atv_add_regression(scf h2_full_grid_alda)
    atv_add_regression(scf h2_analytic_tda_force)
    atv_add_regression(scf h2_tda_force)
    atv_add_regression(scf h2_tda_nac)
    atv_add_regression(scf h2_nac_vectors)
    atv_add_regression(scf h2_alda_vs_cis)
    atv_add_regression(scf h2o_mp2_sto3g)
    atv_add_regression(scf h2o_ri_mp2)
    atv_add_regression(scf h2o_dimer_cp_hf631gs)
    atv_add_regression(scf x2c_argon)
endif()
if(ATV_BUILD_PATH)
    atv_add_regression(path neb_muller_brown)
    atv_add_regression(path remd_harmonic)
    atv_add_regression(path pimd_free_particle)
    atv_add_regression(path dimer_saddle)
    atv_add_regression(path string_muller_brown)
    atv_add_regression(path irc_muller_brown)
    atv_add_regression(path pimd_free_particle_estimators)
endif()
atv_add_regression(architecture harmonic_driver_matrix)
atv_add_regression(architecture module_compile_defs)
if(ATV_BUILD_MINIMIZER)
    atv_add_regression(minimizer lj_dimer)
    atv_add_regression(minimizer fire_lj_dimer)
    atv_add_regression(minimizer cell_lj_pressure)
    atv_add_regression(minimizer joint_lj_zero_p)
    atv_add_regression(minimizer aniso_lj_box)
    atv_add_regression(minimizer triclinic_lj_shear)
endif()
if(ATV_BUILD_POTENTIALS)
    atv_add_regression(forcefield tip3p_two_water)
    atv_add_regression(forcefield gaff_peo_fragment)
    atv_add_regression(forcefield lj_soft_dimer)
    atv_add_regression(forcefield coul_soft_dimer)
    atv_add_regression(forcefield yukawa_dimer)
endif()
if(ATV_BUILD_PATH)
    atv_add_regression(sampling umbrella_harmonic)
    atv_add_regression(sampling umbrella_angle)
    atv_add_regression(sampling umbrella_dihedral)
    atv_add_regression(sampling fep_soft_dimer)
    atv_add_regression(sampling metad_double_well)
    atv_add_regression(sampling metad_2d)
    atv_add_regression(sampling wham_double_well)
    atv_add_regression(path hremd_soft_lambda)
endif()
if(ATV_BUILD_MD)
    atv_add_regression(md brownian_einstein)
    atv_add_regression(md settle_tip3p)
    atv_add_regression(md fene_wca_chain)
    atv_add_regression(md dpd_groot_warren)
    atv_add_regression(md dpd_ab_mix)
endif()
if(ATV_BUILD_SEMIEMPIRICAL)
    atv_add_regression(semiempirical am1_water_force)
    atv_add_regression(semiempirical am1_water_minimize)
    atv_add_regression(semiempirical am1_water_bomd)
    atv_add_regression(semiempirical rm1_water_energy)
    atv_add_regression(semiempirical h2_dftb2)
    atv_add_regression(semiempirical h2o_dftb2_mio)
    atv_add_regression(semiempirical h2o_gfn0)
    atv_add_regression(multiscale qmmm_h2_ar)
    atv_add_regression(multiscale qmmm_es_h2)
    atv_add_regression(multiscale qmmm_link_ch3oh)
    atv_add_regression(multiscale qmmm_pol_h2o)
endif()
if(ATV_BUILD_KSPACE AND ATV_BUILD_MD)
    atv_add_regression(md npt_pme_nacl)
endif()
# FMM operator math is always-on (facade longrange/); no product-module gate.
atv_add_regression(longrange fmm_point_charges)

# ── pure-Python import (python-on-demand-install-02) ──────────────────────
# No C++ link: isolated pure sources only. Needs a Python interpreter.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_custom_target(run_reg_python_pure_import
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/python/pure_import/reg_pure_import.sh
                ${Python3_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/pure_import
        COMMENT "python pure_import: soft import without C++ extension"
        VERBATIM)
    add_dependencies(regressions run_reg_python_pure_import)

    add_custom_target(run_reg_python_install_cli
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/python/install_cli/reg_install_cli.sh
                ${Python3_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/install_cli
        COMMENT "python install_cli: help / mutual exclusion / hf / version lock"
        VERBATIM)
    add_dependencies(regressions run_reg_python_install_cli)
endif()

# ── launcher parity (mpi-launcher-embedded-python-04) ─────────────────────
# Not an atv_add_regression() case: it needs mpirun + a Python environment
# with atomiverse/molrs importable (the ~/work/.x86_64 venv per the spec).
# The launcher binary and interpreter come from the configured build.
if(ATV_BUILD_LAUNCHER)
    # regressions/ may be processed before app/atomiverse/ — resolve the
    # interpreter here too (cached, so this is a no-op when already found).
    find_package(Python3 COMPONENTS Interpreter QUIET)
    add_custom_target(run_reg_launcher_lj_parity
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/launcher/lj_parity/reg_lj_parity.sh
                ${Python3_EXECUTABLE}
                $<TARGET_FILE:atomiverse_launcher>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/launcher/lj_parity
        COMMENT "launcher lj_parity: serial-python vs 2-rank launcher (LJ + AM1)"
        VERBATIM)
    add_dependencies(run_reg_launcher_lj_parity atomiverse_launcher)
    add_dependencies(regressions run_reg_launcher_lj_parity)
endif()
