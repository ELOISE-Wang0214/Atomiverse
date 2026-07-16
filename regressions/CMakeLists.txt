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

atv_add_regression(md lj_argon)
atv_add_regression(mc gcmc_lj)
atv_add_regression(mc npt_lj)
atv_add_regression(md eam_cu)
atv_add_regression(md sw_si)
atv_add_regression(md rigid_water)
atv_add_regression(path neb_muller_brown)
atv_add_regression(architecture harmonic_driver_matrix)
atv_add_regression(minimizer lj_dimer)
atv_add_regression(forcefield tip3p_two_water)
atv_add_regression(path remd_harmonic)
atv_add_regression(sampling umbrella_harmonic)
atv_add_regression(semiempirical am1_water_force)
atv_add_regression(semiempirical am1_water_minimize)
atv_add_regression(semiempirical am1_water_bomd)
if(ATV_BUILD_KSPACE)
    atv_add_regression(md npt_pme_nacl)
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
