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
atv_add_regression(md eam_cu)
atv_add_regression(md sw_si)
atv_add_regression(md rigid_water)
atv_add_regression(path neb_muller_brown)
