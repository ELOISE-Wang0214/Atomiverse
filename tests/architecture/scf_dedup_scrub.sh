#!/usr/bin/env bash
# Grep gate (scf-review-04-meanfield-dedup): the two mean-field dedups landed —
# the FD nuclear-Hessian skeleton lives in exactly ONE place (the protected
# CRTP template MeanFieldMethod<Derived>::central_difference_nuclear_hessian in
# mean_field_method.cpp), and the three factory-named gradient-glue helpers are
# gone.
#
#   (a) rks.cpp / uks.cpp / rohf.cpp hessian() bodies no longer carry the FD
#       displacement-loop sentinels — the hoisted template owns the unit
#       constants and the central-difference denominator:
#         hess_to_ev_ang2    f_to_ha_bohr    (2.0 * h_bohr)
#   (b) NO make_result / make_force_result anywhere under src/cpu/scf — the
#       Ha/Bohr -> eV/Angstrom conversion + sign flip is now the tag ctor
#       GradientResult(FromHartreeBohr, grad_ha_bohr, energy_ha) (make_result
#       and make_force_result also violate the factory-function ban).
#
# Pure grep, seconds, every build config. Usage: scf_dedup_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: scf_dedup_scrub.sh <repo-root>}"

SCF_SRC="${ROOT}/src/cpu/scf"

# Drop comment-only match lines ("path:line:   // ...") so the gate tracks real
# code, mirroring the posthf_dedup_scrub.sh / no_factory_scrub.sh convention.
strip_comments() {
    grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*|/\*|#)' || true
}

status=0

# ── (a) no FD displacement-loop sentinels in rks/uks/rohf.cpp ──────────────
for f in rks uks rohf; do
    src="${SCF_SRC}/${f}.cpp"
    if [ ! -f "${src}" ]; then
        echo "FAIL: missing ${src}"
        status=1
        continue
    fi
    hits="$(grep -nE 'hess_to_ev_ang2|f_to_ha_bohr|\(2\.0 \* h_bohr\)' \
        "${src}" 2>/dev/null | strip_comments)"
    if [ -n "${hits}" ]; then
        echo "FAIL: FD displacement-loop sentinel still in ${f}.cpp — hoist the"
        echo "      skeleton into"
        echo "      MeanFieldMethod<Derived>::central_difference_nuclear_hessian:"
        echo "${hits}"
        status=1
    fi
done

# ── (b) no make_result / make_force_result under src/cpu/scf ───────────────
hits="$(grep -rnE '\<make_result\>|\<make_force_result\>' "${SCF_SRC}" \
    2>/dev/null | strip_comments)"
if [ -n "${hits}" ]; then
    echo "FAIL: factory-named gradient glue still present (must become the tag"
    echo "      ctor GradientResult(FromHartreeBohr, grad_ha_bohr, energy_ha)):"
    echo "${hits}"
    status=1
fi

if [ "${status}" -eq 0 ]; then
    echo "PASS: FD-Hessian skeleton hoisted once; make_result/make_force_result gone."
fi
exit "${status}"
