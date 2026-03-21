/// @file cpu/scf/rhf_solve.cpp
/// @brief CPU RHF solver — thin wrapper around unified scf_solve().

#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/scf_solve.h"

namespace ATV {
namespace cpu {
namespace scf {

RHFResult rhf_solve(const IntegralPack& integrals, Real tol, int max_iter) {
    auto config = rhf_config(integrals, tol);
    config.max_iter = max_iter;
    auto r = scf_solve(integrals, config);

    RHFResult result;
    result.converged = r.converged;
    result.iterations = r.iterations;
    result.energy = r.energy;
    result.electronic_energy = r.energy - integrals.V_nuc;
    result.residual = r.residual;
    result.eigenvalues = r.eigenvalues_alpha;
    result.P = r.P_alpha;
    result.C = r.C_alpha;
    result.trace = r.trace;
    return result;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
