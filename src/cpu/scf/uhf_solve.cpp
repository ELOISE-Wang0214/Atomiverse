/// @file cpu/scf/uhf_solve.cpp
/// @brief CPU UHF solver — thin wrapper around unified scf_solve().

#include "atomiverse/cpu/scf/uhf_solve.h"
#include "atomiverse/cpu/scf/scf_solve.h"

namespace ATV {
namespace cpu {
namespace scf {

UHFResult uhf_solve(const IntegralPack& integrals, int n_alpha, int n_beta,
                    Real tol, int max_iter) {
    auto config = uhf_config(integrals, n_alpha, n_beta, tol);
    config.max_iter = max_iter;
    auto r = scf_solve(integrals, config);

    UHFResult result;
    result.converged = r.converged;
    result.iterations = r.iterations;
    result.energy = r.energy;
    result.electronic_energy = r.energy - integrals.V_nuc;
    result.residual = r.residual;
    result.n_alpha = r.n_alpha;
    result.n_beta = r.n_beta;
    result.spin_contamination = r.spin_contamination;
    result.eigenvalues_alpha = r.eigenvalues_alpha;
    result.eigenvalues_beta = r.eigenvalues_beta;
    result.P_alpha = r.P_alpha;
    result.P_beta = r.P_beta;
    result.C_alpha = r.C_alpha;
    result.C_beta = r.C_beta;
    result.trace = r.trace;
    return result;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
