/// @file cpu/scf/uks_solve.cpp
/// @brief CPU UKS solver — thin wrapper around unified scf_solve().

#include "atomiverse/cpu/scf/uks_solve.h"
#include "atomiverse/cpu/scf/scf_solve.h"

namespace ATV {
namespace cpu {
namespace scf {

UKSResult uks_solve(const IntegralPack& integrals,
                    int n_alpha, int n_beta,
                    XCType xc_type,
                    int n_radial, int n_angular,
                    double tol, int max_iter) {
    auto config = uks_config(integrals, n_alpha, n_beta, xc_type, tol);
    config.max_iter = max_iter;
    config.n_radial = n_radial;
    config.n_angular = n_angular;
    auto r = scf_solve(integrals, config);

    UKSResult result;
    result.converged = r.converged;
    result.iterations = r.iterations;
    result.energy = r.energy;
    result.electronic_energy = r.energy - integrals.V_nuc;
    result.one_electron_energy = r.one_electron_energy;
    result.coulomb_energy = r.coulomb_energy;
    result.xc_energy = r.xc_energy;
    result.residual = r.residual;
    result.n_alpha = r.n_alpha;
    result.n_beta = r.n_beta;
    result.integrated_n_alpha = r.integrated_n_alpha;
    result.integrated_n_beta = r.integrated_n_beta;
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
