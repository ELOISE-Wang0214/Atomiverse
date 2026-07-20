/// @file cpu/scf/scf_result.h
/// @brief Public SCF result bag — `SCFResult` only. Transitional: scheduled
///        for deletion by scf-slater-determinant-03-unbag (wavefunction moves
///        to `SlaterDeterminant`, scalars to Frame meta, loop state to the
///        Driver). Enums live in scf_enums.h, the initial guess in
///        scf_guess.h, the gradient result in gradient_result.h.
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/scf_enums.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Converged (or last) mean-field result.
struct SCFResult {
    bool converged = false;
    int iterations = 0;
    double energy = 0; ///< Total energy [Hartree]
    double residual = 0;

    double one_electron_energy = 0;
    double coulomb_energy = 0;
    double exchange_energy = 0;
    double xc_energy = 0;
    double nuclear_repulsion = 0;

    int n_alpha = 0;
    int n_beta = 0;
    SpinType spin = SpinType::
        Restricted; ///< Spin treatment used (restricted/unrestricted).
    double spin_contamination = 0;

    xt::xtensor<double, 1> eigenvalues_alpha;
    xt::xtensor<double, 1> eigenvalues_beta;
    xt::xtensor<double, 2> P_alpha;
    xt::xtensor<double, 2> P_beta;
    xt::xtensor<double, 2> C_alpha;
    xt::xtensor<double, 2> C_beta;

    double integrated_n_alpha = 0;
    double integrated_n_beta = 0;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
