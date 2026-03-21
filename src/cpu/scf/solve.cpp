/// @file cpu/scf/solve.cpp
/// @brief CPU SCF solver: xtensor-blas reference implementation.

#include "atomiverse/cpu/scf/solve.h"

#include <cmath>

#include <xtensor/xnorm.hpp>
#include <xtensor/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

namespace ATV {
namespace cpu {
namespace scf {

SCFResult solve(Int n, Int n_occ, Real alpha, Real tol, int max_iter) {
    // Build H0: deterministic diagonally-dominant symmetric matrix
    xt::xtensor<Real, 2> H0 = xt::zeros<Real>({n, n});
    for (Int i = 0; i < n; ++i) {
        H0(i, i) = Real(n + 1 - i);
        for (Int j = i + 1; j < n; ++j) {
            Real val = Real(1) / Real(1 + (j - i));
            H0(i, j) = val;
            H0(j, i) = val;
        }
    }

    // Initialize density P = 0
    xt::xtensor<Real, 2> P = xt::zeros<Real>({n, n});
    xt::xtensor<Real, 2> P_old = xt::zeros<Real>({n, n});

    SCFResult result{};

    for (int iter = 0; iter < max_iter; ++iter) {
        // H = H0 + alpha * P
        xt::xtensor<Real, 2> H = H0 + alpha * P;

        // Diagonalize H -> eigenvalues + eigenvectors
        auto [eig, C] = xt::linalg::eigh(H);

        // Save P_old, rebuild P = C_occ * C_occ^T
        P_old = P;
        auto C_occ = xt::view(C, xt::all(), xt::range(0, n_occ));
        P = xt::linalg::dot(C_occ, xt::transpose(C_occ));

        // Convergence: ||P - P_old||_F
        Real norm = xt::norm_l2(P - P_old)();

        // Energy = trace(H * P) = sum of element-wise product
        Real energy = xt::sum(H * P)();

        result.residual = norm;
        result.energy = energy;
        result.eigenvalues = eig;
        result.density = P;
        result.iterations = iter + 1;

        if (iter > 0 && norm < tol) {
            result.converged = true;
            break;
        }
    }

    return result;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
