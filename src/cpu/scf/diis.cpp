/// @file cpu/scf/diis.cpp
/// @brief Combined DIIS accelerator implementation.

#include "atomiverse/cpu/scf/diis.h"

#include <xtensor/xnorm.hpp>

namespace ATV {
namespace cpu {
namespace scf {

std::pair<RealMatrix, RealMatrix> DIIS::extrapolate(
    const RealMatrix& F_alpha, const RealMatrix& e_alpha,
    const RealMatrix& F_beta,  const RealMatrix& e_beta) {

    // Store in history
    Fa_hist_.push_back(F_alpha);
    Fb_hist_.push_back(F_beta);
    ea_hist_.push_back(e_alpha);
    eb_hist_.push_back(e_beta);

    // Trim to max size
    while (static_cast<int>(Fa_hist_.size()) > max_vectors_) {
        Fa_hist_.pop_front(); Fb_hist_.pop_front();
        ea_hist_.pop_front(); eb_hist_.pop_front();
    }

    int m = static_cast<int>(Fa_hist_.size());
    if (m < 2) return {Fa_hist_.back(), Fb_hist_.back()};

    Int n = static_cast<Int>(F_alpha.shape(0));

    // Build combined B matrix: B_ij = Tr(e^α_i · e^α_j) + Tr(e^β_i · e^β_j)
    int dim = m + 1;
    RealMatrix B = xt::zeros<Real>({(Int)dim, (Int)dim});
    for (int i = 0; i < m; ++i)
        for (int j = 0; j <= i; ++j) {
            Real val = xt::sum(ea_hist_[i] * ea_hist_[j])()
                     + xt::sum(eb_hist_[i] * eb_hist_[j])();
            B(i, j) = val;
            B(j, i) = val;
        }
    // Lagrange constraint
    for (int i = 0; i < m; ++i) {
        B(m, i) = Real(-1);
        B(i, m) = Real(-1);
    }
    B(m, m) = Real(0);

    // RHS: [0, 0, ..., 0, -1]
    RealMatrix rhs = xt::zeros<Real>({(Int)dim, (Int)1});
    rhs(m, 0) = Real(-1);

    // Solve B * c = rhs
    try {
        auto c = xt::linalg::solve(B, rhs);
        RealMatrix Fa_diis = xt::zeros<Real>({n, n});
        RealMatrix Fb_diis = xt::zeros<Real>({n, n});
        for (int i = 0; i < m; ++i) {
            Fa_diis += Real(c(i, 0)) * Fa_hist_[i];
            Fb_diis += Real(c(i, 0)) * Fb_hist_[i];
        }
        return {Fa_diis, Fb_diis};
    } catch (...) {
        // DIIS linear solve failed — fall back to latest F
        return {Fa_hist_.back(), Fb_hist_.back()};
    }
}

RealMatrix DIIS::extrapolate(const RealMatrix& F, const RealMatrix& e) {
    auto [Fa, Fb] = extrapolate(F, e, F, e);
    return Fa;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
