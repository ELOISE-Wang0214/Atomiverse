/// @file src/cpu/gpw/xl_bomd.cpp
/// @brief XlBomd implementation — Niklasson dissipative XL-BOMD recursion
///        (see header; coefficients: JCP 130, 214109 (2009), Table I).

#include "atomiverse/cpu/gpw/xl_bomd.h"

#include <stdexcept>
#include <string>

namespace ATV {
namespace cpu {
namespace gpw {

XlBomd::XlBomd(int k) : k_(k) {
    // Table I of Niklasson et al., J. Chem. Phys. 130, 214109 (2009):
    // dissipative (κ, α, c_0..c_K) rows for K = 3..9. Each row satisfies
    // Σ c_k = 0 and Σ k c_k = 0 (no bias on a constant or linearly drifting
    // history), pinned by unit test.
    switch (k) {
    case 3:
        kappa_ = 1.69;
        alpha_ = 0.150;
        c_ = {-2.0, 3.0, 0.0, -1.0};
        break;
    case 4:
        kappa_ = 1.75;
        alpha_ = 0.057;
        c_ = {-3.0, 6.0, -2.0, -2.0, 1.0};
        break;
    case 5:
        kappa_ = 1.82;
        alpha_ = 0.018;
        c_ = {-6.0, 14.0, -8.0, -3.0, 4.0, -1.0};
        break;
    case 6:
        kappa_ = 1.84;
        alpha_ = 0.0055;
        c_ = {-14.0, 36.0, -27.0, -2.0, 12.0, -6.0, 1.0};
        break;
    case 7:
        kappa_ = 1.86;
        alpha_ = 0.0016;
        c_ = {-36.0, 99.0, -88.0, 11.0, 32.0, -25.0, 8.0, -1.0};
        break;
    case 8:
        kappa_ = 1.88;
        alpha_ = 0.00044;
        c_ = {-99.0, 286.0, -286.0, 78.0, 78.0, -90.0, 42.0, -10.0, 1.0};
        break;
    case 9:
        kappa_ = 1.89;
        alpha_ = 0.00012;
        c_ = {-286.0, 858.0, -936.0, 364.0, 168.0, -300.0, 184.0, -63.0,
              12.0, -1.0};
        break;
    default:
        throw std::invalid_argument(
            "XlBomd: dissipation order K must be in [3, 9] (Table I of "
            "Niklasson et al., JCP 130, 214109 (2009)); got K="
            + std::to_string(k));
    }
}

const linalg::RealMatrix& XlBomd::aux() const {
    if (!ready_)
        throw std::logic_error(
            "XlBomd::aux: no auxiliary density yet — call advance() with the "
            "first SCF density before asking for a guess");
    return aux_;
}

void XlBomd::advance(const linalg::RealMatrix& d_scf) {
    if (d_scf.size() == 0)
        throw std::invalid_argument(
            "XlBomd::advance: empty SCF density (nothing to propagate)");
    if (!ready_) {
        // Stationary seed: a constant history at D is the recursion's fixed
        // point (Σ c_k = 0), so the auxiliary dynamics starts at rest on the
        // Born-Oppenheimer surface.
        hist_.assign(static_cast<std::size_t>(k_) + 1, d_scf);
        aux_ = d_scf;
        ready_ = true;
        return;
    }
    if (d_scf.shape(0) != hist_.front().shape(0)
        || d_scf.shape(1) != hist_.front().shape(1))
        throw std::invalid_argument(
            "XlBomd::advance: SCF density shape changed since seeding — "
            "reset() on a basis/geometry discontinuity");

    // P_{n+1} = 2 P_n − P_{n−1} + κ (D_n − P_n) + α Σ_k c_k P_{n−k}
    linalg::RealMatrix p_next =
        2.0 * hist_[0] - hist_[1] + kappa_ * (d_scf - hist_[0]);
    for (std::size_t k = 0; k < c_.size(); ++k)
        p_next += alpha_ * c_[k] * hist_[k];

    hist_.push_front(std::move(p_next));
    hist_.pop_back();
    aux_ = hist_.front();
}

void XlBomd::reset() {
    hist_.clear();
    aux_ = linalg::RealMatrix{};
    ready_ = false;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
