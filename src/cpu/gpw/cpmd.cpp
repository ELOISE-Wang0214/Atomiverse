/// @file src/cpu/gpw/cpmd.cpp
/// @brief Cpmd implementation — Car–Parrinello electronic Verlet + Löwdin.

#include "atomiverse/cpu/gpw/cpmd.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xmath.hpp>

#include "atomiverse/cpu/gpw/internal/ot.h"

namespace ATV {
namespace cpu {
namespace gpw {

namespace {
/// Atomic unit of time in fs (ħ / E_h), CODATA-consistent with common AIMD.
constexpr double kAuTimeInFs = 0.024188843265857;

double frobenius_sq(const linalg::RealMatrix& a) {
    return xt::sum(a * a)();
}
} // namespace

Cpmd::Cpmd(double mu_au) : mu_(mu_au) {
    if (!(mu_ > 0.0) || !std::isfinite(mu_))
        throw std::invalid_argument(
            "Cpmd: fictitious mass mu must be finite and > 0 (got "
            + std::to_string(mu_au) + " a.u.)");
}

const linalg::RealMatrix& Cpmd::C() const {
    if (!ready_)
        throw std::logic_error(
            "Cpmd::C: no coefficients yet — call seed() after the cold SCF");
    return C_;
}

void Cpmd::seed(const linalg::RealMatrix& C_occ) {
    if (C_occ.size() == 0
        || static_cast<int>(C_occ.shape(0)) < 1
        || static_cast<int>(C_occ.shape(1)) < 1)
        throw std::invalid_argument(
            "Cpmd::seed: occupied coefficient matrix is empty");
    C_ = C_occ;
    C_prev_ = C_occ;
    ke_ha_ = 0.0;
    ready_ = true;
}

void Cpmd::align_to_overlap(const linalg::RealMatrix& S) {
    if (!ready_)
        throw std::logic_error(
            "Cpmd::align_to_overlap: not seeded");
    const int nb = static_cast<int>(C_.shape(0));
    if (static_cast<int>(S.shape(0)) != nb
        || static_cast<int>(S.shape(1)) != nb)
        throw std::invalid_argument(
            "Cpmd::align_to_overlap: S shape does not match C rows");
    C_ = internal::ot_orthonormalize(C_, S);
    C_prev_ = internal::ot_orthonormalize(C_prev_, S);
}

void Cpmd::advance(const linalg::RealMatrix& H, const linalg::RealMatrix& S,
                   double dt_fs) {
    if (!ready_)
        throw std::logic_error(
            "Cpmd::advance: not seeded — run a cold SCF and call seed() first");
    if (!(dt_fs > 0.0) || !std::isfinite(dt_fs))
        throw std::invalid_argument(
            "Cpmd::advance: dt_fs must be finite and > 0 (got "
            + std::to_string(dt_fs) + " fs)");

    const int nb = static_cast<int>(C_.shape(0));
    const int no = static_cast<int>(C_.shape(1));
    if (static_cast<int>(H.shape(0)) != nb
        || static_cast<int>(H.shape(1)) != nb)
        throw std::invalid_argument(
            "Cpmd::advance: H shape does not match C rows (nb="
            + std::to_string(nb) + ")");
    if (static_cast<int>(S.shape(0)) != nb
        || static_cast<int>(S.shape(1)) != nb)
        throw std::invalid_argument(
            "Cpmd::advance: S shape does not match C rows (nb="
            + std::to_string(nb) + ")");
    (void)no;

    // Projected electronic force (S-tangent): F = −(I − S C Cᵀ)(4 H C).
    // The unconstrained −4 H C still has a large occupied-space component
    // (Lagrange multipliers); at the SCF fixed point H C = S C ε that
    // component cancels only after the projection, so F → 0 on the BO surface.
    // Matches −ot_gradient (closed-shell P = 2 C Cᵀ).
    const linalg::RealMatrix F = -internal::ot_gradient(H, S, C_);
    const double dt_au = dt_fs / kAuTimeInFs;
    const double dt2_over_mu = (dt_au * dt_au) / mu_;

    // Position Verlet: C_raw = 2 C - C_prev + dt² F/μ
    linalg::RealMatrix C_raw = 2.0 * C_ - C_prev_ + dt2_over_mu * F;
    // Löwdin re-orthonormalization onto the current metric S.
    linalg::RealMatrix C_new = internal::ot_orthonormalize(C_raw, S);

    // Central-difference electronic velocity for T_e [Ha].
    const linalg::RealMatrix Cdot = (C_new - C_prev_) / (2.0 * dt_au);
    ke_ha_ = 0.5 * mu_ * frobenius_sq(Cdot);
    if (!std::isfinite(ke_ha_) || ke_ha_ < 0.0)
        ke_ha_ = 0.0;

    C_prev_ = std::move(C_);
    C_ = std::move(C_new);
}

void Cpmd::reset() {
    C_ = linalg::RealMatrix{};
    C_prev_ = linalg::RealMatrix{};
    ke_ha_ = 0.0;
    ready_ = false;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
