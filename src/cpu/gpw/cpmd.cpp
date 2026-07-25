/// @file src/cpu/gpw/cpmd.cpp
/// @brief Cpmd implementation — Car–Parrinello electronic VV + RATTLE.

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
/// Position-RATTLE residual tolerance on ‖Cᵀ S C − I‖_∞.
constexpr double kRattlePosTol = 1e-12;
/// Max SHAKE iterations for the matrix-valued constraint.
constexpr int kRattlePosMaxIter = 40;

double frobenius_sq(const linalg::RealMatrix& a) { return xt::sum(a * a)(); }

double max_abs(const linalg::RealMatrix& a) {
    double m = 0.0;
    for (double v : a) {
        const double av = std::fabs(v);
        if (av > m)
            m = av;
    }
    return m;
}

/// Identity [n × n].
linalg::RealMatrix eye(int n) {
    linalg::RealMatrix I = xt::zeros<double>(
        {static_cast<std::size_t>(n), static_cast<std::size_t>(n)});
    for (int i = 0; i < n; ++i)
        I(i, i) = 1.0;
    return I;
}

/// Velocity RATTLE: kill sym(Cᵀ S Ċ) so Cᵀ S Ċ is antisymmetric (σ̇ = 0).
/// Ċ' = Ċ + C Γ with Γ = −sym(Cᵀ S Ċ), using Cᵀ S C = I.
linalg::RealMatrix rattle_velocity(const linalg::RealMatrix& Cdot,
                                   const linalg::RealMatrix& C,
                                   const linalg::RealMatrix& S) {
    using xt::linalg::dot;
    const linalg::RealMatrix W = dot(xt::transpose(C), dot(S, Cdot));
    const linalg::RealMatrix Gamma = -0.5 * (W + xt::transpose(W));
    return Cdot + dot(C, Gamma);
}

/// Position RATTLE (matrix SHAKE): C' = C_raw + C_ref Λ with C'ᵀ S C' = I.
/// Iterative linearization: M = C_refᵀ S C, residual = Cᵀ S C − I,
/// dΛ = −½ M⁻¹ residual (Andersen-style for near-identity M).
/// @return {C_constrained, iterations, residual_inf}
struct PosRattleResult {
    linalg::RealMatrix C;
    int iters;
    double residual_inf;
};

PosRattleResult rattle_position(const linalg::RealMatrix& C_raw,
                                const linalg::RealMatrix& C_ref,
                                const linalg::RealMatrix& S) {
    using xt::linalg::dot;
    const int no = static_cast<int>(C_raw.shape(1));
    const linalg::RealMatrix I = eye(no);

    linalg::RealMatrix C = C_raw;
    double res_inf = 0.0;
    int it = 0;
    for (; it < kRattlePosMaxIter; ++it) {
        const linalg::RealMatrix G = dot(xt::transpose(C), dot(S, C));
        const linalg::RealMatrix residual = G - I;
        res_inf = max_abs(residual);
        if (res_inf <= kRattlePosTol)
            return {std::move(C), it, res_inf};

        // M = C_refᵀ S C  (≈ I for small steps). Solve M dΛ + (M dΛ)ᵀ ≈ −res
        // by taking the symmetric particular solution Y = −½ res, dΛ = M⁻¹ Y.
        const linalg::RealMatrix M = dot(xt::transpose(C_ref), dot(S, C));
        linalg::RealMatrix Minv;
        try {
            Minv = xt::linalg::inv(M);
        } catch (...) {
            // Ill-conditioned M: fall back to Löwdin on the current C and stop.
            C = internal::ot_orthonormalize(C, S);
            const linalg::RealMatrix G2 = dot(xt::transpose(C), dot(S, C));
            res_inf = max_abs(G2 - I);
            return {std::move(C), it + 1, res_inf};
        }
        const linalg::RealMatrix dL = dot(Minv, -0.5 * residual);
        C = C + dot(C_ref, dL);
    }

    // Did not meet tol: final Löwdin clean-up so the propagator stays safe.
    C = internal::ot_orthonormalize(C, S);
    const linalg::RealMatrix Gf = dot(xt::transpose(C), dot(S, C));
    res_inf = max_abs(Gf - I);
    return {std::move(C), it, res_inf};
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

const linalg::RealMatrix& Cpmd::Cdot() const {
    if (!ready_)
        throw std::logic_error(
            "Cpmd::Cdot: no coefficients yet — call seed() after the cold SCF");
    return Cdot_;
}

void Cpmd::seed(const linalg::RealMatrix& C_occ) {
    if (C_occ.size() == 0 || static_cast<int>(C_occ.shape(0)) < 1
        || static_cast<int>(C_occ.shape(1)) < 1)
        throw std::invalid_argument(
            "Cpmd::seed: occupied coefficient matrix is empty");
    C_ = C_occ;
    Cdot_ = xt::zeros<double>(C_occ.shape());
    ke_ha_ = 0.0;
    last_pos_iters_ = 0;
    last_sigma_inf_ = 0.0;
    ready_ = true;
}

void Cpmd::align_to_overlap(const linalg::RealMatrix& S) {
    if (!ready_)
        throw std::logic_error("Cpmd::align_to_overlap: not seeded");
    const int nb = static_cast<int>(C_.shape(0));
    if (static_cast<int>(S.shape(0)) != nb
        || static_cast<int>(S.shape(1)) != nb)
        throw std::invalid_argument(
            "Cpmd::align_to_overlap: S shape does not match C rows");
    // Geometry change is not a dynamics step: Löwdin transport of C, then
    // project Ċ onto the new velocity-constraint manifold.
    C_ = internal::ot_orthonormalize(C_, S);
    Cdot_ = rattle_velocity(Cdot_, C_, S);
    using xt::linalg::dot;
    const int no = static_cast<int>(C_.shape(1));
    const linalg::RealMatrix G = dot(xt::transpose(C_), dot(S, C_));
    last_sigma_inf_ = max_abs(G - eye(no));
    last_pos_iters_ = 0;
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

    // Projected electronic force (S-tangent): F = −(I − S C Cᵀ)(4 H C).
    const linalg::RealMatrix F = -internal::ot_gradient(H, S, C_);
    const double dt_au = dt_fs / kAuTimeInFs;
    const double dt_over_2mu = dt_au / (2.0 * mu_);

    // 1) Velocity half-step (unconstrained).
    linalg::RealMatrix Cdot_half = Cdot_ + dt_over_2mu * F;
    // 2) Unconstrained position step.
    const linalg::RealMatrix C_raw = C_ + dt_au * Cdot_half;
    // 3) Position RATTLE: C' = C_raw + C_ref Λ, C'ᵀ S C' = I.
    auto pos = rattle_position(C_raw, /*C_ref=*/C_, S);
    last_pos_iters_ = pos.iters;
    last_sigma_inf_ = pos.residual_inf;
    if (!(last_sigma_inf_ <= 1e-8) || !std::isfinite(last_sigma_inf_))
        throw std::runtime_error(
            "Cpmd::advance: position RATTLE residual too large ("
            + std::to_string(last_sigma_inf_) + ")");

    // 4) Redefine half-step from constrained displacement (Andersen).
    Cdot_half = (pos.C - C_) / dt_au;
    // 5) Complete velocity with the same force (single H evaluation).
    linalg::RealMatrix Cdot_unc = Cdot_half + dt_over_2mu * F;
    // 6) Velocity RATTLE: sym(C'ᵀ S Ċ) = 0.
    Cdot_ = rattle_velocity(Cdot_unc, pos.C, S);
    C_ = std::move(pos.C);

    ke_ha_ = 0.5 * mu_ * frobenius_sq(Cdot_);
    if (!std::isfinite(ke_ha_) || ke_ha_ < 0.0)
        ke_ha_ = 0.0;
}

void Cpmd::reset() {
    C_ = linalg::RealMatrix{};
    Cdot_ = linalg::RealMatrix{};
    ke_ha_ = 0.0;
    last_pos_iters_ = 0;
    last_sigma_inf_ = 0.0;
    ready_ = false;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
