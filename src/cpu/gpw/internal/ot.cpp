/// @file src/cpu/gpw/internal/ot.cpp
/// @brief Implementation of the OT inner minimizer (see header).

#include "atomiverse/cpu/gpw/internal/ot.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/views/xview.hpp>

namespace ATV {
namespace cpu {
namespace gpw {
namespace internal {

namespace {

using xt::linalg::dot;

/// Frobenius inner product Tr(A^T B) = sum_ij A_ij B_ij.
double frob(const RealMatrix& a, const RealMatrix& b) {
    return xt::sum(a * b)();
}

/// sin(u)/u with the u->0 Taylor limit (sinc(0) = 1).
double sinc(double u) {
    if (std::abs(u) < 1e-4) {
        const double u2 = u * u;
        return 1.0 - u2 / 6.0 + u2 * u2 / 120.0;
    }
    return std::sin(u) / u;
}

} // namespace

RealMatrix ot_parameterize(const RealMatrix& C0, const RealMatrix& x,
                           const RealMatrix& S) {
    const RealMatrix G = dot(xt::transpose(x), dot(S, x)); // no x no, SPD
    const RealMatrix cosU = linalg::sym_function(
        G, [](double t) { return std::cos(std::sqrt(std::max(t, 0.0))); });
    const RealMatrix sincU = linalg::sym_function(
        G, [](double t) { return sinc(std::sqrt(std::max(t, 0.0))); });
    return dot(C0, cosU) + dot(x, sincU);
}

OtPrecond::OtPrecond(const RealMatrix& T, const RealMatrix& S) {
    const RealMatrix M = T + S; // FULL_KINETIC, SPD, geometry-fixed
    m_inv_ = linalg::sym_function(M, [](double t) { return 1.0 / t; });
}

RealMatrix OtPrecond::apply(const RealMatrix& g) const {
    return dot(m_inv_, g);
}

OtPrecondSingleInverse::OtPrecondSingleInverse(const RealMatrix& H,
                                               const RealMatrix& S,
                                               const RealMatrix& C) {
    // Lowest occupied orbital energy from the small no x no Rayleigh block
    // eps = C^T H C (its eigenvalues are the occupied orbital energies at
    // self-consistency; the minimum is the most negative (H,S) eigenvalue for
    // a gapped system).
    const RealMatrix eps = dot(xt::transpose(C), dot(H, C));
    auto [theta, V] = xt::linalg::eigh(eps);
    double lambda_lo = theta(0);
    for (double t : theta)
        lambda_lo = std::min(lambda_lo, t);
    // Shift so M = H + a S has the occupied region at ~+margin in the S-metric
    // (positive-definite -> M^{-1} g is a descent direction). margin = 1 Ha.
    const double a = -lambda_lo + 1.0;
    const RealMatrix M = H + a * S;
    // M^{-1} via the symmetric eigendecomposition; guard tiny/negative ordinary
    // eigenvalues (PSD pseudo-inverse) so a stray mode can't flip the direction.
    m_inv_ = linalg::sym_function(
        M, [](double t) { return t > 1e-8 ? 1.0 / t : 0.0; });
}

RealMatrix OtPrecondSingleInverse::apply(const RealMatrix& g) const {
    return dot(m_inv_, g);
}

OtFixedResult ot_minimize_fixed(const RealMatrix& H, const RealMatrix& S,
                                const RealMatrix& T, const RealMatrix& C0,
                                double tol, int max_iter) {
    const OtPrecond precond(T, S);
    RealMatrix C = linalg::lowdin_orthonormalize(C0, S);

    auto band = [&](const RealMatrix& Cc) {
        return 2.0 * frob(Cc, dot(H, Cc)); // 2 Tr(C^T H C)
    };

    RealMatrix d, pz_prev;
    double denom_prev = 1.0;
    bool have_prev = false;

    OtFixedResult res;
    res.iters = 0;
    res.converged = false;

    for (int it = 0; it < max_iter; ++it) {
        res.iters = it + 1;
        const RealMatrix g = linalg::band_gradient(H, S, C);
        double gnorm = 0.0;
        for (double v : g)
            gnorm = std::max(gnorm, std::abs(v));
        if (gnorm < tol) {
            res.converged = true;
            break;
        }
        RealMatrix pz = linalg::s_tangent(
            C, S, precond.apply(g)); // preconditioned, tangent

        if (!have_prev) {
            d = -pz;
        } else {
            const double num = frob(g, pz - pz_prev); // Polak-Ribiere+
            const double beta = std::max(0.0, num / denom_prev);
            d = -pz + beta * d;
            d = linalg::s_tangent(C, S, d);
        }
        pz_prev = pz;
        denom_prev = frob(g, pz);
        have_prev = true;

        // Parabolic line search along the geodesic C(alpha d).
        const double E0 = band(C);
        const double slope = frob(g, d); // dE/dalpha at 0
        double alpha;
        if (slope >= 0.0) {
            d = -pz; // not a descent direction: reset to steepest descent
            alpha = 0.5;
        } else {
            const double at = 1.0;
            const RealMatrix Ct = ot_parameterize(C, at * d, S);
            const double Et = band(Ct);
            const double denom = 2.0 * (Et - E0 - slope * at);
            if (denom > 1e-300) {
                alpha = -slope * at * at / denom;
            } else {
                alpha = at; // non-positive curvature: take the trial step
            }
            // Backtrack if the parabola overshoots to a non-decreasing point.
            for (int bt = 0; bt < 30; ++bt) {
                const RealMatrix Ca = ot_parameterize(C, alpha * d, S);
                if (band(Ca) < E0)
                    break;
                alpha *= 0.5;
            }
        }
        C = linalg::lowdin_orthonormalize(ot_parameterize(C, alpha * d, S), S);
    }

    res.C = C;
    res.band_energy = band(C);
    return res;
}

// ── Complex (Bloch) OT ─────────────────────────────────────────────────────

namespace {

using Cplx = std::complex<double>;

// Re Tr(A† B) Frobenius for complex CG.
double c_frob(const linalg::CplxMatrix& a, const linalg::CplxMatrix& b) {
    double s = 0.0;
    const int n0 = static_cast<int>(a.shape(0));
    const int n1 = static_cast<int>(a.shape(1));
    for (int i = 0; i < n0; ++i)
        for (int j = 0; j < n1; ++j)
            s += (std::conj(a(i, j)) * b(i, j)).real();
    return s;
}

double c_sinc(double u) {
    if (std::abs(u) < 1e-4) {
        const double u2 = u * u;
        return 1.0 - u2 / 6.0 + u2 * u2 / 120.0;
    }
    return std::sin(u) / u;
}

} // namespace

linalg::CplxMatrix ot_parameterize_k(const linalg::CplxMatrix& C0,
                                     const linalg::CplxMatrix& x,
                                     const linalg::CplxMatrix& S) {
    // G = x† S x (no x no, Hermitian PSD).
    const linalg::CplxMatrix G = linalg::c_dagger_mul(x, linalg::c_mul(S, x));
    const linalg::CplxMatrix cosU = linalg::herm_function(
        G, [](double t) { return std::cos(std::sqrt(std::max(t, 0.0))); });
    const linalg::CplxMatrix sincU = linalg::herm_function(
        G, [](double t) { return c_sinc(std::sqrt(std::max(t, 0.0))); });
    return linalg::c_mul(C0, cosU) + linalg::c_mul(x, sincU);
}

OtPrecondSingleInverseK::OtPrecondSingleInverseK(const linalg::CplxMatrix& H,
                                                 const linalg::CplxMatrix& S,
                                                 const linalg::CplxMatrix& C) {
    const linalg::CplxMatrix eps = linalg::c_dagger_mul(C, linalg::c_mul(H, C));
    std::vector<double> theta;
    linalg::CplxMatrix V;
    linalg::eigh_hermitian(eps, theta, V);
    double lambda_lo = theta.empty() ? 0.0 : theta[0];
    for (double t : theta)
        lambda_lo = std::min(lambda_lo, t);
    const double a = -lambda_lo + 1.0;
    const linalg::CplxMatrix M = H + a * S;
    m_inv_ = linalg::herm_function(
        M, [](double t) { return t > 1e-8 ? 1.0 / t : 0.0; });
}

linalg::CplxMatrix
OtPrecondSingleInverseK::apply(const linalg::CplxMatrix& g) const {
    return linalg::c_mul(m_inv_, g);
}

OtFixedResultK ot_minimize_fixed_k(const linalg::CplxMatrix& H,
                                   const linalg::CplxMatrix& S,
                                   const linalg::CplxMatrix& T,
                                   const linalg::CplxMatrix& C0, double tol,
                                   int max_iter) {
    // FULL_KINETIC M = T + S.
    const linalg::CplxMatrix M = T + S;
    const linalg::CplxMatrix Minv = linalg::herm_function(
        M, [](double t) { return t > 1e-12 ? 1.0 / t : 0.0; });

    linalg::CplxMatrix C = linalg::lowdin_orthonormalize_k(C0, S);
    auto band = [&](const linalg::CplxMatrix& Cc) {
        return 2.0 * c_frob(Cc, linalg::c_mul(H, Cc)); // 2 Re Tr(C† H C)
    };

    linalg::CplxMatrix d, pz_prev;
    double denom_prev = 1.0;
    bool have_prev = false;

    OtFixedResultK res;
    for (int it = 0; it < max_iter; ++it) {
        res.iters = it + 1;
        const linalg::CplxMatrix g = linalg::band_gradient_k(H, S, C);
        double gnorm = 0.0;
        for (const Cplx v : g)
            gnorm = std::max(gnorm, std::abs(v));
        if (gnorm < tol) {
            res.converged = true;
            break;
        }
        linalg::CplxMatrix pz =
            linalg::s_tangent_k(C, S, linalg::c_mul(Minv, g));
        if (!have_prev) {
            d = -pz;
        } else {
            const double num = c_frob(g, pz - pz_prev);
            const double beta = std::max(0.0, num / denom_prev);
            d = linalg::s_tangent_k(C, S, -pz + beta * d);
        }
        pz_prev = pz;
        denom_prev = c_frob(g, pz);
        have_prev = true;

        const double E0 = band(C);
        double slope = c_frob(g, d);
        double alpha;
        if (slope >= 0.0) {
            d = -pz;
            slope = c_frob(g, d);
            alpha = 0.5;
        } else {
            const double at = 1.0;
            const double Et = band(ot_parameterize_k(C, at * d, S));
            const double denom = 2.0 * (Et - E0 - slope * at);
            alpha = (denom > 1e-300) ? -slope * at * at / denom : at;
            for (int bt = 0; bt < 30; ++bt) {
                if (band(ot_parameterize_k(C, alpha * d, S)) < E0)
                    break;
                alpha *= 0.5;
            }
        }
        C = linalg::lowdin_orthonormalize_k(ot_parameterize_k(C, alpha * d, S),
                                            S);
    }
    res.C = C;
    res.band_energy = band(C);
    return res;
}

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
