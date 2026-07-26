/// @file src/cpu/gpw/internal/ot.cpp
/// @brief Implementation of the OT inner minimizer (see header).

#include "atomiverse/cpu/gpw/internal/ot.h"

#include "atomiverse/cpu/gpw/internal/kpoints.h" // eigh_hermitian

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
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

/// Symmetric matrix function f applied to an SPD/PSD matrix G = V diag(theta) V^T:
/// returns V diag(f(theta)) V^T. Used for cos(U), sinc(U), and M^{-1}, M^{-1/2}.
template <typename F>
RealMatrix sym_function(const RealMatrix& G, F f) {
    auto [theta, V] = xt::linalg::eigh(G);
    const std::size_t n = theta.size();
    RealMatrix D = xt::zeros<double>({n, n});
    for (std::size_t k = 0; k < n; ++k)
        D(k, k) = f(theta(k));
    return dot(dot(V, D), xt::transpose(V));
}

} // namespace

RealMatrix ot_parameterize(const RealMatrix& C0, const RealMatrix& x,
                           const RealMatrix& S) {
    const RealMatrix G = dot(xt::transpose(x), dot(S, x)); // no x no, SPD
    const RealMatrix cosU = sym_function(
        G, [](double t) { return std::cos(std::sqrt(std::max(t, 0.0))); });
    const RealMatrix sincU = sym_function(
        G, [](double t) { return sinc(std::sqrt(std::max(t, 0.0))); });
    return dot(C0, cosU) + dot(x, sincU);
}

RealMatrix ot_tangent(const RealMatrix& C0, const RealMatrix& S,
                      const RealMatrix& d) {
    // (I - C0 C0^T S) d  =>  C0^T S (result) = 0  (uses C0^T S C0 = I).
    return d - dot(C0, dot(xt::transpose(C0), dot(S, d)));
}

RealMatrix ot_gradient(const RealMatrix& H, const RealMatrix& S,
                       const RealMatrix& C0, double occ) {
    // Factor 2*occ: closed-shell occ=2 → 4 H C (P=2CCᵀ); UKS occ=1 → 2 H C.
    const RealMatrix W = (2.0 * occ) * dot(H, C0);
    // (I - S C0 C0^T) W
    return W - dot(S, dot(C0, dot(xt::transpose(C0), W)));
}

RealMatrix ot_orthonormalize(const RealMatrix& C, const RealMatrix& S) {
    const RealMatrix G = dot(xt::transpose(C), dot(S, C)); // no x no, SPD
    const RealMatrix Ginvhalf =
        sym_function(G, [](double t) { return 1.0 / std::sqrt(t); });
    return dot(C, Ginvhalf);
}

OtPrecond::OtPrecond(const RealMatrix& T, const RealMatrix& S) {
    const RealMatrix M = T + S; // FULL_KINETIC, SPD, geometry-fixed
    m_inv_ = sym_function(M, [](double t) { return 1.0 / t; });
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
    m_inv_ = sym_function(M, [](double t) { return t > 1e-8 ? 1.0 / t : 0.0; });
}

RealMatrix OtPrecondSingleInverse::apply(const RealMatrix& g) const {
    return dot(m_inv_, g);
}

OtFixedResult ot_minimize_fixed(const RealMatrix& H, const RealMatrix& S,
                                const RealMatrix& T, const RealMatrix& C0,
                                double tol, int max_iter) {
    const OtPrecond precond(T, S);
    RealMatrix C = ot_orthonormalize(C0, S);

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
        const RealMatrix g = ot_gradient(H, S, C);
        double gnorm = 0.0;
        for (double v : g)
            gnorm = std::max(gnorm, std::abs(v));
        if (gnorm < tol) {
            res.converged = true;
            break;
        }
        RealMatrix pz =
            ot_tangent(C, S, precond.apply(g)); // preconditioned, tangent

        if (!have_prev) {
            d = -pz;
        } else {
            const double num = frob(g, pz - pz_prev); // Polak-Ribiere+
            const double beta = std::max(0.0, num / denom_prev);
            d = -pz + beta * d;
            d = ot_tangent(C, S, d);
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
        C = ot_orthonormalize(ot_parameterize(C, alpha * d, S), S);
    }

    res.C = C;
    res.band_energy = band(C);
    return res;
}

// ── Complex (Bloch) OT ─────────────────────────────────────────────────────

namespace {

using Cplx = std::complex<double>;

// Row-major complex products (same shape as kpoints.cpp; no xtensor-blas).
OtCplxMatrix c_mul(const OtCplxMatrix& A, const OtCplxMatrix& B) {
    const int n = static_cast<int>(A.shape(0));
    const int m = static_cast<int>(A.shape(1));
    const int p = static_cast<int>(B.shape(1));
    OtCplxMatrix C = xt::zeros<Cplx>({n, p});
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < m; ++k) {
            const Cplx a = A(i, k);
            if (a == Cplx(0, 0))
                continue;
            for (int j = 0; j < p; ++j)
                C(i, j) += a * B(k, j);
        }
    return C;
}

// A† B
OtCplxMatrix c_dagger_mul(const OtCplxMatrix& A, const OtCplxMatrix& B) {
    const int n = static_cast<int>(A.shape(0));
    const int m = static_cast<int>(A.shape(1));
    const int p = static_cast<int>(B.shape(1));
    OtCplxMatrix C = xt::zeros<Cplx>({m, p});
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < n; ++k) {
            const Cplx a = std::conj(A(k, i));
            if (a == Cplx(0, 0))
                continue;
            for (int j = 0; j < p; ++j)
                C(i, j) += a * B(k, j);
        }
    return C;
}

// Re Tr(A† B) Frobenius for complex CG.
double c_frob(const OtCplxMatrix& a, const OtCplxMatrix& b) {
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

// Hermitian matrix function f(G) = V diag(f(θ)) V† via zheevd.
template <typename F>
OtCplxMatrix herm_function(const OtCplxMatrix& G, F f) {
    std::vector<double> theta;
    CplxMatrix V;
    // OtCplxMatrix and CplxMatrix are the same xtensor layout/type.
    eigh_hermitian(G, theta, V);
    const int n = static_cast<int>(theta.size());
    OtCplxMatrix out = xt::zeros<Cplx>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            Cplx s(0, 0);
            for (int k = 0; k < n; ++k) {
                const double fk = f(theta[static_cast<std::size_t>(k)]);
                s += V(i, k) * fk * std::conj(V(j, k));
            }
            out(i, j) = s;
        }
    return out;
}

} // namespace

OtCplxMatrix ot_parameterize_k(const OtCplxMatrix& C0, const OtCplxMatrix& x,
                               const OtCplxMatrix& S) {
    // G = x† S x (no x no, Hermitian PSD).
    const OtCplxMatrix G = c_dagger_mul(x, c_mul(S, x));
    const OtCplxMatrix cosU = herm_function(
        G, [](double t) { return std::cos(std::sqrt(std::max(t, 0.0))); });
    const OtCplxMatrix sincU = herm_function(
        G, [](double t) { return c_sinc(std::sqrt(std::max(t, 0.0))); });
    return c_mul(C0, cosU) + c_mul(x, sincU);
}

OtCplxMatrix ot_tangent_k(const OtCplxMatrix& C0, const OtCplxMatrix& S,
                          const OtCplxMatrix& d) {
    // (I − C0 C0† S) d
    return d - c_mul(C0, c_dagger_mul(C0, c_mul(S, d)));
}

OtCplxMatrix ot_gradient_k(const OtCplxMatrix& H, const OtCplxMatrix& S,
                           const OtCplxMatrix& C0, double occ) {
    const OtCplxMatrix W = (2.0 * occ) * c_mul(H, C0);
    // (I − S C0 C0†) W
    return W - c_mul(S, c_mul(C0, c_dagger_mul(C0, W)));
}

OtCplxMatrix ot_orthonormalize_k(const OtCplxMatrix& C, const OtCplxMatrix& S) {
    const OtCplxMatrix G = c_dagger_mul(C, c_mul(S, C));
    const OtCplxMatrix Ginvhalf = herm_function(G, [](double t) {
        if (!(t > 1e-16))
            throw std::runtime_error(
                "ot_orthonormalize_k: non-positive metric eigenvalue");
        return 1.0 / std::sqrt(t);
    });
    return c_mul(C, Ginvhalf);
}

OtPrecondSingleInverseK::OtPrecondSingleInverseK(const OtCplxMatrix& H,
                                                 const OtCplxMatrix& S,
                                                 const OtCplxMatrix& C) {
    const OtCplxMatrix eps = c_dagger_mul(C, c_mul(H, C));
    std::vector<double> theta;
    CplxMatrix V;
    eigh_hermitian(eps, theta, V);
    double lambda_lo = theta.empty() ? 0.0 : theta[0];
    for (double t : theta)
        lambda_lo = std::min(lambda_lo, t);
    const double a = -lambda_lo + 1.0;
    const OtCplxMatrix M = H + a * S;
    m_inv_ =
        herm_function(M, [](double t) { return t > 1e-8 ? 1.0 / t : 0.0; });
}

OtCplxMatrix OtPrecondSingleInverseK::apply(const OtCplxMatrix& g) const {
    return c_mul(m_inv_, g);
}

OtFixedResultK ot_minimize_fixed_k(const OtCplxMatrix& H, const OtCplxMatrix& S,
                                   const OtCplxMatrix& T, const OtCplxMatrix& C0,
                                   double tol, int max_iter) {
    // FULL_KINETIC M = T + S.
    const OtCplxMatrix M = T + S;
    const OtCplxMatrix Minv =
        herm_function(M, [](double t) { return t > 1e-12 ? 1.0 / t : 0.0; });

    OtCplxMatrix C = ot_orthonormalize_k(C0, S);
    auto band = [&](const OtCplxMatrix& Cc) {
        return 2.0 * c_frob(Cc, c_mul(H, Cc)); // 2 Re Tr(C† H C)
    };

    OtCplxMatrix d, pz_prev;
    double denom_prev = 1.0;
    bool have_prev = false;

    OtFixedResultK res;
    for (int it = 0; it < max_iter; ++it) {
        res.iters = it + 1;
        const OtCplxMatrix g = ot_gradient_k(H, S, C);
        double gnorm = 0.0;
        for (const Cplx v : g)
            gnorm = std::max(gnorm, std::abs(v));
        if (gnorm < tol) {
            res.converged = true;
            break;
        }
        OtCplxMatrix pz = ot_tangent_k(C, S, c_mul(Minv, g));
        if (!have_prev) {
            d = -pz;
        } else {
            const double num = c_frob(g, pz - pz_prev);
            const double beta = std::max(0.0, num / denom_prev);
            d = ot_tangent_k(C, S, -pz + beta * d);
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
        C = ot_orthonormalize_k(ot_parameterize_k(C, alpha * d, S), S);
    }
    res.C = C;
    res.band_energy = band(C);
    return res;
}

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
