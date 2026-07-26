/// @file src/cpu/linalg/orbital_tangent.cpp
/// @brief Implementation of the S-metric orbital primitives (see header).

#include "atomiverse/cpu/linalg/orbital_tangent.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/core/xmath.hpp>

namespace ATV {
namespace cpu {
namespace linalg {

namespace {

using xt::linalg::dot;
using Cplx = std::complex<double>;

} // namespace

RealMatrix s_tangent(const RealMatrix& C, const RealMatrix& S,
                     const RealMatrix& d) {
    // (I - C C^T S) d  =>  C^T S (result) = 0  (uses C^T S C = I).
    return d - dot(C, dot(xt::transpose(C), dot(S, d)));
}

RealMatrix band_gradient(const RealMatrix& H, const RealMatrix& S,
                         const RealMatrix& C, double occ) {
    // Factor 2*occ: closed-shell occ=2 → 4 H C (P=2CCᵀ); UKS occ=1 → 2 H C.
    const RealMatrix W = (2.0 * occ) * dot(H, C);
    // (I - S C C^T) W
    return W - dot(S, dot(C, dot(xt::transpose(C), W)));
}

RealMatrix lowdin_orthonormalize(const RealMatrix& C, const RealMatrix& S) {
    const RealMatrix G = dot(xt::transpose(C), dot(S, C)); // no x no, SPD
    const RealMatrix Ginvhalf =
        sym_function(G, [](double t) { return 1.0 / std::sqrt(t); });
    return dot(C, Ginvhalf);
}

// ── Complex (Bloch) twins ──────────────────────────────────────────────────

// Row-major complex products (same shape as kpoints.cpp; no xtensor-blas).
CplxMatrix c_mul(const CplxMatrix& A, const CplxMatrix& B) {
    const int n = static_cast<int>(A.shape(0));
    const int m = static_cast<int>(A.shape(1));
    const int p = static_cast<int>(B.shape(1));
    CplxMatrix C = xt::zeros<Cplx>({n, p});
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
CplxMatrix c_dagger_mul(const CplxMatrix& A, const CplxMatrix& B) {
    const int n = static_cast<int>(A.shape(0));
    const int m = static_cast<int>(A.shape(1));
    const int p = static_cast<int>(B.shape(1));
    CplxMatrix C = xt::zeros<Cplx>({m, p});
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

CplxMatrix s_tangent_k(const CplxMatrix& C, const CplxMatrix& S,
                       const CplxMatrix& d) {
    // (I − C C† S) d
    return d - c_mul(C, c_dagger_mul(C, c_mul(S, d)));
}

CplxMatrix band_gradient_k(const CplxMatrix& H, const CplxMatrix& S,
                           const CplxMatrix& C, double occ) {
    const CplxMatrix W = (2.0 * occ) * c_mul(H, C);
    // (I − S C C†) W
    return W - c_mul(S, c_mul(C, c_dagger_mul(C, W)));
}

CplxMatrix band_gradient_k(const CplxMatrix& H, const CplxMatrix& S,
                           const CplxMatrix& C, const std::vector<double>& f) {
    const int nb = static_cast<int>(C.shape(0));
    const int no = static_cast<int>(C.shape(1));
    CplxMatrix HC = c_mul(H, C); // [nb x no]
    // W[:,i] = 2 f_i (H C)[:,i]  — fractional band occupations in [0, 2].
    for (int o = 0; o < no; ++o) {
        const double fo = (o < static_cast<int>(f.size()))
                              ? f[static_cast<std::size_t>(o)]
                              : 0.0;
        const double s = 2.0 * fo;
        for (int mu = 0; mu < nb; ++mu)
            HC(mu, o) *= s;
    }
    return HC - c_mul(S, c_mul(C, c_dagger_mul(C, HC)));
}

CplxMatrix lowdin_orthonormalize_k(const CplxMatrix& C, const CplxMatrix& S) {
    const CplxMatrix G = c_dagger_mul(C, c_mul(S, C));
    const CplxMatrix Ginvhalf = herm_function(G, [](double t) {
        if (!(t > 1e-16))
            throw std::runtime_error(
                "lowdin_orthonormalize_k: non-positive metric eigenvalue");
        return 1.0 / std::sqrt(t);
    });
    return c_mul(C, Ginvhalf);
}

} // namespace linalg
} // namespace cpu
} // namespace ATV
