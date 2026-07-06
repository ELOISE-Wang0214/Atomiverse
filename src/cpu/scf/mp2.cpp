/// @file cpu/scf/mp2.cpp
/// @brief Implementation of ao2mo_ovov + closed-shell MP2 (see mp2.h).

#include "atomiverse/cpu/scf/mp2.h"

#include <cstddef>
#include <stdexcept>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/views/xview.hpp>

namespace ATV {
namespace cpu {
namespace scf {

xt::xtensor<double, 4> ao2mo_ovov(const xt::xtensor<double, 4>& ERI,
                                  const xt::xtensor<double, 2>& C, int n_occ,
                                  int n_mo) {
    const std::size_t n = ERI.shape(0);
    const std::size_t no = static_cast<std::size_t>(n_occ);
    const std::size_t nv = static_cast<std::size_t>(n_mo - n_occ);

    // Occupied / virtual MO coefficient blocks, materialized [n_ao x no|nv].
    xt::xtensor<double, 2> Co = xt::view(C, xt::all(), xt::range(0, n_occ));
    xt::xtensor<double, 2> Cv = xt::view(C, xt::all(), xt::range(n_occ, n_mo));

    // Quarter 1: mu -> i (occupied).  (mu | nu lambda sigma) -> (i | nu lambda sigma)
    xt::xtensor<double, 2> eri2 = xt::reshape_view(ERI, {n, n * n * n});
    xt::xtensor<double, 2> q1 = xt::linalg::dot(xt::transpose(Co), eri2);
    xt::xtensor<double, 4> q1r = xt::reshape_view(q1, {no, n, n, n}); // (i,nu,lam,sig)

    // Quarter 2: nu -> a (virtual).  bring nu to the front, contract, restore.
    xt::xtensor<double, 4> q1t = xt::transpose(q1r, {1, 0, 2, 3}); // (nu,i,lam,sig)
    xt::xtensor<double, 2> q1tm = xt::reshape_view(q1t, {n, no * n * n});
    xt::xtensor<double, 2> q2 = xt::linalg::dot(xt::transpose(Cv), q1tm);
    xt::xtensor<double, 4> q2r = xt::reshape_view(q2, {nv, no, n, n}); // (a,i,lam,sig)
    xt::xtensor<double, 4> q2t = xt::transpose(q2r, {1, 0, 2, 3});     // (i,a,lam,sig)

    // Quarter 3: lambda -> j (occupied).  bring lambda to the front.
    xt::xtensor<double, 4> q2l = xt::transpose(q2t, {2, 0, 1, 3}); // (lam,i,a,sig)
    xt::xtensor<double, 2> q2lm = xt::reshape_view(q2l, {n, no * nv * n});
    xt::xtensor<double, 2> q3 = xt::linalg::dot(xt::transpose(Co), q2lm);
    xt::xtensor<double, 4> q3r = xt::reshape_view(q3, {no, no, nv, n}); // (j,i,a,sig)
    xt::xtensor<double, 4> q3t = xt::transpose(q3r, {1, 2, 0, 3});      // (i,a,j,sig)

    // Quarter 4: sigma -> b (virtual).  sigma is the trailing index already.
    xt::xtensor<double, 2> q3m = xt::reshape_view(q3t, {no * nv * no, n});
    xt::xtensor<double, 2> q4 = xt::linalg::dot(q3m, Cv);
    xt::xtensor<double, 4> ovov = xt::reshape_view(q4, {no, nv, no, nv}); // (i,a,j,b)
    return ovov;
}

MP2Result mp2_energy(const IntegralPack& ints, const SCFResult& scf) {
    // No-fallback guards -- throw, never return a sentinel.
    if (ints.ERI.size() == 0) {
        throw std::invalid_argument(
            "mp2_energy: ERI is not materialized; call build_integrals(..., "
            "compute_eri=true).");
    }
    if (!scf.converged) {
        throw std::invalid_argument("mp2_energy: SCF result has not converged.");
    }
    if (scf.n_alpha != scf.n_beta) {
        throw std::invalid_argument(
            "mp2_energy: closed-shell reference only (n_alpha must equal "
            "n_beta).");
    }

    const int n_occ = scf.n_alpha;
    const int n_mo = static_cast<int>(scf.C_alpha.shape(1));
    const int n_vir = n_mo - n_occ;

    // No virtuals (e.g. He/STO-3G): correlation is exactly zero. Return a finite
    // 0.0 rather than entering the loop and dividing by an empty denominator.
    if (n_vir == 0) {
        return {0.0, scf.energy};
    }

    const xt::xtensor<double, 4> ovov =
        ao2mo_ovov(ints.ERI, scf.C_alpha, n_occ, n_mo);
    const xt::xtensor<double, 1>& eps = scf.eigenvalues_alpha;

    double e_corr = 0.0;
    for (int i = 0; i < n_occ; ++i) {
        for (int j = 0; j < n_occ; ++j) {
            for (int a = 0; a < n_vir; ++a) {
                for (int b = 0; b < n_vir; ++b) {
                    const double iajb = ovov(i, a, j, b);
                    const double ibja = ovov(i, b, j, a);
                    const double denom = eps(i) + eps(j) -
                                         eps(n_occ + a) - eps(n_occ + b);
                    // Closed-shell aufbau reference: every virtual lies above
                    // every occupied orbital, so the denominator is strictly
                    // negative. A non-negative value signals a degenerate or
                    // non-aufbau reference -- reject rather than emit garbage.
                    if (denom >= 0.0) {
                        throw std::invalid_argument(
                            "mp2_energy: non-negative orbital-energy "
                            "denominator; reference is degenerate or "
                            "non-aufbau.");
                    }
                    e_corr += iajb * (2.0 * iajb - ibja) / denom;
                }
            }
        }
    }

    return {e_corr, scf.energy + e_corr};
}

} // namespace scf
} // namespace cpu
} // namespace ATV
