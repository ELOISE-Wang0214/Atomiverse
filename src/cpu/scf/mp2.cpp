/// @file cpu/scf/mp2.cpp
/// @brief Implementation of ao2mo_ovov + closed-shell / RI MP2 (see mp2.h).

#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/internal/density_fitting.h"
#include "atomiverse/cpu/scf/internal/mp2_pack.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

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
    xt::xtensor<double, 4> q1r =
        xt::reshape_view(q1, {no, n, n, n}); // (i,nu,lam,sig)

    // Quarter 2: nu -> a (virtual).  bring nu to the front, contract, restore.
    xt::xtensor<double, 4> q1t =
        xt::transpose(q1r, {1, 0, 2, 3}); // (nu,i,lam,sig)
    xt::xtensor<double, 2> q1tm = xt::reshape_view(q1t, {n, no * n * n});
    xt::xtensor<double, 2> q2 = xt::linalg::dot(xt::transpose(Cv), q1tm);
    xt::xtensor<double, 4> q2r =
        xt::reshape_view(q2, {nv, no, n, n}); // (a,i,lam,sig)
    xt::xtensor<double, 4> q2t =
        xt::transpose(q2r, {1, 0, 2, 3}); // (i,a,lam,sig)

    // Quarter 3: lambda -> j (occupied).  bring lambda to the front.
    xt::xtensor<double, 4> q2l =
        xt::transpose(q2t, {2, 0, 1, 3}); // (lam,i,a,sig)
    xt::xtensor<double, 2> q2lm = xt::reshape_view(q2l, {n, no * nv * n});
    xt::xtensor<double, 2> q3 = xt::linalg::dot(xt::transpose(Co), q2lm);
    xt::xtensor<double, 4> q3r =
        xt::reshape_view(q3, {no, no, nv, n}); // (j,i,a,sig)
    xt::xtensor<double, 4> q3t =
        xt::transpose(q3r, {1, 2, 0, 3}); // (i,a,j,sig)

    // Quarter 4: sigma -> b (virtual).  sigma is the trailing index already.
    xt::xtensor<double, 2> q3m = xt::reshape_view(q3t, {no * nv * no, n});
    xt::xtensor<double, 2> q4 = xt::linalg::dot(q3m, Cv);
    xt::xtensor<double, 4> ovov =
        xt::reshape_view(q4, {no, nv, no, nv}); // (i,a,j,b)
    return ovov;
}

xt::xtensor<double, 4> ao2mo_full(const xt::xtensor<double, 4>& ERI,
                                  const xt::xtensor<double, 2>& C, int n_mo) {
    const std::size_t n = ERI.shape(0);
    const std::size_t m = static_cast<std::size_t>(n_mo);

    // MO coefficient block, materialized [n_ao x m].
    xt::xtensor<double, 2> Cm = xt::view(C, xt::all(), xt::range(0, n_mo));

    // Quarter 1: mu -> p.  (mu | nu lam sig) -> (p | nu lam sig)
    xt::xtensor<double, 2> eri2 = xt::reshape_view(ERI, {n, n * n * n});
    xt::xtensor<double, 2> q1 = xt::linalg::dot(xt::transpose(Cm), eri2);
    xt::xtensor<double, 4> q1r =
        xt::reshape_view(q1, {m, n, n, n}); // (p,nu,lam,sig)

    // Quarter 2: nu -> q.  bring nu to the front, contract, restore.
    xt::xtensor<double, 4> q1t =
        xt::transpose(q1r, {1, 0, 2, 3}); // (nu,p,lam,sig)
    xt::xtensor<double, 2> q1tm = xt::reshape_view(q1t, {n, m * n * n});
    xt::xtensor<double, 2> q2 = xt::linalg::dot(xt::transpose(Cm), q1tm);
    xt::xtensor<double, 4> q2r =
        xt::reshape_view(q2, {m, m, n, n}); // (q,p,lam,sig)
    xt::xtensor<double, 4> q2t =
        xt::transpose(q2r, {1, 0, 2, 3}); // (p,q,lam,sig)

    // Quarter 3: lambda -> r.  bring lambda to the front.
    xt::xtensor<double, 4> q2l =
        xt::transpose(q2t, {2, 0, 1, 3}); // (lam,p,q,sig)
    xt::xtensor<double, 2> q2lm = xt::reshape_view(q2l, {n, m * m * n});
    xt::xtensor<double, 2> q3 = xt::linalg::dot(xt::transpose(Cm), q2lm);
    xt::xtensor<double, 4> q3r =
        xt::reshape_view(q3, {m, m, m, n}); // (r,p,q,sig)
    xt::xtensor<double, 4> q3t =
        xt::transpose(q3r, {1, 2, 0, 3}); // (p,q,r,sig)

    // Quarter 4: sigma -> s.  sigma is the trailing index already.
    xt::xtensor<double, 2> q3m = xt::reshape_view(q3t, {m * m * m, n});
    xt::xtensor<double, 2> q4 = xt::linalg::dot(q3m, Cm);
    xt::xtensor<double, 4> full =
        xt::reshape_view(q4, {m, m, m, m}); // (p,q,r,s)
    return full;
}

xt::xtensor<double, 2> mo_hcore(const xt::xtensor<double, 2>& Hcore,
                                const xt::xtensor<double, 2>& C, int n_mo) {
    // MO coefficient block, materialized [n_ao x n_mo].
    xt::xtensor<double, 2> Cm = xt::view(C, xt::all(), xt::range(0, n_mo));
    xt::xtensor<double, 2> tmp = xt::linalg::dot(Hcore, Cm); // [n_ao x n_mo]
    return xt::linalg::dot(xt::transpose(Cm), tmp);          // [n_mo x n_mo]
}

MP2Result mp2_energy(const IntegralPack& ints, const SlaterDeterminant& det,
                     double e_scf_ha) {
    using Spin = SlaterDeterminant::Spin;
    // No-fallback guards -- throw, never return a sentinel. (Convergence of
    // the reference is guarded by the caller, not here.)
    if (ints.ERI.size() == 0) {
        throw std::invalid_argument(
            "mp2_energy: ERI is not materialized; call IntegralPack{..., "
            "compute_eri=true}.");
    }
    if (det.n_electrons(Spin::alpha) != det.n_electrons(Spin::beta)) {
        throw std::invalid_argument(
            "mp2_energy: closed-shell reference only (n_alpha must equal "
            "n_beta).");
    }

    // Row-major copy of the alpha MO coefficients (hoisted once).
    const xt::xtensor<double, 2> C = det.coefficients(Spin::alpha);
    const int n_occ = det.n_electrons(Spin::alpha);
    const int n_mo = static_cast<int>(C.shape(1));
    const int n_vir = n_mo - n_occ;

    // No virtuals (e.g. He/STO-3G): correlation is exactly zero. Return a finite
    // 0.0 rather than entering the loop and dividing by an empty denominator.
    if (n_vir == 0) {
        return {0.0, e_scf_ha};
    }

    const xt::xtensor<double, 4> ovov = ao2mo_ovov(ints.ERI, C, n_occ, n_mo);
    const xt::xtensor<double, 1>& eps = det.orbital_energies(Spin::alpha);

    double e_corr = 0.0;
    for (int i = 0; i < n_occ; ++i) {
        for (int j = 0; j < n_occ; ++j) {
            for (int a = 0; a < n_vir; ++a) {
                for (int b = 0; b < n_vir; ++b) {
                    const double iajb = ovov(i, a, j, b);
                    const double ibja = ovov(i, b, j, a);
                    const double denom =
                        eps(i) + eps(j) - eps(n_occ + a) - eps(n_occ + b);
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

    return {e_corr, e_scf_ha + e_corr, /*ri=*/false};
}

MP2Result ri_mp2_energy(const IntegralPack& ints, const SlaterDeterminant& det,
                        double e_scf_ha, AuxBasis aux) {
    using Spin = SlaterDeterminant::Spin;
    if (det.n_electrons(Spin::alpha) != det.n_electrons(Spin::beta))
        throw std::invalid_argument(
            "ri_mp2_energy: closed-shell reference only (n_alpha must equal "
            "n_beta).");
    if (ints.shells.empty())
        throw std::invalid_argument(
            "ri_mp2_energy: IntegralPack has no shells for density fitting.");

    const RealMatrix& C_rm = det.coefficients(Spin::alpha);
    // Working AO coefficients may be spherical; DensityFitting lives in
    // Cartesian engine AO. Bridge with pack.c2s when present.
    const int n_occ = det.n_electrons(Spin::alpha);
    const int n_mo = static_cast<int>(C_rm.shape()[1]);
    const int n_vir = n_mo - n_occ;
    if (n_vir == 0)
        return {0.0, e_scf_ha, /*ri=*/true};

    internal::DensityFitting df(ints.shells, ints.atoms, aux);
    const int n_cart = df.n_ao();
    const int n_aux = df.n_aux();
    const RealMatrix& B = df.b_tensor(); // (n_cart² × n_aux), column-major

    RealMatrix C_cart;
    if (ints.spherical()) {
        C_cart = xt::linalg::dot(ints.c2s, C_rm); // [n_cart × n_mo]
    } else {
        if (static_cast<int>(C_rm.shape()[0]) != n_cart)
            throw std::invalid_argument(
                "ri_mp2_energy: AO dimension mismatch with DensityFitting");
        C_cart = C_rm;
    }

    // B_ia^P = Σ_μν C_μi C_νa B_μν^P
    // Layout: for each P, M_P[μ,ν] = B(μ*n_cart+ν, P) is column-major flat.
    const int nov = n_occ * n_vir;
    // Store B_ov[ia, P]
    xt::xtensor<double, 2> B_ov = xt::zeros<double>(
        {static_cast<std::size_t>(nov), static_cast<std::size_t>(n_aux)});

    for (int P = 0; P < n_aux; ++P) {
        // M[μ,ν] = B(μ*n + ν, P)
        RealMatrix M = xt::zeros<double>(
            {static_cast<std::size_t>(n_cart), static_cast<std::size_t>(n_cart)});
        for (int mu = 0; mu < n_cart; ++mu)
            for (int nu = 0; nu < n_cart; ++nu)
                M(mu, nu) =
                    B(static_cast<std::size_t>(mu * n_cart + nu),
                      static_cast<std::size_t>(P));
        // (i a) block: Co^T M Cv
        RealMatrix Co = xt::view(C_cart, xt::all(), xt::range(0, n_occ));
        RealMatrix Cv =
            xt::view(C_cart, xt::all(), xt::range(n_occ, n_mo));
        RealMatrix tmp = xt::linalg::dot(xt::transpose(Co), M);
        RealMatrix ia = xt::linalg::dot(tmp, Cv); // [n_occ × n_vir]
        for (int i = 0; i < n_occ; ++i)
            for (int a = 0; a < n_vir; ++a)
                B_ov(static_cast<std::size_t>(i * n_vir + a),
                     static_cast<std::size_t>(P)) = ia(i, a);
    }

    // (ia|jb) ≈ Σ_P B_ov(ia,P) B_ov(jb,P)
    const RealVector& eps = det.orbital_energies(Spin::alpha);
    double e_corr = 0.0;
    for (int i = 0; i < n_occ; ++i) {
        for (int j = 0; j < n_occ; ++j) {
            for (int a = 0; a < n_vir; ++a) {
                for (int b = 0; b < n_vir; ++b) {
                    const int ia = i * n_vir + a;
                    const int jb = j * n_vir + b;
                    const int ib = i * n_vir + b;
                    const int ja = j * n_vir + a;
                    double iajb = 0.0, ibja = 0.0;
                    for (int P = 0; P < n_aux; ++P) {
                        iajb += B_ov(ia, P) * B_ov(jb, P);
                        ibja += B_ov(ib, P) * B_ov(ja, P);
                    }
                    const double denom =
                        eps(i) + eps(j) - eps(n_occ + a) - eps(n_occ + b);
                    if (denom >= 0.0)
                        throw std::invalid_argument(
                            "ri_mp2_energy: non-negative orbital-energy "
                            "denominator; reference is degenerate or "
                            "non-aufbau.");
                    e_corr += iajb * (2.0 * iajb - ibja) / denom;
                }
            }
        }
    }
    return {e_corr, e_scf_ha + e_corr, /*ri=*/true};
}

} // namespace scf
} // namespace cpu
} // namespace ATV
