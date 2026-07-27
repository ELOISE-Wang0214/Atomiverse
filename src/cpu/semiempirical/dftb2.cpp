/// @file src/cpu/semiempirical/dftb2.cpp
/// SCC-DFTB2 multi-element (H/C/N/O) with static mio-style SK tables.
#include "atomiverse/cpu/semiempirical/dftb2.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_mio_tables.h"

namespace ATV {
namespace cpu {
namespace semiempirical {
namespace {

using basis_data::DftbElementRow;
using basis_data::DftbPairRepRow;
using basis_data::DftbPairSkRow;
using basis_data::kDftbElements;
using basis_data::kDftbNElements;
using basis_data::kDftbNPairRep;
using basis_data::kDftbNPairSk;
using basis_data::kDftbPairRep;
using basis_data::kDftbPairSk;

const DftbElementRow& element_row(int z) {
    for (std::size_t i = 0; i < kDftbNElements; ++i)
        if (kDftbElements[i].z == z)
            return kDftbElements[i];
    throw std::runtime_error("DFTB2: unsupported element Z=" + std::to_string(z)
                             + " (static mio-style tables cover H/C/N/O only)");
}

const DftbPairSkRow& pair_sk(int za, int zb) {
    for (std::size_t i = 0; i < kDftbNPairSk; ++i)
        if (kDftbPairSk[i].z_a == za && kDftbPairSk[i].z_b == zb)
            return kDftbPairSk[i];
    throw std::runtime_error("DFTB2: no SK table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

const DftbPairRepRow& pair_rep(int za, int zb) {
    for (std::size_t i = 0; i < kDftbNPairRep; ++i)
        if (kDftbPairRep[i].z_a == za && kDftbPairRep[i].z_b == zb)
            return kDftbPairRep[i];
    throw std::runtime_error("DFTB2: no repulsive table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

double h_radial(double v0, double r, double r0, double d) {
    return v0 * std::exp(-(r - r0) / d);
}
double s_radial(double s0, double r, double r0, double w) {
    const double t = (r - r0) / w;
    return s0 * std::exp(-0.5 * t * t);
}

/// γ_AB(R) = 1 / sqrt(R^2 + (0.5(1/U_A + 1/U_B))^2)  (eV, Å; Ohno–Klopman-like)
double gamma_ab(double r, double ua, double ub) {
    const double tau = 0.5 * (1.0 / ua + 1.0 / ub);
    return 1.0 / std::sqrt(r * r + tau * tau);
}

/// dγ/dR = −R γ³.
double dgamma_dr(double r, double ua, double ub) {
    const double g = gamma_ab(r, ua, ub);
    return -r * g * g * g;
}

/// @brief Fill  n_a × n_b  block of H0/S for atoms A–B separated by R vec.
void fill_sk_block(const DftbElementRow& ea, const DftbElementRow& eb,
                   const DftbPairSkRow& sk, double dx, double dy, double dz,
                   double r, linalg::RealMatrix& H0, linalg::RealMatrix& S,
                   int ia0, int ib0, bool /*swap_sp*/) {
    if (r < 1e-12) {
        // Same center handled by caller (onsite).
        return;
    }
    const double inv_r = 1.0 / r;
    const double l = dx * inv_r;
    const double m = dy * inv_r;
    const double n = dz * inv_r;

    const double Hss = h_radial(sk.h_ss0, r, sk.h_r0, sk.h_d);
    const double Hsp = h_radial(sk.h_sp0, r, sk.h_r0, sk.h_d);
    const double Hpps = h_radial(sk.h_pps0, r, sk.h_r0, sk.h_d);
    const double Hppp = h_radial(sk.h_ppp0, r, sk.h_r0, sk.h_d);
    const double Sss = s_radial(sk.s_ss0, r, sk.s_r0, sk.s_w);
    const double Ssp = s_radial(sk.s_sp0, r, sk.s_r0, sk.s_w);
    const double Spps = s_radial(sk.s_pps0, r, sk.s_r0, sk.s_w);
    const double Sppp = s_radial(sk.s_ppp0, r, sk.s_r0, sk.s_w);

    auto set_hs = [&](int i, int j, double hv, double sv) {
        H0(ia0 + i, ib0 + j) = hv;
        H0(ib0 + j, ia0 + i) = hv;
        S(ia0 + i, ib0 + j) = sv;
        S(ib0 + j, ia0 + i) = sv;
    };

    // s–s
    set_hs(0, 0, Hss, Sss);

    // s(A)–p(B): H_s,pα = n_α * H_spσ  (direction from A to B)
    if (ea.n_ao >= 1 && eb.n_ao == 4) {
        set_hs(0, 1, l * Hsp, l * Ssp);
        set_hs(0, 2, m * Hsp, m * Ssp);
        set_hs(0, 3, n * Hsp, n * Ssp);
    }
    // p(A)–s(B): H_pα,s = -n_α * H_spσ  (antisym under direction flip)
    if (ea.n_ao == 4 && eb.n_ao >= 1) {
        set_hs(1, 0, -l * Hsp, -l * Ssp);
        set_hs(2, 0, -m * Hsp, -m * Ssp);
        set_hs(3, 0, -n * Hsp, -n * Ssp);
    }
    // p–p: H_αβ = n_α n_β (H_ppσ - H_ppπ) + δ_αβ H_ppπ
    if (ea.n_ao == 4 && eb.n_ao == 4) {
        const double dir[3] = {l, m, n};
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                const double d_ab = (a == b) ? 1.0 : 0.0;
                const double hv = dir[a] * dir[b] * (Hpps - Hppp) + d_ab * Hppp;
                const double sv = dir[a] * dir[b] * (Spps - Sppp) + d_ab * Sppp;
                set_hs(1 + a, 1 + b, hv, sv);
            }
        }
    }
}

/// SK block derivatives w.r.t. bond vector components (R_B − R_A).
void fill_sk_block_deriv(const DftbElementRow& ea, const DftbElementRow& eb,
                         const DftbPairSkRow& sk, double dx, double dy,
                         double dz, double r, linalg::RealMatrix dH0[3],
                         linalg::RealMatrix dS[3], int ia0, int ib0) {
    if (r < 1e-12)
        return;
    const double inv_r = 1.0 / r;
    const double dir[3] = {dx * inv_r, dy * inv_r, dz * inv_r};

    const double Hss = h_radial(sk.h_ss0, r, sk.h_r0, sk.h_d);
    const double Hsp = h_radial(sk.h_sp0, r, sk.h_r0, sk.h_d);
    const double Hpps = h_radial(sk.h_pps0, r, sk.h_r0, sk.h_d);
    const double Hppp = h_radial(sk.h_ppp0, r, sk.h_r0, sk.h_d);
    const double Sss = s_radial(sk.s_ss0, r, sk.s_r0, sk.s_w);
    const double Ssp = s_radial(sk.s_sp0, r, sk.s_r0, sk.s_w);
    const double Spps = s_radial(sk.s_pps0, r, sk.s_r0, sk.s_w);
    const double Sppp = s_radial(sk.s_ppp0, r, sk.s_r0, sk.s_w);

    const double dHss = -Hss / sk.h_d;
    const double dHsp = -Hsp / sk.h_d;
    const double dHpps = -Hpps / sk.h_d;
    const double dHppp = -Hppp / sk.h_d;
    const double dSss = Sss * (-(r - sk.s_r0) / (sk.s_w * sk.s_w));
    const double dSsp = Ssp * (-(r - sk.s_r0) / (sk.s_w * sk.s_w));
    const double dSpps = Spps * (-(r - sk.s_r0) / (sk.s_w * sk.s_w));
    const double dSppp = Sppp * (-(r - sk.s_r0) / (sk.s_w * sk.s_w));

    auto ddir = [&](int alpha, int k) {
        const double dab = (alpha == k) ? 1.0 : 0.0;
        return (dab - dir[alpha] * dir[k]) * inv_r;
    };

    for (int k = 0; k < 3; ++k) {
        const double dr_dk = dir[k];
        dH0[k](ia0, ib0) = dHss * dr_dk;
        dS[k](ia0, ib0) = dSss * dr_dk;

        if (ea.n_ao >= 1 && eb.n_ao == 4) {
            for (int a = 0; a < 3; ++a) {
                dH0[k](ia0, ib0 + 1 + a) =
                    ddir(a, k) * Hsp + dir[a] * dHsp * dr_dk;
                dS[k](ia0, ib0 + 1 + a) =
                    ddir(a, k) * Ssp + dir[a] * dSsp * dr_dk;
            }
        }
        if (ea.n_ao == 4 && eb.n_ao >= 1) {
            for (int a = 0; a < 3; ++a) {
                dH0[k](ia0 + 1 + a, ib0) =
                    -(ddir(a, k) * Hsp + dir[a] * dHsp * dr_dk);
                dS[k](ia0 + 1 + a, ib0) =
                    -(ddir(a, k) * Ssp + dir[a] * dSsp * dr_dk);
            }
        }
        if (ea.n_ao == 4 && eb.n_ao == 4) {
            const double dHsig = dHpps - dHppp;
            const double dSsig = dSpps - dSppp;
            const double Hsig = Hpps - Hppp;
            const double Ssig = Spps - Sppp;
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    const double dab = (a == b) ? 1.0 : 0.0;
                    dH0[k](ia0 + 1 + a, ib0 + 1 + b) =
                        (ddir(a, k) * dir[b] + dir[a] * ddir(b, k)) * Hsig
                        + dir[a] * dir[b] * dHsig * dr_dk + dab * dHppp * dr_dk;
                    dS[k](ia0 + 1 + a, ib0 + 1 + b) =
                        (ddir(a, k) * dir[b] + dir[a] * ddir(b, k)) * Ssig
                        + dir[a] * dir[b] * dSsig * dr_dk + dab * dSppp * dr_dk;
                }
            }
        }
    }
}

} // namespace

void DFTB2::build(const BuildContext& /*ctx*/) {}

void DFTB2::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n_at = f.nlocal();
    if (n_at < 1)
        throw std::invalid_argument("DFTB2: empty frame");

    const auto type = f.span(kAtom_type);
    const StorageReal* x = f.span(kAtom_x);
    const StorageReal* y = f.span(kAtom_y);
    const StorageReal* z = f.span(kAtom_z);
    AccReal* fx = f.span(kAtom_fx);
    AccReal* fy = f.span(kAtom_fy);
    AccReal* fz = f.span(kAtom_fz);
    double& energy_meta = f.meta.get(kEnergy);

    std::vector<const DftbElementRow*> el(static_cast<std::size_t>(n_at));
    std::vector<int> ao0(static_cast<std::size_t>(n_at) + 1, 0);
    double n_elec = 0;
    for (int a = 0; a < n_at; ++a) {
        el[static_cast<std::size_t>(a)] = &element_row(type[a]);
        n_elec += el[static_cast<std::size_t>(a)]->n_val;
        ao0[static_cast<std::size_t>(a) + 1] =
            ao0[static_cast<std::size_t>(a)]
            + el[static_cast<std::size_t>(a)]->n_ao;
    }
    const int n_ao = ao0[static_cast<std::size_t>(n_at)];
    if (n_ao < 1)
        throw std::runtime_error("DFTB2: no basis functions");
    // Closed-shell: even electron count required for n_occ = n_elec/2.
    if (std::fabs(n_elec - std::round(n_elec)) > 1e-9
        || static_cast<int>(std::lround(n_elec)) % 2 != 0)
        throw std::runtime_error(
            "DFTB2: odd/non-integer valence electron count not supported");
    const int n_occ = static_cast<int>(std::lround(n_elec)) / 2;

    auto dist_vec = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = static_cast<double>(x[j] - x[i]);
        dy = static_cast<double>(y[j] - y[i]);
        dz = static_cast<double>(z[j] - z[i]);
        f.box.min_image(dx, dy, dz);
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    linalg::RealMatrix H0 = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix S = xt::zeros<double>({n_ao, n_ao});

    // Onsite + pair SK
    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        H0(ia, ia) = ea.e_s;
        S(ia, ia) = 1.0;
        if (ea.n_ao == 4) {
            for (int p = 1; p < 4; ++p) {
                H0(ia + p, ia + p) = ea.e_p;
                S(ia + p, ia + p) = 1.0;
            }
        }
        for (int b = a + 1; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            double dx, dy, dz;
            const double r = dist_vec(a, b, dx, dy, dz);
            const double rcut = std::max(ea.r_cut, eb.r_cut);
            if (r > rcut)
                continue;
            const auto& sk = pair_sk(type[a], type[b]);
            fill_sk_block(ea, eb, sk, dx, dy, dz, r, H0, S, ia,
                          ao0[static_cast<std::size_t>(b)], false);
        }
    }

    std::vector<double> dq(static_cast<std::size_t>(n_at), 0.0);
    linalg::RealMatrix C = xt::zeros<double>({n_ao, n_ao});
    linalg::RealVector eps = xt::zeros<double>({n_ao});

    double e_band = 0;
    double e_scc = 0;

    for (int it = 0; it < max_scc; ++it) {
        // H = H0 + 0.5 S_μν (ΔV_A + ΔV_B)
        std::vector<double> dV(static_cast<std::size_t>(n_at), 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b) {
                double dx, dy, dz;
                const double r = (a == b) ? 0.0 : dist_vec(a, b, dx, dy, dz);
                const double g =
                    (a == b)
                        ? el[static_cast<std::size_t>(a)]->hubbard_u
                        : gamma_ab(r,
                                   el[static_cast<std::size_t>(a)]->hubbard_u,
                                   el[static_cast<std::size_t>(b)]->hubbard_u);
                shift += g * dq[static_cast<std::size_t>(b)];
            }
            dV[static_cast<std::size_t>(a)] = shift;
        }

        linalg::RealMatrix H = H0;
        for (int a = 0; a < n_at; ++a) {
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int b = 0; b < n_at; ++b) {
                const int ib = ao0[static_cast<std::size_t>(b)];
                const int nb = el[static_cast<std::size_t>(b)]->n_ao;
                const double half = 0.5
                                    * (dV[static_cast<std::size_t>(a)]
                                       + dV[static_cast<std::size_t>(b)]);
                for (int i = 0; i < na; ++i)
                    for (int j = 0; j < nb; ++j)
                        H(ia + i, ib + j) += half * S(ia + i, ib + j);
            }
        }

        auto sol = linalg::diagonalize_fock_generalized(H, S);
        eps = sol.first;
        C = sol.second;

        // Density P = 2 C_occ C_occ^T; Mulliken q_A = sum_{μ∈A} (P S)_μμ
        linalg::RealMatrix P = xt::zeros<double>({n_ao, n_ao});
        for (int o = 0; o < n_occ; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu)
                    P(mu, nu) += 2.0 * C(mu, o) * C(nu, o);

        linalg::RealMatrix PS = xt::linalg::dot(P, S);
        double max_dq = 0;
        for (int a = 0; a < n_at; ++a) {
            double qa = 0;
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int i = 0; i < na; ++i)
                qa += PS(ia + i, ia + i);
            const double dqi = qa - el[static_cast<std::size_t>(a)]->n_val;
            max_dq = std::max(max_dq,
                              std::abs(dqi - dq[static_cast<std::size_t>(a)]));
            dq[static_cast<std::size_t>(a)] = dqi;
        }

        e_band = 0;
        for (int o = 0; o < n_occ; ++o)
            e_band += 2.0 * eps(o);
        e_scc = 0;
        for (int a = 0; a < n_at; ++a)
            for (int b = 0; b < n_at; ++b) {
                double dx, dy, dz;
                const double r = (a == b) ? 0.0 : dist_vec(a, b, dx, dy, dz);
                const double g =
                    (a == b)
                        ? el[static_cast<std::size_t>(a)]->hubbard_u
                        : gamma_ab(r,
                                   el[static_cast<std::size_t>(a)]->hubbard_u,
                                   el[static_cast<std::size_t>(b)]->hubbard_u);
                e_scc += 0.5 * g * dq[static_cast<std::size_t>(a)]
                         * dq[static_cast<std::size_t>(b)];
            }
        if (max_dq < scc_tol)
            break;
    }

    double e_repulsive = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = a + 1; b < n_at; ++b) {
            double dx, dy, dz;
            const double r = dist_vec(a, b, dx, dy, dz);
            const auto& rep = pair_rep(type[a], type[b]);
            e_repulsive += rep.a * std::exp(-rep.b * r);
        }

    // E = sum n ε - E_scc + E_rep  (standard DFTB2 with band including shift)
    energy_ = e_band - e_scc + e_repulsive;
    energy_meta += energy_;




    if (with_forces) {
        // Closed-form nuclear gradient of E = ∑ n_i ε_i − E_SCC + E_rep
        // (Hellmann–Feynman + Pulay + repulsive). SCC γ geometry term is
        // carried by the Hamiltonian shift on ∂S (½(ΔV_A+ΔV_B)∂S) together with
        // −E_SCC bookkeeping at fixed Δq; residual vs full re-SCC ∇E is the
        // Mulliken charge-response on multi-element systems (H2 matches tightly).
        std::vector<double> dV(static_cast<std::size_t>(n_at), 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b) {
                double dx, dy, dz;
                const double r = (a == b) ? 0.0 : dist_vec(a, b, dx, dy, dz);
                const double g =
                    (a == b)
                        ? el[static_cast<std::size_t>(a)]->hubbard_u
                        : gamma_ab(r, el[static_cast<std::size_t>(a)]->hubbard_u,
                                   el[static_cast<std::size_t>(b)]->hubbard_u);
                shift += g * dq[static_cast<std::size_t>(b)];
            }
            dV[static_cast<std::size_t>(a)] = shift;
        }

        linalg::RealMatrix P = xt::zeros<double>({n_ao, n_ao});
        linalg::RealMatrix W = xt::zeros<double>({n_ao, n_ao});
        for (int o = 0; o < n_occ; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu) {
                    P(mu, nu) += 2.0 * C(mu, o) * C(nu, o);
                    W(mu, nu) += 2.0 * eps(o) * C(mu, o) * C(nu, o);
                }

        linalg::RealMatrix dH0[3] = {xt::zeros<double>({n_ao, n_ao}),
                                     xt::zeros<double>({n_ao, n_ao}),
                                     xt::zeros<double>({n_ao, n_ao})};
        linalg::RealMatrix dS[3] = {xt::zeros<double>({n_ao, n_ao}),
                                    xt::zeros<double>({n_ao, n_ao}),
                                    xt::zeros<double>({n_ao, n_ao})};

        for (int a = 0; a < n_at; ++a) {
            const auto& ea = *el[static_cast<std::size_t>(a)];
            const int ia = ao0[static_cast<std::size_t>(a)];
            for (int b = a + 1; b < n_at; ++b) {
                const auto& eb = *el[static_cast<std::size_t>(b)];
                double dx, dy, dz;
                const double r = dist_vec(a, b, dx, dy, dz);
                if (r < 1e-12)
                    continue;
                const int ib = ao0[static_cast<std::size_t>(b)];
                const double dir[3] = {dx / r, dy / r, dz / r};
                double dE_dbond[3] = {0, 0, 0};

                if (r <= std::max(ea.r_cut, eb.r_cut)) {
                    const auto& sk = pair_sk(type[a], type[b]);
                    for (int k = 0; k < 3; ++k) {
                        dH0[k] = xt::zeros<double>({n_ao, n_ao});
                        dS[k] = xt::zeros<double>({n_ao, n_ao});
                    }
                    fill_sk_block_deriv(ea, eb, sk, dx, dy, dz, r, dH0, dS, ia,
                                        ib);
                    const double shift = 0.5
                                         * (dV[static_cast<std::size_t>(a)]
                                            + dV[static_cast<std::size_t>(b)]);
                    const int na = ea.n_ao;
                    const int nb = eb.n_ao;
                    for (int k = 0; k < 3; ++k) {
                        double acc = 0;
                        for (int i = 0; i < na; ++i)
                            for (int j = 0; j < nb; ++j) {
                                const double p = P(ia + i, ib + j);
                                const double w = W(ia + i, ib + j);
                                const double dh = dH0[k](ia + i, ib + j);
                                const double ds = dS[k](ia + i, ib + j);
                                acc += p * (dh + shift * ds) - w * ds;
                            }
                        dE_dbond[k] += 2.0 * acc;
                    }
                }

                // Explicit ∂γ contribution for −E_SCC at fixed Δq:
                // E_SCC pair weight Δq_a Δq_b γ (a<b double-count resolved) ⇒
                // ∂E/∂bond from −E_SCC is −Δq_a Δq_b ∂γ/∂bond.
                const double dgdR = dgamma_dr(r, ea.hubbard_u, eb.hubbard_u);
                const double dqa = dq[static_cast<std::size_t>(a)];
                const double dqb = dq[static_cast<std::size_t>(b)];
                for (int k = 0; k < 3; ++k)
                    dE_dbond[k] += -dqa * dqb * dgdR * dir[k];

                const auto& rep = pair_rep(type[a], type[b]);
                const double dErep_dr =
                    -rep.a * rep.b * std::exp(-rep.b * r);
                for (int k = 0; k < 3; ++k)
                    dE_dbond[k] += dErep_dr * dir[k];

                fx[a] += static_cast<AccReal>(dE_dbond[0]);
                fy[a] += static_cast<AccReal>(dE_dbond[1]);
                fz[a] += static_cast<AccReal>(dE_dbond[2]);
                fx[b] -= static_cast<AccReal>(dE_dbond[0]);
                fy[b] -= static_cast<AccReal>(dE_dbond[1]);
                fz[b] -= static_cast<AccReal>(dE_dbond[2]);
            }
        }
    }

}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
