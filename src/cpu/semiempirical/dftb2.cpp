/// @file src/cpu/semiempirical/dftb2.cpp
/// SCC-DFTB2 multi-element: mio sketch (spin + analytic F) + 3ob-3-1 SK.
#include "atomiverse/cpu/semiempirical/dftb2.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_3ob_tables.h"
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

using basis_data::Dftb3obElement;
using basis_data::Dftb3obPair;
using basis_data::Dftb3obSkPoint;
using basis_data::kDftb3obElements;
using basis_data::kDftb3obNElements;
using basis_data::kDftb3obNPairs;
using basis_data::kDftb3obPairs;

// Hartree → eV; Bohr → Å (CODATA-aligned with common DFTB pins).
constexpr double kHa2eV = 27.211386245988;
constexpr double kBohr2A = 0.529177210903;
constexpr double kA2Bohr = 1.0 / kBohr2A;

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

/// @brief Fill  n_a × n_b  block of H0/S for atoms A–B separated by R vec.
void fill_sk_block(const DftbElementRow& ea, const DftbElementRow& eb,
                   const DftbPairSkRow& sk, double dx, double dy, double dz,
                   double r, linalg::RealMatrix& H0, linalg::RealMatrix& S,
                   int ia0, int ib0) {
    if (r < 1e-12)
        return;
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

    set_hs(0, 0, Hss, Sss);

    if (ea.n_ao >= 1 && eb.n_ao == 4) {
        set_hs(0, 1, l * Hsp, l * Ssp);
        set_hs(0, 2, m * Hsp, m * Ssp);
        set_hs(0, 3, n * Hsp, n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao >= 1) {
        set_hs(1, 0, -l * Hsp, -l * Ssp);
        set_hs(2, 0, -m * Hsp, -m * Ssp);
        set_hs(3, 0, -n * Hsp, -n * Ssp);
    }
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

/// dγ/dR = −R γ³.
double dgamma_dr(double r, double ua, double ub) {
    const double g = gamma_ab(r, ua, ub);
    return -r * g * g * g;
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

struct Occupancy {
    int n_alpha = 0;
    int n_beta = 0;
    bool unrestricted = false;
};

Occupancy resolve_occupancy(double n_elec, int multiplicity) {
    if (multiplicity < 1)
        throw std::runtime_error("DFTB2: multiplicity must be >= 1");
    if (std::fabs(n_elec - std::round(n_elec)) > 1e-9)
        throw std::runtime_error("DFTB2: non-integer valence electron count");
    const int n_e = static_cast<int>(std::lround(n_elec));
    if (n_e < 0)
        throw std::runtime_error("DFTB2: negative electron count after charge");
    const int n_unpaired = multiplicity - 1;
    if ((n_e - n_unpaired) % 2 != 0 || n_e < n_unpaired)
        throw std::runtime_error(
            "DFTB2: charge/multiplicity inconsistent with valence electron count");
    Occupancy occ;
    occ.n_alpha = (n_e + n_unpaired) / 2;
    occ.n_beta = (n_e - n_unpaired) / 2;
    occ.unrestricted = (multiplicity > 1) || (n_e % 2 != 0);
    return occ;
}

struct SccResult {
    double energy = 0;
    std::vector<double> dq;
    std::vector<double> spin; // m_A = Nα − Nβ
    // Restricted-path orbitals for closed-form forces.
    linalg::RealMatrix C;
    linalg::RealVector eps;
    int n_occ = 0;
    // Unrestricted-path orbitals (α/β) for collinear open-shell analytic F.
    linalg::RealMatrix Ca;
    linalg::RealMatrix Cb;
    linalg::RealVector eps_a;
    linalg::RealVector eps_b;
    int n_alpha = 0;
    int n_beta = 0;
    bool unrestricted = false;
};

/// Build H0/S once for current geometry.
void build_h0_s(int n_at, int n_ao, const std::vector<const DftbElementRow*>& el,
                const std::vector<int>& ao0, const int* type,
                const StorageReal* x, const StorageReal* y, const StorageReal* z,
                const cpu::Box& box, linalg::RealMatrix& H0,
                linalg::RealMatrix& S) {
    H0 = xt::zeros<double>({n_ao, n_ao});
    S = xt::zeros<double>({n_ao, n_ao});
    auto dist_vec = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = static_cast<double>(x[j] - x[i]);
        dy = static_cast<double>(y[j] - y[i]);
        dz = static_cast<double>(z[j] - z[i]);
        box.min_image(dx, dy, dz);
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };
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
            if (r > std::max(ea.r_cut, eb.r_cut))
                continue;
            const auto& sk = pair_sk(type[a], type[b]);
            fill_sk_block(ea, eb, sk, dx, dy, dz, r, H0, S, ia,
                          ao0[static_cast<std::size_t>(b)]);
        }
    }
}

double repulsive_energy(int n_at, const int* type, const StorageReal* x,
                        const StorageReal* y, const StorageReal* z,
                        const cpu::Box& box) {
    double e_rep = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = a + 1; b < n_at; ++b) {
            double dx = static_cast<double>(x[b] - x[a]);
            double dy = static_cast<double>(y[b] - y[a]);
            double dz = static_cast<double>(z[b] - z[a]);
            box.min_image(dx, dy, dz);
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            const auto& rep = pair_rep(type[a], type[b]);
            e_rep += rep.a * std::exp(-rep.b * r);
        }
    return e_rep;
}

void apply_shift(int n_at, const std::vector<const DftbElementRow*>& el,
                 const std::vector<int>& ao0, const linalg::RealMatrix& H0,
                 const linalg::RealMatrix& S, const std::vector<double>& eps_at,
                 linalg::RealMatrix& H) {
    H = H0;
    for (int a = 0; a < n_at; ++a) {
        const int ia = ao0[static_cast<std::size_t>(a)];
        const int na = el[static_cast<std::size_t>(a)]->n_ao;
        for (int b = 0; b < n_at; ++b) {
            const int ib = ao0[static_cast<std::size_t>(b)];
            const int nb = el[static_cast<std::size_t>(b)]->n_ao;
            const double half = 0.5
                                * (eps_at[static_cast<std::size_t>(a)]
                                   + eps_at[static_cast<std::size_t>(b)]);
            for (int i = 0; i < na; ++i)
                for (int j = 0; j < nb; ++j)
                    H(ia + i, ib + j) += half * S(ia + i, ib + j);
        }
    }
}

double gamma_pair(int a, int b, const std::vector<const DftbElementRow*>& el,
                  const StorageReal* x, const StorageReal* y,
                  const StorageReal* z, const cpu::Box& box) {
    if (a == b)
        return el[static_cast<std::size_t>(a)]->hubbard_u;
    double dx = static_cast<double>(x[b] - x[a]);
    double dy = static_cast<double>(y[b] - y[a]);
    double dz = static_cast<double>(z[b] - z[a]);
    box.min_image(dx, dy, dz);
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    return gamma_ab(r, el[static_cast<std::size_t>(a)]->hubbard_u,
                    el[static_cast<std::size_t>(b)]->hubbard_u);
}

SccResult scc_energy_restricted(
    int n_at, int n_ao, int n_occ, const std::vector<const DftbElementRow*>& el,
    const std::vector<int>& ao0, const int* type, const StorageReal* x,
    const StorageReal* y, const StorageReal* z, const cpu::Box& box,
    double scc_tol, int max_scc) {
    linalg::RealMatrix H0, S;
    build_h0_s(n_at, n_ao, el, ao0, type, x, y, z, box, H0, S);

    std::vector<double> dq(static_cast<std::size_t>(n_at), 0.0);
    linalg::RealMatrix C = xt::zeros<double>({n_ao, n_ao});
    linalg::RealVector eps = xt::zeros<double>({n_ao});
    double e_band = 0, e_scc = 0;

    for (int it = 0; it < max_scc; ++it) {
        std::vector<double> dV(static_cast<std::size_t>(n_at), 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b)
                shift += gamma_pair(a, b, el, x, y, z, box)
                         * dq[static_cast<std::size_t>(b)];
            dV[static_cast<std::size_t>(a)] = shift;
        }
        linalg::RealMatrix H;
        apply_shift(n_at, el, ao0, H0, S, dV, H);
        auto sol = linalg::diagonalize_fock_generalized(H, S);
        eps = sol.first;
        C = sol.second;

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
            for (int b = 0; b < n_at; ++b)
                e_scc += 0.5 * gamma_pair(a, b, el, x, y, z, box)
                         * dq[static_cast<std::size_t>(a)]
                         * dq[static_cast<std::size_t>(b)];
        if (max_dq < scc_tol)
            break;
    }

    SccResult out;
    out.energy = e_band - e_scc + repulsive_energy(n_at, type, x, y, z, box);
    out.dq = std::move(dq);
    out.spin.assign(static_cast<std::size_t>(n_at), 0.0);
    out.C = std::move(C);
    out.eps = std::move(eps);
    out.n_occ = n_occ;
    return out;
}

SccResult scc_energy_unrestricted(
    int n_at, int n_ao, int n_alpha, int n_beta,
    const std::vector<const DftbElementRow*>& el, const std::vector<int>& ao0,
    const int* type, const StorageReal* x, const StorageReal* y,
    const StorageReal* z, const cpu::Box& box, double scc_tol, int max_scc) {
    if (n_alpha > n_ao || n_beta > n_ao)
        throw std::runtime_error(
            "DFTB2: n_alpha/n_beta exceeds basis size");

    linalg::RealMatrix H0, S;
    build_h0_s(n_at, n_ao, el, ao0, type, x, y, z, box, H0, S);

    std::vector<double> dq(static_cast<std::size_t>(n_at), 0.0);
    std::vector<double> mspin(static_cast<std::size_t>(n_at), 0.0);
    // Seed high-spin magnetization so the first iteration is not spin-null.
    if (n_alpha != n_beta) {
        const double m_tot = static_cast<double>(n_alpha - n_beta);
        for (int a = 0; a < n_at; ++a)
            mspin[static_cast<std::size_t>(a)] =
                m_tot / static_cast<double>(n_at);
    }

    double e_band = 0, e_scc = 0, e_spin = 0;
    linalg::RealMatrix CA = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix CB = xt::zeros<double>({n_ao, n_ao});
    linalg::RealVector epsA = xt::zeros<double>({n_ao});
    linalg::RealVector epsB = xt::zeros<double>({n_ao});

    for (int it = 0; it < max_scc; ++it) {
        std::vector<double> dVc(static_cast<std::size_t>(n_at), 0.0);
        std::vector<double> dVs(static_cast<std::size_t>(n_at), 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b)
                shift += gamma_pair(a, b, el, x, y, z, box)
                         * dq[static_cast<std::size_t>(b)];
            dVc[static_cast<std::size_t>(a)] = shift;
            dVs[static_cast<std::size_t>(a)] =
                el[static_cast<std::size_t>(a)]->spin_w
                * mspin[static_cast<std::size_t>(a)];
        }

        // ε^σ_A = dVc + σ dVs  (σ = +1 α, −1 β)
        std::vector<double> eps_shift_a(static_cast<std::size_t>(n_at));
        std::vector<double> eps_shift_b(static_cast<std::size_t>(n_at));
        for (int a = 0; a < n_at; ++a) {
            eps_shift_a[static_cast<std::size_t>(a)] =
                dVc[static_cast<std::size_t>(a)]
                + dVs[static_cast<std::size_t>(a)];
            eps_shift_b[static_cast<std::size_t>(a)] =
                dVc[static_cast<std::size_t>(a)]
                - dVs[static_cast<std::size_t>(a)];
        }

        linalg::RealMatrix Ha, Hb;
        apply_shift(n_at, el, ao0, H0, S, eps_shift_a, Ha);
        apply_shift(n_at, el, ao0, H0, S, eps_shift_b, Hb);

        auto sol_a = linalg::diagonalize_fock_generalized(Ha, S);
        auto sol_b = linalg::diagonalize_fock_generalized(Hb, S);
        epsA = sol_a.first;
        CA = sol_a.second;
        epsB = sol_b.first;
        CB = sol_b.second;

        linalg::RealMatrix Pa = xt::zeros<double>({n_ao, n_ao});
        linalg::RealMatrix Pb = xt::zeros<double>({n_ao, n_ao});
        for (int o = 0; o < n_alpha; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu)
                    Pa(mu, nu) += CA(mu, o) * CA(nu, o);
        for (int o = 0; o < n_beta; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu)
                    Pb(mu, nu) += CB(mu, o) * CB(nu, o);

        linalg::RealMatrix PSa = xt::linalg::dot(Pa, S);
        linalg::RealMatrix PSb = xt::linalg::dot(Pb, S);

        double max_res = 0;
        for (int a = 0; a < n_at; ++a) {
            double qa = 0, qb = 0;
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int i = 0; i < na; ++i) {
                qa += PSa(ia + i, ia + i);
                qb += PSb(ia + i, ia + i);
            }
            const double dqi =
                (qa + qb) - el[static_cast<std::size_t>(a)]->n_val;
            const double mi = qa - qb;
            max_res = std::max(
                max_res, std::abs(dqi - dq[static_cast<std::size_t>(a)]));
            max_res = std::max(
                max_res, std::abs(mi - mspin[static_cast<std::size_t>(a)]));
            dq[static_cast<std::size_t>(a)] = dqi;
            mspin[static_cast<std::size_t>(a)] = mi;
        }

        e_band = 0;
        for (int o = 0; o < n_alpha; ++o)
            e_band += epsA(o);
        for (int o = 0; o < n_beta; ++o)
            e_band += epsB(o);
        e_scc = 0;
        for (int a = 0; a < n_at; ++a)
            for (int b = 0; b < n_at; ++b)
                e_scc += 0.5 * gamma_pair(a, b, el, x, y, z, box)
                         * dq[static_cast<std::size_t>(a)]
                         * dq[static_cast<std::size_t>(b)];
        e_spin = 0;
        for (int a = 0; a < n_at; ++a)
            e_spin += 0.5 * el[static_cast<std::size_t>(a)]->spin_w
                      * mspin[static_cast<std::size_t>(a)]
                      * mspin[static_cast<std::size_t>(a)];
        if (max_res < scc_tol)
            break;
    }

    SccResult out;
    out.energy = e_band - e_scc - e_spin
                 + repulsive_energy(n_at, type, x, y, z, box);
    out.dq = std::move(dq);
    out.spin = std::move(mspin);
    out.Ca = std::move(CA);
    out.Cb = std::move(CB);
    out.eps_a = std::move(epsA);
    out.eps_b = std::move(epsB);
    out.n_alpha = n_alpha;
    out.n_beta = n_beta;
    out.unrestricted = true;
    return out;
}


/// Closed-form nuclear gradient for **collinear unrestricted** SCC-DFTB2
/// (fixed Δq, m; dual α/β densities). Spin W is geometry-independent.
void analytic_forces_unrestricted(
    int n_at, int n_ao, int n_alpha, int n_beta,
    const std::vector<const DftbElementRow*>& el, const std::vector<int>& ao0,
    const int* type, const StorageReal* x, const StorageReal* y,
    const StorageReal* z, const cpu::Box& box, const std::vector<double>& dq,
    const std::vector<double>& mspin, const linalg::RealMatrix& Ca,
    const linalg::RealMatrix& Cb, const linalg::RealVector& eps_a_orb,
    const linalg::RealVector& eps_b_orb, AccReal* fx, AccReal* fy,
    AccReal* fz) {
    auto dist_vec = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = static_cast<double>(x[j] - x[i]);
        dy = static_cast<double>(y[j] - y[i]);
        dz = static_cast<double>(z[j] - z[i]);
        box.min_image(dx, dy, dz);
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    // Atom shifts ε^σ_A = Σ_B γ_AB Δq_B ± W_A m_A
    std::vector<double> eps_at_a(static_cast<std::size_t>(n_at), 0.0);
    std::vector<double> eps_at_b(static_cast<std::size_t>(n_at), 0.0);
    for (int a = 0; a < n_at; ++a) {
        double shift = 0;
        for (int b = 0; b < n_at; ++b)
            shift += gamma_pair(a, b, el, x, y, z, box)
                     * dq[static_cast<std::size_t>(b)];
        const double vs = el[static_cast<std::size_t>(a)]->spin_w
                          * mspin[static_cast<std::size_t>(a)];
        eps_at_a[static_cast<std::size_t>(a)] = shift + vs;
        eps_at_b[static_cast<std::size_t>(a)] = shift - vs;
    }

    linalg::RealMatrix Pa = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix Pb = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix Wa = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix Wb = xt::zeros<double>({n_ao, n_ao});
    for (int o = 0; o < n_alpha; ++o)
        for (int mu = 0; mu < n_ao; ++mu)
            for (int nu = 0; nu < n_ao; ++nu) {
                Pa(mu, nu) += Ca(mu, o) * Ca(nu, o);
                Wa(mu, nu) += eps_a_orb(o) * Ca(mu, o) * Ca(nu, o);
            }
    for (int o = 0; o < n_beta; ++o)
        for (int mu = 0; mu < n_ao; ++mu)
            for (int nu = 0; nu < n_ao; ++nu) {
                Pb(mu, nu) += Cb(mu, o) * Cb(nu, o);
                Wb(mu, nu) += eps_b_orb(o) * Cb(mu, o) * Cb(nu, o);
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
                fill_sk_block_deriv(ea, eb, sk, dx, dy, dz, r, dH0, dS, ia, ib);
                const double shift_a =
                    0.5 * (eps_at_a[static_cast<std::size_t>(a)]
                           + eps_at_a[static_cast<std::size_t>(b)]);
                const double shift_b =
                    0.5 * (eps_at_b[static_cast<std::size_t>(a)]
                           + eps_at_b[static_cast<std::size_t>(b)]);
                const int na = ea.n_ao;
                const int nb = eb.n_ao;
                for (int k = 0; k < 3; ++k) {
                    double acc = 0;
                    for (int i = 0; i < na; ++i)
                        for (int j = 0; j < nb; ++j) {
                            const double dh = dH0[k](ia + i, ib + j);
                            const double ds = dS[k](ia + i, ib + j);
                            acc += Pa(ia + i, ib + j) * (dh + shift_a * ds)
                                   - Wa(ia + i, ib + j) * ds;
                            acc += Pb(ia + i, ib + j) * (dh + shift_b * ds)
                                   - Wb(ia + i, ib + j) * ds;
                        }
                    // Factor 2: only the A–B upper block is contracted (Hermitian).
                    dE_dbond[k] += 2.0 * acc;
                }
            }

            const double dgdR = dgamma_dr(r, ea.hubbard_u, eb.hubbard_u);
            const double dqa = dq[static_cast<std::size_t>(a)];
            const double dqb = dq[static_cast<std::size_t>(b)];
            for (int k = 0; k < 3; ++k)
                dE_dbond[k] += -dqa * dqb * dgdR * dir[k];

            // Spin energy −½ Σ W m²: W constant, m fixed ⇒ no geometry term.

            const auto& rep = pair_rep(type[a], type[b]);
            const double dErep_dr = -rep.a * rep.b * std::exp(-rep.b * r);
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

/// Closed-form nuclear gradient for **restricted** SCC-DFTB2 (fixed Δq).
void analytic_forces_restricted(
    int n_at, int n_ao, int n_occ, const std::vector<const DftbElementRow*>& el,
    const std::vector<int>& ao0, const int* type, const StorageReal* x,
    const StorageReal* y, const StorageReal* z, const cpu::Box& box,
    const std::vector<double>& dq, const linalg::RealMatrix& C,
    const linalg::RealVector& eps, AccReal* fx, AccReal* fy, AccReal* fz) {
    auto dist_vec = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = static_cast<double>(x[j] - x[i]);
        dy = static_cast<double>(y[j] - y[i]);
        dz = static_cast<double>(z[j] - z[i]);
        box.min_image(dx, dy, dz);
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    std::vector<double> dV(static_cast<std::size_t>(n_at), 0.0);
    for (int a = 0; a < n_at; ++a) {
        double shift = 0;
        for (int b = 0; b < n_at; ++b)
            shift += gamma_pair(a, b, el, x, y, z, box)
                     * dq[static_cast<std::size_t>(b)];
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
                fill_sk_block_deriv(ea, eb, sk, dx, dy, dz, r, dH0, dS, ia, ib);
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

            const double dgdR = dgamma_dr(r, ea.hubbard_u, eb.hubbard_u);
            const double dqa = dq[static_cast<std::size_t>(a)];
            const double dqb = dq[static_cast<std::size_t>(b)];
            for (int k = 0; k < 3; ++k)
                dE_dbond[k] += -dqa * dqb * dgdR * dir[k];

            const auto& rep = pair_rep(type[a], type[b]);
            const double dErep_dr = -rep.a * rep.b * std::exp(-rep.b * r);
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


SccResult scc_total(int n_at, int n_ao, const Occupancy& occ,
                    const std::vector<const DftbElementRow*>& el,
                    const std::vector<int>& ao0, const int* type,
                    const StorageReal* x, const StorageReal* y,
                    const StorageReal* z, const cpu::Box& box, double scc_tol,
                    int max_scc) {
    if (!occ.unrestricted)
        return scc_energy_restricted(n_at, n_ao, occ.n_alpha, el, ao0, type, x,
                                     y, z, box, scc_tol, max_scc);
    return scc_energy_unrestricted(n_at, n_ao, occ.n_alpha, occ.n_beta, el, ao0,
                                   type, x, y, z, box, scc_tol, max_scc);
}


double launch_mio(cpu::Frame& f, double scc_tol, int max_scc, bool with_forces,
                  int charge, int multiplicity, std::vector<double>& out_dq,
                  std::vector<double>& out_spin) {
    const int n_at = f.nlocal();
    const int* type = f.span(kAtom_type).data();
    StorageReal* x = f.span(kAtom_x);
    StorageReal* y = f.span(kAtom_y);
    StorageReal* z = f.span(kAtom_z);
    AccReal* fx = f.span(kAtom_fx);
    AccReal* fy = f.span(kAtom_fy);
    AccReal* fz = f.span(kAtom_fz);
    double& energy_meta = f.meta.get(kEnergy);

    std::vector<const DftbElementRow*> el(static_cast<std::size_t>(n_at));
    std::vector<int> ao0(static_cast<std::size_t>(n_at) + 1, 0);
    double n_elec = -static_cast<double>(charge);
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

    const Occupancy occ = resolve_occupancy(n_elec, multiplicity);
    auto res = scc_total(n_at, n_ao, occ, el, ao0, type, x, y, z, f.box,
                         scc_tol, max_scc);
    out_dq = res.dq;
    out_spin = res.spin;
    const double energy = res.energy;
    energy_meta += energy;

    if (with_forces) {
        if (occ.unrestricted && res.Ca.shape(0) > 0) {
            analytic_forces_unrestricted(
                n_at, n_ao, res.n_alpha, res.n_beta, el, ao0, type, x, y, z,
                f.box, res.dq, res.spin, res.Ca, res.Cb, res.eps_a, res.eps_b,
                fx, fy, fz);
        } else if (!occ.unrestricted && res.n_occ > 0 && res.C.shape(0) > 0) {
            analytic_forces_restricted(n_at, n_ao, res.n_occ, el, ao0, type, x,
                                       y, z, f.box, res.dq, res.C, res.eps, fx,
                                       fy, fz);
        }
        // else: no orbitals (empty system) — leave forces unchanged
    }
    return energy;
}

// ── 3ob-3-1 path (internal Ha / Bohr; public eV / Å) ───────────────────

const Dftb3obElement& element_3ob(int z) {
    for (std::size_t i = 0; i < kDftb3obNElements; ++i)
        if (kDftb3obElements[i].z == z)
            return kDftb3obElements[i];
    throw std::runtime_error("DFTB2(ThreeOb): unsupported element Z="
                             + std::to_string(z) + " (H/C/N/O only)");
}

const Dftb3obPair& pair_3ob(int za, int zb) {
    for (std::size_t i = 0; i < kDftb3obNPairs; ++i)
        if (kDftb3obPairs[i].z_a == za && kDftb3obPairs[i].z_b == zb)
            return kDftb3obPairs[i];
    throw std::runtime_error("DFTB2(ThreeOb): no SK table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

/// Linear interp of SK channels at r [Bohr]. Returns Hss,Hsp,Hpps,Hppp,Sss,...
void interp_3ob_sk(const Dftb3obPair& p, double r_bohr, double out[8]) {
    if (r_bohr <= p.grid) {
        const Dftb3obSkPoint& pt = p.sk[0];
        out[0] = pt.h_ss;
        out[1] = pt.h_sp;
        out[2] = pt.h_pps;
        out[3] = pt.h_ppp;
        out[4] = pt.s_ss;
        out[5] = pt.s_sp;
        out[6] = pt.s_pps;
        out[7] = pt.s_ppp;
        return;
    }
    const double x = r_bohr / p.grid; // 1-based continuous index
    if (x >= static_cast<double>(p.n_grid)) {
        for (int i = 0; i < 8; ++i)
            out[i] = 0;
        return;
    }
    int i0 = static_cast<int>(std::floor(x)) - 1;
    if (i0 < 0)
        i0 = 0;
    if (i0 > p.n_grid - 2)
        i0 = p.n_grid - 2;
    const int i1 = i0 + 1;
    const double t = x - static_cast<double>(i0 + 1);
    const Dftb3obSkPoint& a = p.sk[i0];
    const Dftb3obSkPoint& b = p.sk[i1];
    const double u = 1.0 - t;
    out[0] = u * a.h_ss + t * b.h_ss;
    out[1] = u * a.h_sp + t * b.h_sp;
    out[2] = u * a.h_pps + t * b.h_pps;
    out[3] = u * a.h_ppp + t * b.h_ppp;
    out[4] = u * a.s_ss + t * b.s_ss;
    out[5] = u * a.s_sp + t * b.s_sp;
    out[6] = u * a.s_pps + t * b.s_pps;
    out[7] = u * a.s_ppp + t * b.s_ppp;
}

double rep_3ob(const Dftb3obPair& p, double r_bohr) {
    if (r_bohr >= p.rep_cutoff)
        return 0.0;
    // First interval start
    const double r0 = p.rep_int[0];
    if (r_bohr < r0)
        return std::exp(-p.rep_a1 * r_bohr + p.rep_a2) + p.rep_a3;
    for (int i = 0; i < p.n_rep_int; ++i) {
        const double* iv = p.rep_int + static_cast<std::size_t>(i) * 8;
        const double ra = iv[0];
        const double rb = iv[1];
        if (r_bohr < rb || std::abs(r_bohr - rb) < 1e-14) {
            const double dr = r_bohr - ra;
            double e = 0;
            double pow = 1.0;
            for (int k = 0; k < 6; ++k) {
                e += iv[2 + k] * pow;
                pow *= dr;
            }
            return e;
        }
    }
    return 0.0;
}

void fill_sk_block_3ob(const Dftb3obElement& ea, const Dftb3obElement& eb,
                       int za, int zb, double dx, double dy, double dz, double r,
                       linalg::RealMatrix& H0, linalg::RealMatrix& S, int ia0,
                       int ib0) {
    if (r < 1e-12)
        return;
    const double inv_r = 1.0 / r;
    const double l = dx * inv_r;
    const double m = dy * inv_r;
    const double n = dz * inv_r;

    const Dftb3obPair& p_ab = pair_3ob(za, zb);
    double ch[8];
    interp_3ob_sk(p_ab, r, ch);
    double Hss = ch[0], Hsp = ch[1], Hpps = ch[2], Hppp = ch[3];
    double Sss = ch[4], Ssp = ch[5], Spps = ch[6], Sppp = ch[7];
    // X–H tables store only ss; s–p lives in H–X.
    if (za != 1 && zb == 1) {
        const Dftb3obPair& p_ba = pair_3ob(zb, za);
        double rev[8];
        interp_3ob_sk(p_ba, r, rev);
        Hsp = rev[1];
        Ssp = rev[5];
    }

    auto set_hs = [&](int i, int j, double hv, double sv) {
        H0(ia0 + i, ib0 + j) = hv;
        H0(ib0 + j, ia0 + i) = hv;
        S(ia0 + i, ib0 + j) = sv;
        S(ib0 + j, ia0 + i) = sv;
    };

    set_hs(0, 0, Hss, Sss);
    if (ea.n_ao >= 1 && eb.n_ao == 4) {
        set_hs(0, 1, l * Hsp, l * Ssp);
        set_hs(0, 2, m * Hsp, m * Ssp);
        set_hs(0, 3, n * Hsp, n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao >= 1) {
        set_hs(1, 0, -l * Hsp, -l * Ssp);
        set_hs(2, 0, -m * Hsp, -m * Ssp);
        set_hs(3, 0, -n * Hsp, -n * Ssp);
    }
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

double gamma_ab_ha(double r_bohr, double ua_ha, double ub_ha) {
    if (r_bohr < 1e-14)
        return ua_ha;
    const double tau = 0.5 * (1.0 / ua_ha + 1.0 / ub_ha);
    return 1.0 / std::sqrt(r_bohr * r_bohr + tau * tau);
}

double energy_3ob_ha(const std::vector<const Dftb3obElement*>& el,
                     const std::vector<int>& ao0, const std::vector<int>& Z,
                     const std::vector<double>& x_b,
                     const std::vector<double>& y_b,
                     const std::vector<double>& z_b, int n_at, int n_ao,
                     int n_occ, double scc_tol, int max_scc) {
    auto dist = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = x_b[static_cast<std::size_t>(j)] - x_b[static_cast<std::size_t>(i)];
        dy = y_b[static_cast<std::size_t>(j)] - y_b[static_cast<std::size_t>(i)];
        dz = z_b[static_cast<std::size_t>(j)] - z_b[static_cast<std::size_t>(i)];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    linalg::RealMatrix H0 = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix S = xt::zeros<double>({n_ao, n_ao});

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
            const double r = dist(a, b, dx, dy, dz);
            if (r > std::max(ea.r_cut, eb.r_cut))
                continue;
            fill_sk_block_3ob(ea, eb, Z[static_cast<std::size_t>(a)],
                              Z[static_cast<std::size_t>(b)], dx, dy, dz, r, H0,
                              S, ia, ao0[static_cast<std::size_t>(b)]);
        }
    }

    std::vector<double> dq(static_cast<std::size_t>(n_at), 0.0);
    linalg::RealMatrix P = xt::zeros<double>({n_ao, n_ao});
    double e_scc = 0;

    for (int it = 0; it < max_scc; ++it) {
        std::vector<double> dV(static_cast<std::size_t>(n_at), 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b) {
                double dx, dy, dz;
                const double r = (a == b) ? 0.0 : dist(a, b, dx, dy, dz);
                const double g =
                    (a == b)
                        ? el[static_cast<std::size_t>(a)]->hubbard_u
                        : gamma_ab_ha(r,
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
        auto C = sol.second;

        P = xt::zeros<double>({n_ao, n_ao});
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

        e_scc = 0;
        for (int a = 0; a < n_at; ++a)
            for (int b = 0; b < n_at; ++b) {
                double dx, dy, dz;
                const double r = (a == b) ? 0.0 : dist(a, b, dx, dy, dz);
                const double g =
                    (a == b)
                        ? el[static_cast<std::size_t>(a)]->hubbard_u
                        : gamma_ab_ha(r,
                                      el[static_cast<std::size_t>(a)]->hubbard_u,
                                      el[static_cast<std::size_t>(b)]->hubbard_u);
                e_scc += 0.5 * g * dq[static_cast<std::size_t>(a)]
                         * dq[static_cast<std::size_t>(b)];
            }
        if (max_dq < scc_tol)
            break;
    }

    // E_elec = Tr(P H0) + E_scc  (DFTB2; matches DFTB+ double-counting)
    double e_h0 = 0;
    for (int mu = 0; mu < n_ao; ++mu)
        for (int nu = 0; nu < n_ao; ++nu)
            e_h0 += P(mu, nu) * H0(mu, nu);

    double e_repulsive = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = a + 1; b < n_at; ++b) {
            double dx, dy, dz;
            const double r = dist(a, b, dx, dy, dz);
            const auto& pr = pair_3ob(Z[static_cast<std::size_t>(a)],
                                      Z[static_cast<std::size_t>(b)]);
            e_repulsive += rep_3ob(pr, r);
        }

    return e_h0 + e_scc + e_repulsive;
}

double launch_3ob(DFTB2& self, cpu::Frame& f, double scc_tol, int max_scc,
                  bool with_forces) {
    const int n_at = f.nlocal();
    const auto type = f.span(kAtom_type);
    StorageReal* x = f.span(kAtom_x);
    StorageReal* y = f.span(kAtom_y);
    StorageReal* z = f.span(kAtom_z);
    AccReal* fx = f.span(kAtom_fx);
    AccReal* fy = f.span(kAtom_fy);
    AccReal* fz = f.span(kAtom_fz);
    double& energy_meta = f.meta.get(kEnergy);

    std::vector<const Dftb3obElement*> el(static_cast<std::size_t>(n_at));
    std::vector<int> ao0(static_cast<std::size_t>(n_at) + 1, 0);
    std::vector<int> Z(static_cast<std::size_t>(n_at));
    std::vector<double> xb(static_cast<std::size_t>(n_at));
    std::vector<double> yb(static_cast<std::size_t>(n_at));
    std::vector<double> zb(static_cast<std::size_t>(n_at));
    double n_elec = 0;
    for (int a = 0; a < n_at; ++a) {
        Z[static_cast<std::size_t>(a)] = type[a];
        el[static_cast<std::size_t>(a)] = &element_3ob(type[a]);
        n_elec += el[static_cast<std::size_t>(a)]->n_val;
        ao0[static_cast<std::size_t>(a) + 1] =
            ao0[static_cast<std::size_t>(a)]
            + el[static_cast<std::size_t>(a)]->n_ao;
        // Å → Bohr for SK/rep (molecular; no min-image in 3ob pin path)
        xb[static_cast<std::size_t>(a)] = static_cast<double>(x[a]) * kA2Bohr;
        yb[static_cast<std::size_t>(a)] = static_cast<double>(y[a]) * kA2Bohr;
        zb[static_cast<std::size_t>(a)] = static_cast<double>(z[a]) * kA2Bohr;
    }
    const int n_ao = ao0[static_cast<std::size_t>(n_at)];
    if (n_ao < 1)
        throw std::runtime_error("DFTB2(ThreeOb): no basis functions");
    if (std::fabs(n_elec - std::round(n_elec)) > 1e-9
        || static_cast<int>(std::lround(n_elec)) % 2 != 0)
        throw std::runtime_error(
            "DFTB2(ThreeOb): odd/non-integer valence electron count");
    const int n_occ = static_cast<int>(std::lround(n_elec)) / 2;

    const double e_ha =
        energy_3ob_ha(el, ao0, Z, xb, yb, zb, n_at, n_ao, n_occ, scc_tol,
                      max_scc);
    const double energy = e_ha * kHa2eV;
    energy_meta += energy;

    if (with_forces) {
        constexpr double h_fd_A = 1.0e-4; // Å
        const double h_fd_B = h_fd_A * kA2Bohr;
        for (int a = 0; a < n_at; ++a) {
            for (int k = 0; k < 3; ++k) {
                double* col = (k == 0)   ? xb.data()
                              : (k == 1) ? yb.data()
                                         : zb.data();
                const double saved = col[static_cast<std::size_t>(a)];
                col[static_cast<std::size_t>(a)] = saved + h_fd_B;
                const double ep =
                    energy_3ob_ha(el, ao0, Z, xb, yb, zb, n_at, n_ao, n_occ,
                                  scc_tol, max_scc)
                    * kHa2eV;
                col[static_cast<std::size_t>(a)] = saved - h_fd_B;
                const double em =
                    energy_3ob_ha(el, ao0, Z, xb, yb, zb, n_at, n_ao, n_occ,
                                  scc_tol, max_scc)
                    * kHa2eV;
                col[static_cast<std::size_t>(a)] = saved;
                const double fi = -(ep - em) / (2.0 * h_fd_A);
                if (k == 0)
                    fx[a] += static_cast<AccReal>(fi);
                else if (k == 1)
                    fy[a] += static_cast<AccReal>(fi);
                else
                    fz[a] += static_cast<AccReal>(fi);
            }
        }
    }
    (void)self;
    return energy;
}


} // namespace

void DFTB2::build(const BuildContext& /*ctx*/) {}

void DFTB2::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    if (f.nlocal() < 1)
        throw std::invalid_argument("DFTB2: empty frame");

    if (sk_set == DftbSkSet::ThreeOb) {
        if (multiplicity != 1 || charge != 0)
            throw std::runtime_error(
                "DFTB2(ThreeOb): open-shell/charge not supported on this path");
        energy_ = launch_3ob(*this, f, scc_tol, max_scc, with_forces);
        last_dq_.clear();
        last_spin_.assign(static_cast<std::size_t>(f.nlocal()), 0.0);
    } else {
        energy_ = launch_mio(f, scc_tol, max_scc, with_forces, charge,
                             multiplicity, last_dq_, last_spin_);
    }
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
