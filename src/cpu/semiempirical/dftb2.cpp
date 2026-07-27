/// @file src/cpu/semiempirical/dftb2.cpp
/// SCC-DFTB2 multi-element (H/C/N/O) with static mio-style SK tables + collinear spin.
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
        auto eps = sol.first;
        auto C = sol.second;

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
        std::vector<double> eps_a(static_cast<std::size_t>(n_at));
        std::vector<double> eps_b(static_cast<std::size_t>(n_at));
        for (int a = 0; a < n_at; ++a) {
            eps_a[static_cast<std::size_t>(a)] =
                dVc[static_cast<std::size_t>(a)]
                + dVs[static_cast<std::size_t>(a)];
            eps_b[static_cast<std::size_t>(a)] =
                dVc[static_cast<std::size_t>(a)]
                - dVs[static_cast<std::size_t>(a)];
        }

        linalg::RealMatrix Ha, Hb;
        apply_shift(n_at, el, ao0, H0, S, eps_a, Ha);
        apply_shift(n_at, el, ao0, H0, S, eps_b, Hb);

        auto sol_a = linalg::diagonalize_fock_generalized(Ha, S);
        auto sol_b = linalg::diagonalize_fock_generalized(Hb, S);
        auto epsA = sol_a.first;
        auto CA = sol_a.second;
        auto epsB = sol_b.first;
        auto CB = sol_b.second;

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
    return out;
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

} // namespace

void DFTB2::build(const BuildContext& /*ctx*/) {}

void DFTB2::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n_at = f.nlocal();
    if (n_at < 1)
        throw std::invalid_argument("DFTB2: empty frame");

    const int* type = f.span(kAtom_type).data();
    const StorageReal* x = f.span(kAtom_x);
    const StorageReal* y = f.span(kAtom_y);
    const StorageReal* z = f.span(kAtom_z);
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
    energy_ = res.energy;
    last_dq_ = res.dq;
    last_spin_ = res.spin;
    energy_meta += energy_;

    if (with_forces) {
        constexpr double h_fd = 1.0e-4; // Å
        StorageReal* mx = f.span(kAtom_x);
        StorageReal* my = f.span(kAtom_y);
        StorageReal* mz = f.span(kAtom_z);
        for (int a = 0; a < n_at; ++a) {
            for (int k = 0; k < 3; ++k) {
                StorageReal* col = (k == 0) ? mx : (k == 1) ? my : mz;
                const StorageReal saved = col[a];
                col[a] = saved + static_cast<StorageReal>(h_fd);
                const double ep =
                    scc_total(n_at, n_ao, occ, el, ao0, type, mx, my, mz,
                              f.box, scc_tol, max_scc)
                        .energy;
                col[a] = saved - static_cast<StorageReal>(h_fd);
                const double em =
                    scc_total(n_at, n_ao, occ, el, ao0, type, mx, my, mz,
                              f.box, scc_tol, max_scc)
                        .energy;
                col[a] = saved;
                const double fi = -(ep - em) / (2.0 * h_fd);
                if (k == 0)
                    fx[a] += static_cast<AccReal>(fi);
                else if (k == 1)
                    fy[a] += static_cast<AccReal>(fi);
                else
                    fz[a] += static_cast<AccReal>(fi);
            }
        }
    }
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
