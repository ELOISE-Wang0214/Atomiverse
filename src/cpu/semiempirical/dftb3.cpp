/// @file src/cpu/semiempirical/dftb3.cpp
/// SCC-DFTB3 multi-element (H/C/N/O) with static mio-style SK tables.
/// H–H reads the shared transcribed mio-1-1 grid (`internal/mio_hh_table.h`),
/// the same one DFTB0/DFTB2 read, so DFTB3 reduces to DFTB2 exactly when the
/// charge fluctuations vanish.
#include "atomiverse/cpu/semiempirical/dftb3.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_mio_tables.h"
#include "atomiverse/cpu/semiempirical/internal/dftb3_scc.h"
#include "atomiverse/cpu/semiempirical/internal/mio_hh_table.h"
#include "atomiverse/cpu/semiempirical/internal/scc_mixer.h"

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
    throw std::runtime_error("DFTB3: unsupported element Z=" + std::to_string(z)
                             + " (static mio-style tables cover H/C/N/O only)");
}

const DftbPairSkRow& pair_sk(int za, int zb) {
    for (std::size_t i = 0; i < kDftbNPairSk; ++i)
        if (kDftbPairSk[i].z_a == za && kDftbPairSk[i].z_b == zb)
            return kDftbPairSk[i];
    throw std::runtime_error("DFTB3: no SK table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

const DftbPairRepRow& pair_rep(int za, int zb) {
    for (std::size_t i = 0; i < kDftbNPairRep; ++i)
        if (kDftbPairRep[i].z_a == za && kDftbPairRep[i].z_b == zb)
            return kDftbPairRep[i];
    throw std::runtime_error("DFTB3: no repulsive table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

double h_radial(double v0, double r, double r0, double d) {
    return v0 * std::exp(-(r - r0) / d);
}
double s_radial(double s0, double r, double r0, double w) {
    const double t = (r - r0) / w;
    return s0 * std::exp(-0.5 * t * t);
}

/// H–H reaches as far as the transcribed mio grid; sketch pairs use r_cut [Å].
double pair_rcut_A(const DftbElementRow& ea, const DftbElementRow& eb) {
    if (ea.z == 1 && eb.z == 1)
        return internal::MioHHTable{}.cutoff();
    return std::max(ea.r_cut, eb.r_cut);
}

/// Pair repulsive energy [eV]: mio-1-1 spline for H–H, sketch exponential else.
double pair_repulsive_eV(int za, int zb, double r_A) {
    if (za == 1 && zb == 1)
        return internal::MioHHTable{}.repulsive(r_A);
    const auto& rep = pair_rep(za, zb);
    return rep.a * std::exp(-rep.b * r_A);
}

void fill_sk_block(const DftbElementRow& ea, const DftbElementRow& eb,
                   const DftbPairSkRow& sk, double dx, double dy, double dz,
                   double r, linalg::RealMatrix& H0, linalg::RealMatrix& S,
                   int ia0, int ib0) {
    if (r < 1e-12)
        return;

    auto set_hs = [&](int i, int j, double hv, double sv) {
        H0(ia0 + i, ib0 + j) = hv;
        H0(ib0 + j, ia0 + i) = hv;
        S(ia0 + i, ib0 + j) = sv;
        S(ib0 + j, ia0 + i) = sv;
    };

    // mio-1-1 H–H (s only): shared transcribed table, not the sketch radial.
    if (ea.z == 1 && eb.z == 1) {
        const internal::MioHHTable mio;
        set_hs(0, 0, mio.hopping_ss(r), mio.overlap_ss(r));
        (void)sk;
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

/// One converged DFTB3 solve at a frozen geometry.
struct SccSolution {
    double e_h0 = 0;  ///< Tr(ρH⁰) [eV]
    double e2 = 0;    ///< E^(2) [eV]
    double e3 = 0;    ///< E^(3) [eV]
    double e_rep = 0; ///< pair repulsive [eV]
    double residual = 0;
    int iterations = 0;
    std::vector<double> dq;

    /// @brief Total electronic + repulsive energy [eV].
    double energy() const { return e_h0 + e2 + e3 + e_rep; }
};

/// The one DFTB3 self-consistency loop — energy and its finite-difference
/// gradient both go through here, so they cannot disagree on the functional.
///
/// Charges are Anderson-mixed (`internal::ChargeMixer`): undamped Picard is a
/// period-2 limit cycle for polar systems and never meets `scc_tol`. Not
/// converging is an error, never a silently returned iterate.
SccSolution solve_scc(const cpu::Frame& f, const int* type,
                      const StorageReal* x, const StorageReal* y,
                      const StorageReal* z,
                      const std::vector<const DftbElementRow*>& el,
                      const std::vector<int>& ao0, int n_at, int n_ao,
                      int n_occ, double scc_tol, int max_scc, double mix_beta) {
    auto dist_vec = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = static_cast<double>(x[j] - x[i]);
        dy = static_cast<double>(y[j] - y[i]);
        dz = static_cast<double>(z[j] - z[i]);
        f.box.min_image(dx, dy, dz);
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    linalg::RealMatrix H0 = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix S = xt::zeros<double>({n_ao, n_ao});
    // Row-major n×n minimum-image distances of this geometry [Å].
    std::vector<double> dist(static_cast<std::size_t>(n_at)
                                 * static_cast<std::size_t>(n_at),
                             0.0);
    const std::size_t n = static_cast<std::size_t>(n_at);

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
            dist[static_cast<std::size_t>(a) * n + static_cast<std::size_t>(b)] =
                r;
            dist[static_cast<std::size_t>(b) * n + static_cast<std::size_t>(a)] =
                r;
            if (r > pair_rcut_A(ea, eb))
                continue;
            fill_sk_block(ea, eb, pair_sk(type[a], type[b]), dx, dy, dz, r, H0,
                          S, ia, ao0[static_cast<std::size_t>(b)]);
        }
    }

    std::vector<double> hubbard_u(n);
    std::vector<double> hubbard_du(n);
    for (int a = 0; a < n_at; ++a) {
        hubbard_u[static_cast<std::size_t>(a)] =
            el[static_cast<std::size_t>(a)]->hubbard_u;
        hubbard_du[static_cast<std::size_t>(a)] =
            el[static_cast<std::size_t>(a)]->hubbard_du;
    }
    const internal::Dftb3Scc scc(hubbard_u, hubbard_du, dist);

    std::vector<double> dq(n, 0.0);
    std::vector<double> dq_out(n, 0.0);
    linalg::RealMatrix P = xt::zeros<double>({n_ao, n_ao});
    SccSolution out;
    bool converged = false;
    // Retained iterates in the Anderson charge mixer (dimensionless).
    constexpr int kSccMixDepth = 6;
    internal::ChargeMixer mixer(mix_beta, kSccMixDepth);

    for (int it = 0; it < max_scc; ++it) {
        const std::vector<double> dV = scc.potential(dq);

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
        const auto& C = sol.second;

        P = xt::zeros<double>({n_ao, n_ao});
        for (int o = 0; o < n_occ; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu)
                    P(mu, nu) += 2.0 * C(mu, o) * C(nu, o);

        const linalg::RealMatrix PS = xt::linalg::dot(P, S);
        double residual = 0;
        for (int a = 0; a < n_at; ++a) {
            double qa = 0;
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int i = 0; i < na; ++i)
                qa += PS(ia + i, ia + i);
            dq_out[static_cast<std::size_t>(a)] =
                qa - el[static_cast<std::size_t>(a)]->n_val;
            residual = std::max(residual,
                                std::abs(dq_out[static_cast<std::size_t>(a)]
                                         - dq[static_cast<std::size_t>(a)]));
        }

        out.residual = residual;
        out.iterations = it + 1;
        if (residual < scc_tol) {
            dq = dq_out;
            converged = true;
            break;
        }
        dq = mixer.next(dq, dq_out);
    }

    if (!converged)
        throw std::runtime_error(
            "DFTB3: SCC did not converge in " + std::to_string(max_scc)
            + " iterations (residual " + std::to_string(out.residual)
            + " e > scc_tol " + std::to_string(scc_tol)
            + "); raise max_scc or lower mix_beta");

    for (int mu = 0; mu < n_ao; ++mu)
        for (int nu = 0; nu < n_ao; ++nu)
            out.e_h0 += P(mu, nu) * H0(mu, nu);
    out.e2 = scc.second_order_energy(dq);
    out.e3 = scc.third_order_energy(dq);
    for (int a = 0; a < n_at; ++a)
        for (int b = a + 1; b < n_at; ++b)
            out.e_rep += pair_repulsive_eV(
                type[a], type[b],
                dist[static_cast<std::size_t>(a) * n
                     + static_cast<std::size_t>(b)]);
    out.dq = std::move(dq);
    return out;
}

} // namespace

void DFTB3::build(const BuildContext& /*ctx*/) {}

void DFTB3::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n_at = f.nlocal();
    if (n_at < 1)
        throw std::invalid_argument("DFTB3: empty frame");
    if (!(mix_beta > 0.0) || mix_beta > 1.0)
        throw std::invalid_argument("DFTB3: mix_beta must be in (0, 1]");

    const auto type = f.span(kAtom_type);
    StorageReal* x = f.span(kAtom_x);
    StorageReal* y = f.span(kAtom_y);
    StorageReal* z = f.span(kAtom_z);
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
        throw std::runtime_error("DFTB3: no basis functions");
    if (std::fabs(n_elec - std::round(n_elec)) > 1e-9
        || static_cast<int>(std::lround(n_elec)) % 2 != 0)
        throw std::runtime_error(
            "DFTB3: odd/non-integer valence electron count not supported");
    const int n_occ = static_cast<int>(std::lround(n_elec)) / 2;

    const SccSolution sol = solve_scc(f, type.data(), x, y, z, el, ao0, n_at,
                                      n_ao, n_occ, scc_tol, max_scc, mix_beta);
    e2_ = sol.e2;
    e3_ = sol.e3;
    scc_residual_ = sol.residual;
    scc_iterations_ = sol.iterations;
    // E = Tr(ρH⁰) + E² + E³ + E_rep. The band energy is NOT the starting
    // point: Tr(ρH) = Tr(ρH⁰) + Σ_A V_A q_A carries the GROSS Mulliken
    // population, so subtracting a double-counting built from Δq alone
    // leaves a spurious Σ_A V_A q⁰_A that is linear in Δq.
    energy_ = sol.energy();
    energy_meta += energy_;

    if (with_forces) {
        constexpr double h_fd = 1.0e-4; // Å
        for (int a = 0; a < n_at; ++a) {
            for (int k = 0; k < 3; ++k) {
                StorageReal* col = (k == 0) ? x : (k == 1) ? y : z;
                const StorageReal saved = col[a];
                col[a] = saved + static_cast<StorageReal>(h_fd);
                const double ep =
                    solve_scc(f, type.data(), x, y, z, el, ao0, n_at, n_ao,
                              n_occ, scc_tol, max_scc, mix_beta)
                        .energy();
                col[a] = saved - static_cast<StorageReal>(h_fd);
                const double em =
                    solve_scc(f, type.data(), x, y, z, el, ao0, n_at, n_ao,
                              n_occ, scc_tol, max_scc, mix_beta)
                        .energy();
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
