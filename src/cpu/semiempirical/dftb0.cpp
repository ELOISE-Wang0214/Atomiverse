/// @file src/cpu/semiempirical/dftb0.cpp
/// Non-SCC DFTB0 (H/C/N/O): transcribed mio-1-1 H–H, mio-style sketch elsewhere.
#include "atomiverse/cpu/semiempirical/dftb0.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_mio_tables.h"
#include "atomiverse/cpu/semiempirical/internal/mio_hh_table.h"

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
    throw std::runtime_error("DFTB0: unsupported element Z=" + std::to_string(z)
                             + " (static mio-style tables cover H/C/N/O only)");
}

const DftbPairSkRow& pair_sk(int za, int zb) {
    for (std::size_t i = 0; i < kDftbNPairSk; ++i)
        if (kDftbPairSk[i].z_a == za && kDftbPairSk[i].z_b == zb)
            return kDftbPairSk[i];
    throw std::runtime_error("DFTB0: no SK table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

const DftbPairRepRow& pair_rep(int za, int zb) {
    for (std::size_t i = 0; i < kDftbNPairRep; ++i)
        if (kDftbPairRep[i].z_a == za && kDftbPairRep[i].z_b == zb)
            return kDftbPairRep[i];
    throw std::runtime_error("DFTB0: no repulsive table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

double h_radial(double v0, double r, double r0, double d) {
    return v0 * std::exp(-(r - r0) / d);
}
double s_radial(double s0, double r, double r0, double w) {
    const double t = (r - r0) / w;
    return s0 * std::exp(-0.5 * t * t);
}

/// @brief Fill n_a × n_b block of H0/S for atoms A–B separated by R vec.
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

/// Non-SCC total energy at the current Frame geometry [eV].
double energy_non_scc(const Frame& f, const StorageReal* x, const StorageReal* y,
                      const StorageReal* z, const int* type, int n_at,
                      const std::vector<const DftbElementRow*>& el,
                      const std::vector<int>& ao0, int n_ao, int n_occ) {
    linalg::RealMatrix H0 = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix S = xt::zeros<double>({n_ao, n_ao});

    auto dist_vec = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = static_cast<double>(x[j] - x[i]);
        dy = static_cast<double>(y[j] - y[i]);
        dz = static_cast<double>(z[j] - z[i]);
        f.box.min_image(dx, dy, dz);
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    const internal::MioHHTable mio_hh;

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
            const bool hh = (ea.z == 1 && eb.z == 1);
            double dx, dy, dz;
            const double r = dist_vec(a, b, dx, dy, dz);
            if (r > (hh ? mio_hh.cutoff() : std::max(ea.r_cut, eb.r_cut)))
                continue;
            const int ib = ao0[static_cast<std::size_t>(b)];
            if (hh) {
                // mio-1-1 H–H (s only): tabulated H0/S, no sketch radial form.
                const double Hss = mio_hh.hopping_ss(r);
                const double Sss = mio_hh.overlap_ss(r);
                H0(ia, ib) = Hss;
                H0(ib, ia) = Hss;
                S(ia, ib) = Sss;
                S(ib, ia) = Sss;
                continue;
            }
            const auto& sk = pair_sk(type[a], type[b]);
            fill_sk_block(ea, eb, sk, dx, dy, dz, r, H0, S, ia, ib);
        }
    }

    auto sol = linalg::diagonalize_fock_generalized(H0, S);
    const auto& eps = sol.first;

    double e_band = 0;
    for (int o = 0; o < n_occ; ++o)
        e_band += 2.0 * eps(o);

    double e_repulsive = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = a + 1; b < n_at; ++b) {
            double dx, dy, dz;
            const double r = dist_vec(a, b, dx, dy, dz);
            if (type[a] == 1 && type[b] == 1) {
                // mio-1-1 H–H repulsive spline (zero past its own cutoff).
                e_repulsive += mio_hh.repulsive(r);
                continue;
            }
            const auto& rep = pair_rep(type[a], type[b]);
            e_repulsive += rep.a * std::exp(-rep.b * r);
        }

    return e_band + e_repulsive;
}

} // namespace

void DFTB0::build(const BuildContext& /*ctx*/) {}

void DFTB0::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n_at = f.nlocal();
    if (n_at < 1)
        throw std::invalid_argument("DFTB0: empty frame");

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
        throw std::runtime_error("DFTB0: no basis functions");
    if (std::fabs(n_elec - std::round(n_elec)) > 1e-9
        || static_cast<int>(std::lround(n_elec)) % 2 != 0)
        throw std::runtime_error(
            "DFTB0: odd/non-integer valence electron count not supported");
    const int n_occ = static_cast<int>(std::lround(n_elec)) / 2;

    energy_ = energy_non_scc(f, x, y, z, type, n_at, el, ao0, n_ao, n_occ);
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
                    energy_non_scc(f, mx, my, mz, type, n_at, el, ao0, n_ao,
                                   n_occ);
                col[a] = saved - static_cast<StorageReal>(h_fd);
                const double em =
                    energy_non_scc(f, mx, my, mz, type, n_at, el, ao0, n_ao,
                                   n_occ);
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
