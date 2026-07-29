/// @file src/cpu/semiempirical/mio_sk_block.cpp
/// mio-style sketch Slater–Koster two-centre block + its bond-vector
/// derivative. Moved verbatim out of `dftb2.cpp`'s private helpers so DFTB2 and
/// DFTB3 share one kernel — the derivative must be the derivative of the value
/// that the SCC actually converged.
#include "atomiverse/cpu/semiempirical/internal/mio_sk_block.h"

#include <cmath>

#include "atomiverse/cpu/semiempirical/internal/mio_hh_table.h"

namespace ATV {
namespace cpu {
namespace semiempirical {
namespace internal {
namespace {

/// Sketch Hamiltonian radial channel: v0 exp(−(r−r0)/d) [eV].
double h_radial(double v0, double r, double r0, double d) {
    return v0 * std::exp(-(r - r0) / d);
}

/// Sketch overlap radial channel: s0 exp(−½((r−r0)/w)²) (dimensionless).
double s_radial(double s0, double r, double r0, double w) {
    const double t = (r - r0) / w;
    return s0 * std::exp(-0.5 * t * t);
}

} // namespace

void MioSketchPair::add_block(double dx, double dy, double dz, double r,
                              linalg::RealMatrix& h0, linalg::RealMatrix& s,
                              int ia0, int ib0) const {
    if (r < 1e-12)
        return;

    // mio-1-1 H–H (s only): shared transcribed table, not the sketch radial.
    if (a_.z == 1 && b_.z == 1) {
        const MioHHTable mio;
        h0(ia0, ib0) += mio.hopping_ss(r);
        s(ia0, ib0) += mio.overlap_ss(r);
        return;
    }

    const double inv_r = 1.0 / r;
    const double l = dx * inv_r;
    const double m = dy * inv_r;
    const double n = dz * inv_r;

    const double Hss = h_radial(sk_.h_ss0, r, sk_.h_r0, sk_.h_d);
    const double Hsp = h_radial(sk_.h_sp0, r, sk_.h_r0, sk_.h_d);
    const double Hpps = h_radial(sk_.h_pps0, r, sk_.h_r0, sk_.h_d);
    const double Hppp = h_radial(sk_.h_ppp0, r, sk_.h_r0, sk_.h_d);
    const double Sss = s_radial(sk_.s_ss0, r, sk_.s_r0, sk_.s_w);
    const double Ssp = s_radial(sk_.s_sp0, r, sk_.s_r0, sk_.s_w);
    const double Spps = s_radial(sk_.s_pps0, r, sk_.s_r0, sk_.s_w);
    const double Sppp = s_radial(sk_.s_ppp0, r, sk_.s_r0, sk_.s_w);

    auto add_hs = [&](int i, int j, double hv, double sv) {
        h0(ia0 + i, ib0 + j) += hv;
        s(ia0 + i, ib0 + j) += sv;
    };

    add_hs(0, 0, Hss, Sss);

    if (a_.n_ao >= 1 && b_.n_ao == 4) {
        add_hs(0, 1, l * Hsp, l * Ssp);
        add_hs(0, 2, m * Hsp, m * Ssp);
        add_hs(0, 3, n * Hsp, n * Ssp);
    }
    if (a_.n_ao == 4 && b_.n_ao >= 1) {
        add_hs(1, 0, -l * Hsp, -l * Ssp);
        add_hs(2, 0, -m * Hsp, -m * Ssp);
        add_hs(3, 0, -n * Hsp, -n * Ssp);
    }
    if (a_.n_ao == 4 && b_.n_ao == 4) {
        const double dir[3] = {l, m, n};
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                const double d_ab = (a == b) ? 1.0 : 0.0;
                const double hv = dir[a] * dir[b] * (Hpps - Hppp) + d_ab * Hppp;
                const double sv = dir[a] * dir[b] * (Spps - Sppp) + d_ab * Sppp;
                add_hs(1 + a, 1 + b, hv, sv);
            }
        }
    }
}

void MioSketchPair::set_block_derivative(double dx, double dy, double dz,
                                         double r,
                                         linalg::RealMatrix (&dh0)[3],
                                         linalg::RealMatrix (&ds)[3], int ia0,
                                         int ib0) const {
    if (r < 1e-12)
        return;
    const double inv_r = 1.0 / r;
    const double dir[3] = {dx * inv_r, dy * inv_r, dz * inv_r};

    // mio-1-1 H–H: radial derivative of the shared table, chained to bond vec.
    if (a_.z == 1 && b_.z == 1) {
        const MioHHTable mio;
        const double dHss = mio.d_hopping_ss_dr(r);
        const double dSss = mio.d_overlap_ss_dr(r);
        for (int k = 0; k < 3; ++k) {
            dh0[k](ia0, ib0) = dHss * dir[k];
            ds[k](ia0, ib0) = dSss * dir[k];
        }
        return;
    }

    const double Hss = h_radial(sk_.h_ss0, r, sk_.h_r0, sk_.h_d);
    const double Hsp = h_radial(sk_.h_sp0, r, sk_.h_r0, sk_.h_d);
    const double Hpps = h_radial(sk_.h_pps0, r, sk_.h_r0, sk_.h_d);
    const double Hppp = h_radial(sk_.h_ppp0, r, sk_.h_r0, sk_.h_d);
    const double Sss = s_radial(sk_.s_ss0, r, sk_.s_r0, sk_.s_w);
    const double Ssp = s_radial(sk_.s_sp0, r, sk_.s_r0, sk_.s_w);
    const double Spps = s_radial(sk_.s_pps0, r, sk_.s_r0, sk_.s_w);
    const double Sppp = s_radial(sk_.s_ppp0, r, sk_.s_r0, sk_.s_w);

    const double dHss = -Hss / sk_.h_d;
    const double dHsp = -Hsp / sk_.h_d;
    const double dHpps = -Hpps / sk_.h_d;
    const double dHppp = -Hppp / sk_.h_d;
    const double dSss = Sss * (-(r - sk_.s_r0) / (sk_.s_w * sk_.s_w));
    const double dSsp = Ssp * (-(r - sk_.s_r0) / (sk_.s_w * sk_.s_w));
    const double dSpps = Spps * (-(r - sk_.s_r0) / (sk_.s_w * sk_.s_w));
    const double dSppp = Sppp * (-(r - sk_.s_r0) / (sk_.s_w * sk_.s_w));

    auto ddir = [&](int alpha, int k) {
        const double dab = (alpha == k) ? 1.0 : 0.0;
        return (dab - dir[alpha] * dir[k]) * inv_r;
    };

    for (int k = 0; k < 3; ++k) {
        const double dr_dk = dir[k];
        dh0[k](ia0, ib0) = dHss * dr_dk;
        ds[k](ia0, ib0) = dSss * dr_dk;

        if (a_.n_ao >= 1 && b_.n_ao == 4) {
            for (int a = 0; a < 3; ++a) {
                dh0[k](ia0, ib0 + 1 + a) =
                    ddir(a, k) * Hsp + dir[a] * dHsp * dr_dk;
                ds[k](ia0, ib0 + 1 + a) =
                    ddir(a, k) * Ssp + dir[a] * dSsp * dr_dk;
            }
        }
        if (a_.n_ao == 4 && b_.n_ao >= 1) {
            for (int a = 0; a < 3; ++a) {
                dh0[k](ia0 + 1 + a, ib0) =
                    -(ddir(a, k) * Hsp + dir[a] * dHsp * dr_dk);
                ds[k](ia0 + 1 + a, ib0) =
                    -(ddir(a, k) * Ssp + dir[a] * dSsp * dr_dk);
            }
        }
        if (a_.n_ao == 4 && b_.n_ao == 4) {
            const double dHsig = dHpps - dHppp;
            const double dSsig = dSpps - dSppp;
            const double Hsig = Hpps - Hppp;
            const double Ssig = Spps - Sppp;
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    const double dab = (a == b) ? 1.0 : 0.0;
                    dh0[k](ia0 + 1 + a, ib0 + 1 + b) =
                        (ddir(a, k) * dir[b] + dir[a] * ddir(b, k)) * Hsig
                        + dir[a] * dir[b] * dHsig * dr_dk + dab * dHppp * dr_dk;
                    ds[k](ia0 + 1 + a, ib0 + 1 + b) =
                        (ddir(a, k) * dir[b] + dir[a] * ddir(b, k)) * Ssig
                        + dir[a] * dir[b] * dSsig * dr_dk + dab * dSppp * dr_dk;
                }
            }
        }
    }
}

} // namespace internal
} // namespace semiempirical
} // namespace cpu
} // namespace ATV
