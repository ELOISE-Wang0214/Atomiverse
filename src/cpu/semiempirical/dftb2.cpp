/// @file src/cpu/semiempirical/dftb2.cpp
/// SCC-DFTB2 multi-element: mio-1-1 H–H pin + sketch (spin + analytic F) + 3ob.
/// Γ-point PBC (MioSketch): lattice-summed H0/S, screened γ, and repulsive.
#include "atomiverse/cpu/semiempirical/dftb2.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_3ob_tables.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_mio_hh_sk.h"
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
using basis_data::kAToBohr;
using basis_data::kBohrToA;
using basis_data::kHaToEV;
using basis_data::kMioHH_Hss_Ha;
using basis_data::kMioHH_Sss;
using basis_data::kMioHH_distFudge_Bohr;
using basis_data::kMioHH_expCoeffs;
using basis_data::kMioHH_gridDist_Bohr;
using basis_data::kMioHH_nGrid;
using basis_data::kMioHH_nSpline;
using basis_data::kMioHH_repCutoff_Bohr;
using basis_data::kMioHH_spCoeffs;
using basis_data::kMioHH_spLastCoeffs;
using basis_data::kMioHH_xStart_Bohr;

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

/// Neville polynomial interpolation on a uniform abscissa (DFTB+ polyInterU1).
double poly_interp(const double* xp, const double* yp, int nn, double xx) {
    std::vector<double> X(static_cast<std::size_t>(nn + 1));
    std::vector<double> C(static_cast<std::size_t>(nn + 1));
    std::vector<double> D(static_cast<std::size_t>(nn + 1));
    for (int i = 1; i <= nn; ++i) {
        X[static_cast<std::size_t>(i)] = xp[i - 1];
        C[static_cast<std::size_t>(i)] = yp[i - 1];
        D[static_cast<std::size_t>(i)] = yp[i - 1];
    }
    int iCl = static_cast<int>(
        std::ceil((xx - X[1]) / std::abs(X[2] - X[1])));
    if (iCl < 1)
        iCl = 1;
    if (iCl > nn)
        iCl = nn;
    double yy = yp[iCl - 1];
    iCl = iCl - 1;
    for (int mm = 1; mm <= nn - 1; ++mm) {
        for (int ii = 1; ii <= nn - mm; ++ii) {
            double rTmp = X[static_cast<std::size_t>(ii)]
                          - X[static_cast<std::size_t>(ii + mm)];
            rTmp = (C[static_cast<std::size_t>(ii + 1)]
                    - D[static_cast<std::size_t>(ii)])
                   / rTmp;
            C[static_cast<std::size_t>(ii)] =
                (X[static_cast<std::size_t>(ii)] - xx) * rTmp;
            D[static_cast<std::size_t>(ii)] =
                (X[static_cast<std::size_t>(ii + mm)] - xx) * rTmp;
        }
        double dyy;
        if (2 * iCl < nn - mm) {
            dyy = C[static_cast<std::size_t>(iCl + 1)];
        } else {
            dyy = D[static_cast<std::size_t>(iCl)];
            iCl = iCl - 1;
        }
        yy += dyy;
    }
    return yy;
}

/// DFTB+ “new” SK grid interpolation for one integral column (Ha or dimless).
double mio_hh_interp(const double* tab, double r_bohr) {
    constexpr int nInterNew = 8;
    constexpr int nRightInterNew = 4;
    constexpr double deltaR = 1.0e-5;
    const int leng = static_cast<int>(kMioHH_nGrid);
    const double incr = kMioHH_gridDist_Bohr;
    const double rMax =
        static_cast<double>(leng) * incr + kMioHH_distFudge_Bohr;
    const int ind = static_cast<int>(std::floor(r_bohr / incr));
    if (r_bohr >= rMax)
        return 0.0;
    if (ind < leng) {
        int iLast = std::min(leng, ind + nRightInterNew);
        iLast = std::max(iLast, nInterNew);
        double xa[nInterNew];
        double ya[nInterNew];
        for (int ii = 1; ii <= nInterNew; ++ii) {
            xa[ii - 1] =
                static_cast<double>(iLast - nInterNew + ii) * incr;
            ya[ii - 1] = tab[iLast - nInterNew + ii - 1];
        }
        return poly_interp(xa, ya, nInterNew, r_bohr);
    }
    // Beyond last grid point: 5th-order taper to zero over distFudge.
    const double dr = r_bohr - rMax;
    const int iLast = leng;
    double xa[nInterNew];
    double ya[nInterNew];
    for (int ii = 1; ii <= nInterNew; ++ii) {
        xa[ii - 1] = static_cast<double>(iLast - nInterNew + ii) * incr;
        ya[ii - 1] = tab[iLast - nInterNew + ii - 1];
    }
    const double y0 =
        poly_interp(xa, ya, nInterNew, xa[nInterNew - 1] - deltaR);
    const double y2 =
        poly_interp(xa, ya, nInterNew, xa[nInterNew - 1] + deltaR);
    const double y1 = ya[nInterNew - 1];
    const double y1p = (y2 - y0) / (2.0 * deltaR);
    const double y1pp = (y2 + y0 - 2.0 * y1) / (deltaR * deltaR);
    const double dx = -kMioHH_distFudge_Bohr;
    const double invdx = -1.0 / kMioHH_distFudge_Bohr;
    const double dx1 = y1p * dx;
    const double dx2 = y1pp * dx * dx;
    const double dd = 10.0 * y1 - 4.0 * dx1 + 0.5 * dx2;
    const double ee = -15.0 * y1 + 7.0 * dx1 - 1.0 * dx2;
    const double ff = 6.0 * y1 - 3.0 * dx1 + 0.5 * dx2;
    const double xr = dr * invdx;
    return ((ff * xr + ee) * xr + dd) * xr * xr * xr;
}

/// mio-1-1 H–H repulsive energy [eV] at separation r [Å].
double mio_hh_repulsive_eV(double r_A) {
    const double rr = r_A * kAToBohr;
    if (rr < 1e-12 || rr >= kMioHH_repCutoff_Bohr)
        return 0.0;
    if (rr < kMioHH_xStart_Bohr[0]) {
        const double e_ha = std::exp(-kMioHH_expCoeffs[0] * rr
                                     + kMioHH_expCoeffs[1])
                            + kMioHH_expCoeffs[2];
        return e_ha * kHaToEV;
    }
    std::size_t lo = 0;
    std::size_t hi = kMioHH_nSpline;
    while (hi - lo > 1) {
        const std::size_t mid = (lo + hi) / 2;
        if (kMioHH_xStart_Bohr[mid] <= rr)
            lo = mid;
        else
            hi = mid;
    }
    const std::size_t iMatch = lo;
    const double dr = rr - kMioHH_xStart_Bohr[iMatch];
    double e_ha = 0.0;
    if (iMatch + 1 < kMioHH_nSpline) {
        const double* c = &kMioHH_spCoeffs[4 * iMatch];
        e_ha = c[0] + c[1] * dr + c[2] * dr * dr + c[3] * dr * dr * dr;
    } else {
        double xh = 1.0;
        for (int ii = 0; ii < 6; ++ii) {
            e_ha += kMioHH_spLastCoeffs[ii] * xh;
            xh *= dr;
        }
    }
    return e_ha * kHaToEV;
}

/// dE_rep/dr [eV/Å] for mio-1-1 H–H.
double mio_hh_repulsive_dE_dr(double r_A) {
    const double rr = r_A * kAToBohr;
    if (rr < 1e-12 || rr >= kMioHH_repCutoff_Bohr)
        return 0.0;
    double de_drr_ha = 0.0;
    if (rr < kMioHH_xStart_Bohr[0]) {
        de_drr_ha = -kMioHH_expCoeffs[0]
                    * std::exp(-kMioHH_expCoeffs[0] * rr
                               + kMioHH_expCoeffs[1]);
    } else {
        std::size_t lo = 0;
        std::size_t hi = kMioHH_nSpline;
        while (hi - lo > 1) {
            const std::size_t mid = (lo + hi) / 2;
            if (kMioHH_xStart_Bohr[mid] <= rr)
                lo = mid;
            else
                hi = mid;
        }
        const std::size_t iMatch = lo;
        const double dr = rr - kMioHH_xStart_Bohr[iMatch];
        if (iMatch + 1 < kMioHH_nSpline) {
            const double* c = &kMioHH_spCoeffs[4 * iMatch];
            de_drr_ha = c[1] + 2.0 * c[2] * dr + 3.0 * c[3] * dr * dr;
        } else {
            // d/dr (c0 + c1 dr + ... + c5 dr^5)
            de_drr_ha = kMioHH_spLastCoeffs[1]
                        + 2.0 * kMioHH_spLastCoeffs[2] * dr
                        + 3.0 * kMioHH_spLastCoeffs[3] * dr * dr
                        + 4.0 * kMioHH_spLastCoeffs[4] * dr * dr * dr
                        + 5.0 * kMioHH_spLastCoeffs[5] * dr * dr * dr * dr;
        }
    }
    return de_drr_ha * kHaToEV * kAToBohr;
}

/// H–H mio cutoff [Å]; sketch pairs use element r_cut.
double pair_rcut_A(const DftbElementRow& ea, const DftbElementRow& eb, int za,
                   int zb) {
    if (za == 1 && zb == 1) {
        const double r_bohr = static_cast<double>(kMioHH_nGrid)
                                  * kMioHH_gridDist_Bohr
                              + kMioHH_distFudge_Bohr;
        return r_bohr * kBohrToA;
    }
    return std::max(ea.r_cut, eb.r_cut);
}

double pair_repulsive_eV(int za, int zb, double r_A) {
    if (za == 1 && zb == 1)
        return mio_hh_repulsive_eV(r_A);
    const auto& rep = pair_rep(za, zb);
    return rep.a * std::exp(-rep.b * r_A);
}

double pair_repulsive_dE_dr(int za, int zb, double r_A) {
    if (za == 1 && zb == 1)
        return mio_hh_repulsive_dE_dr(r_A);
    const auto& rep = pair_rep(za, zb);
    return -rep.a * rep.b * std::exp(-rep.b * r_A);
}

/// @brief Write directed A→B block of H0/S for bond vector R = R_B − R_A (accumulate).
///
/// Does **not** write the transpose — lattice-sum callers walk all ordered pairs (a,b,T).
/// Min-image molecular forces keep a separate assign-both-triangles helper below.
void fill_sk_block_directed(const DftbElementRow& ea, const DftbElementRow& eb,
                            const DftbPairSkRow& sk, double dx, double dy,
                            double dz, double r, linalg::RealMatrix& H0,
                            linalg::RealMatrix& S, int ia0, int ib0) {
    if (r < 1e-12)
        return;

    // mio-1-1 H–H (s only): table interpolation in Bohr / Hartree.
    if (ea.z == 1 && eb.z == 1) {
        const double r_bohr = r * kAToBohr;
        const double Hss = mio_hh_interp(kMioHH_Hss_Ha, r_bohr) * kHaToEV;
        const double Sss = mio_hh_interp(kMioHH_Sss, r_bohr);
        H0(ia0, ib0) += Hss;
        S(ia0, ib0) += Sss;
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

    auto add_hs = [&](int i, int j, double hv, double sv) {
        H0(ia0 + i, ib0 + j) += hv;
        S(ia0 + i, ib0 + j) += sv;
    };

    add_hs(0, 0, Hss, Sss);

    if (ea.n_ao >= 1 && eb.n_ao == 4) {
        add_hs(0, 1, l * Hsp, l * Ssp);
        add_hs(0, 2, m * Hsp, m * Ssp);
        add_hs(0, 3, n * Hsp, n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao >= 1) {
        add_hs(1, 0, -l * Hsp, -l * Ssp);
        add_hs(2, 0, -m * Hsp, -m * Ssp);
        add_hs(3, 0, -n * Hsp, -n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao == 4) {
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


struct LatVec {
    double tx, ty, tz;
};

/// @brief Enumerate lattice translations with |T| ≲ max_cut + cell edge.
/// Columns of `cell` are lattice vectors (matches Box). Includes T = 0.
std::vector<LatVec> enumerate_lattice(const Mat3& cell, double max_cut) {
    double col_len[3];
    for (int d = 0; d < 3; ++d) {
        const double cx = static_cast<double>(cell(0, d));
        const double cy = static_cast<double>(cell(1, d));
        const double cz = static_cast<double>(cell(2, d));
        col_len[d] = std::sqrt(cx * cx + cy * cy + cz * cz);
        if (col_len[d] < 1e-9)
            throw std::runtime_error("DFTB2: degenerate cell edge");
    }
    int nmax[3];
    for (int d = 0; d < 3; ++d)
        nmax[d] = std::max(
            1, static_cast<int>(std::ceil(max_cut / col_len[d])) + 1);

    std::vector<LatVec> out;
    out.reserve(static_cast<std::size_t>((2 * nmax[0] + 1) * (2 * nmax[1] + 1)
                                         * (2 * nmax[2] + 1)));
    for (int n0 = -nmax[0]; n0 <= nmax[0]; ++n0) {
        for (int n1 = -nmax[1]; n1 <= nmax[1]; ++n1) {
            for (int n2 = -nmax[2]; n2 <= nmax[2]; ++n2) {
                const double tx = n0 * static_cast<double>(cell(0, 0))
                                  + n1 * static_cast<double>(cell(0, 1))
                                  + n2 * static_cast<double>(cell(0, 2));
                const double ty = n0 * static_cast<double>(cell(1, 0))
                                  + n1 * static_cast<double>(cell(1, 1))
                                  + n2 * static_cast<double>(cell(1, 2));
                const double tz = n0 * static_cast<double>(cell(2, 0))
                                  + n1 * static_cast<double>(cell(2, 1))
                                  + n2 * static_cast<double>(cell(2, 2));
                out.push_back({tx, ty, tz});
            }
        }
    }
    return out;
}

std::vector<LatVec>
lattice_images_for(const cpu::Box& box,
                   const std::vector<const DftbElementRow*>& el) {
    double max_cut = 0;
    for (const auto* e : el)
        max_cut = std::max(max_cut, e->r_cut);
    // mio-1-1 H–H table may extend beyond sketch r_cut.
    {
        const double r_bohr = static_cast<double>(kMioHH_nGrid)
                                  * kMioHH_gridDist_Bohr
                              + kMioHH_distFudge_Bohr;
        max_cut = std::max(max_cut, r_bohr * kBohrToA);
    }
    if (max_cut < 1e-12)
        max_cut = 1.0;
    return enumerate_lattice(box.cell(), max_cut);
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

    // mio-1-1 H–H: d/dr via central FD of table interp, then chain to bond vec.
    if (ea.z == 1 && eb.z == 1) {
        constexpr double h = 1.0e-5; // Å
        const double r_bohr = r * kAToBohr;
        const double rp = (r + h) * kAToBohr;
        const double rm = std::max(r - h, 1.0e-8) * kAToBohr;
        const double dHss =
            (mio_hh_interp(kMioHH_Hss_Ha, rp) - mio_hh_interp(kMioHH_Hss_Ha, rm))
            * kHaToEV / (2.0 * h);
        const double dSss =
            (mio_hh_interp(kMioHH_Sss, rp) - mio_hh_interp(kMioHH_Sss, rm))
            / (2.0 * h);
        for (int k = 0; k < 3; ++k) {
            dH0[k](ia0, ib0) = dHss * dir[k];
            dS[k](ia0, ib0) = dSss * dir[k];
        }
        (void)sk;
        (void)r_bohr;
        return;
    }

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
    // Restricted-path orbitals for closed-form forces (empty if unrestricted).
    linalg::RealMatrix C;
    linalg::RealVector eps;
    int n_occ = 0;
};

/// Build lattice-summed H0/S at the Γ-point (molecular limit when L ≫ r_cut).
void build_h0_s(int n_at, int n_ao, const std::vector<const DftbElementRow*>& el,
                const std::vector<int>& ao0, const int* type,
                const StorageReal* x, const StorageReal* y, const StorageReal* z,
                const cpu::Box& box, linalg::RealMatrix& H0,
                linalg::RealMatrix& S) {
    H0 = xt::zeros<double>({n_ao, n_ao});
    S = xt::zeros<double>({n_ao, n_ao});
    const auto images = lattice_images_for(box, el);

    // Onsite (T = 0, a = b): diagonal e_s / e_p and S = 1.
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
    }

    // Lattice sum over ordered pairs (a,b) and images T of B.
    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        const double xa = static_cast<double>(x[a]);
        const double ya = static_cast<double>(y[a]);
        const double za = static_cast<double>(z[a]);
        for (int b = 0; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            const int ib = ao0[static_cast<std::size_t>(b)];
            const double xb = static_cast<double>(x[b]);
            const double yb = static_cast<double>(y[b]);
            const double zb = static_cast<double>(z[b]);
            const double rcut = pair_rcut_A(ea, eb, type[a], type[b]);
            const auto& sk = pair_sk(type[a], type[b]);
            for (const auto& T : images) {
                const double dx = xb + T.tx - xa;
                const double dy = yb + T.ty - ya;
                const double dz = zb + T.tz - za;
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < 1e-12 || r > rcut)
                    continue;
                fill_sk_block_directed(ea, eb, sk, dx, dy, dz, r, H0, S, ia,
                                       ib);
            }
        }
    }
}

double repulsive_energy(int n_at, const int* type, const StorageReal* x,
                        const StorageReal* y, const StorageReal* z,
                        const cpu::Box& box,
                        const std::vector<const DftbElementRow*>& el) {
    const auto images = lattice_images_for(box, el);
    double e_rep = 0;
    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const double xa = static_cast<double>(x[a]);
        const double ya = static_cast<double>(y[a]);
        const double za = static_cast<double>(z[a]);
        for (int b = 0; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            const double xb = static_cast<double>(x[b]);
            const double yb = static_cast<double>(y[b]);
            const double zb = static_cast<double>(z[b]);
            const double rcut = pair_rcut_A(ea, eb, type[a], type[b]);
            for (const auto& T : images) {
                const double dx = xb + T.tx - xa;
                const double dy = yb + T.ty - ya;
                const double dz = zb + T.tz - za;
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < 1e-12 || r > rcut)
                    continue;
                // 1/2 for double-counting over ordered (a,b,T).
                e_rep += 0.5 * pair_repulsive_eV(type[a], type[b], r);
            }
        }
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
    const double ua = el[static_cast<std::size_t>(a)]->hubbard_u;
    const double ub = el[static_cast<std::size_t>(b)]->hubbard_u;
    const double rcut = pair_rcut_A(*el[static_cast<std::size_t>(a)],
                                    *el[static_cast<std::size_t>(b)],
                                    el[static_cast<std::size_t>(a)]->z,
                                    el[static_cast<std::size_t>(b)]->z);
    const auto images = lattice_images_for(box, el);
    const double xa = static_cast<double>(x[a]);
    const double ya = static_cast<double>(y[a]);
    const double za = static_cast<double>(z[a]);
    const double xb = static_cast<double>(x[b]);
    const double yb = static_cast<double>(y[b]);
    const double zb = static_cast<double>(z[b]);
    double g = 0;
    for (const auto& T : images) {
        const double dx = xb + T.tx - xa;
        const double dy = yb + T.ty - ya;
        const double dz = zb + T.tz - za;
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (r < 1e-12) {
            if (a == b)
                g += ua; // onsite Hubbard
            continue;
        }
        if (r > rcut)
            continue;
        g += gamma_ab(r, ua, ub);
    }
    return g;
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
    out.energy = e_band - e_scc + repulsive_energy(n_at, type, x, y, z, box, el);
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
                 + repulsive_energy(n_at, type, x, y, z, box, el);
    out.dq = std::move(dq);
    out.spin = std::move(mspin);
    return out;
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

            if (r <= pair_rcut_A(ea, eb, type[a], type[b])) {
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

            const double dErep_dr = pair_repulsive_dE_dr(type[a], type[b], r);
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
    // The electronic per-term share of the total, same operator, same value
    // (esmm-contract-05): DFTB2 already accumulated into kEnergy.
    f.meta.get(kEnergyElectronic) += energy;

    if (with_forces) {
        if (!occ.unrestricted && res.n_occ > 0 && res.C.shape(0) > 0) {
            analytic_forces_restricted(n_at, n_ao, res.n_occ, el, ao0, type, x,
                                       y, z, f.box, res.dq, res.C, res.eps, fx,
                                       fy, fz);
        } else {
            // Open-shell: re-SCC central finite difference of E.
            constexpr double h_fd = 1.0e-4; // Å
            for (int a = 0; a < n_at; ++a) {
                for (int k = 0; k < 3; ++k) {
                    StorageReal* col = (k == 0) ? x : (k == 1) ? y : z;
                    const StorageReal saved = col[a];
                    col[a] = saved + static_cast<StorageReal>(h_fd);
                    const double ep =
                        scc_total(n_at, n_ao, occ, el, ao0, type, x, y, z,
                                  f.box, scc_tol, max_scc)
                            .energy;
                    col[a] = saved - static_cast<StorageReal>(h_fd);
                    const double em =
                        scc_total(n_at, n_ao, occ, el, ao0, type, x, y, z,
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
    // The electronic per-term share of the total, same operator, same value
    // (esmm-contract-05): DFTB2 already accumulated into kEnergy.
    f.meta.get(kEnergyElectronic) += energy;

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
