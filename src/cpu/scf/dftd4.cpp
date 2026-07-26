/// @file cpu/scf/dftd4.cpp
/// @brief DFT-D4 two-body + optional ATM three-body energy/forces
///        (molecular + orthorhombic periodic).
///
/// Algorithm matches simple-dftd4 / multicharge (Caldeweyher et al.): EN-weighted
/// erf CNs, EEQ2019 charges (Ewald A-matrix when periodic), BJ-damped pair sum
/// over lattice images. ATM (`s9`) follows dftd4 `damping/atm.f90` with neutral
/// C6 reweight (dftd4 `disp.f90`).

#include "atomiverse/cpu/scf/dftd4.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/internal/dftd4_tables.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace {

using internal::d4_idx2;
using internal::d4_idx_aiw;
using internal::d4_idx_sec;
using internal::kD4Alphaiw;
using internal::kD4Ascale;
using internal::kD4Clsh;
using internal::kD4Clsq;
using internal::kD4EeqChi;
using internal::kD4EeqEta;
using internal::kD4EeqKcn;
using internal::kD4EeqRad;
using internal::kD4En;
using internal::kD4Hardness;
using internal::kD4Hcount;
using internal::kD4MaxElem;
using internal::kD4MaxRef;
using internal::kD4NFreq;
using internal::kD4NRef;
using internal::kD4R4r2;
using internal::kD4Rcov;
using internal::kD4RefCn;
using internal::kD4RefCovCn;
using internal::kD4RefSys;
using internal::kD4SecAiw;
using internal::kD4Sscale;
using internal::kD4Zeff;

constexpr double kSqrtPi = 1.7724538509055160273;
constexpr double kSqrt2Pi = 0.79788456080286535588; // sqrt(2/pi)
constexpr double kThoPi = 3.0 / units::pi;
constexpr double kKcn = 7.5;
constexpr double kCnMaxEeq = 8.0;
constexpr double kK4 = 4.10451;
constexpr double kK5 = 19.08857;
constexpr double kK6 = 2.0 * 11.28174 * 11.28174;

// Casimir–Polder quadrature (dftd4 trapzd weights).
constexpr std::array<double, kD4NFreq> kFreq = {
    0.000001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.2,      1.4,  1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0};

std::array<double, kD4NFreq> trapzd_weights() {
    std::array<double, kD4NFreq> w{};
    w[0] = 0.5 * (kFreq[1] - kFreq[0]);
    for (int i = 1; i < kD4NFreq - 1; ++i)
        w[static_cast<std::size_t>(i)] =
            0.5
            * ((kFreq[static_cast<std::size_t>(i)]
                - kFreq[static_cast<std::size_t>(i - 1)])
               + (kFreq[static_cast<std::size_t>(i + 1)]
                  - kFreq[static_cast<std::size_t>(i)]));
    w[kD4NFreq - 1] = 0.5 * (kFreq[kD4NFreq - 1] - kFreq[kD4NFreq - 2]);
    return w;
}

const std::array<double, kD4NFreq> kTrapW = trapzd_weights();

double trapzd(const std::array<double, kD4NFreq>& pol) {
    double s = 0.0;
    for (int k = 0; k < kD4NFreq; ++k)
        s += pol[static_cast<std::size_t>(k)]
             * kTrapW[static_cast<std::size_t>(k)];
    return s;
}

double zeta(double a, double c, double qref, double qmod) {
    if (qmod < 0.0)
        return std::exp(a);
    return std::exp(a * (1.0 - std::exp(c * (1.0 - qref / qmod))));
}

double weight_cn(double wf, double cn, double cnref) {
    const double d = cn - cnref;
    return std::exp(-wf * d * d);
}

double erf_count(double r, double rc, double kcn) {
    return 0.5 * (1.0 + std::erf(-kcn * (r - rc) / rc));
}

double en_factor(double eni, double enj) {
    const double d = std::abs(eni - enj) + kK5;
    return kK4 * std::exp(-d * d / kK6);
}

double log_cn_cut(double cn, double cnmax) {
    return std::log(1.0 + std::exp(cnmax))
           - std::log(1.0 + std::exp(cnmax - cn));
}

void require_z(int z) {
    if (z < 1 || z > kD4MaxElem || kD4NRef[static_cast<std::size_t>(z)] <= 0)
        throw std::runtime_error(
            "D4: element Z=" + std::to_string(z)
            + " is not in the in-tree D4 table (no fallback)");
}

int nref_of(int z) { return kD4NRef[static_cast<std::size_t>(z)]; }

void get_ngw(int z, std::vector<int>& ngw) {
    const int nref = nref_of(z);
    constexpr int max_cn = 19;
    int cnc[max_cn + 1] = {1};
    for (int ir = 1; ir <= nref; ++ir) {
        const int icn =
            std::min(static_cast<int>(std::lround(
                         kD4RefCn[static_cast<std::size_t>(d4_idx2(ir, z))])),
                     max_cn);
        cnc[icn] += 1;
    }
    ngw.assign(static_cast<std::size_t>(nref + 1), 1);
    for (int ir = 1; ir <= nref; ++ir) {
        const int icn = cnc[std::min(
            static_cast<int>(std::lround(
                kD4RefCn[static_cast<std::size_t>(d4_idx2(ir, z))])),
            max_cn)];
        ngw[static_cast<std::size_t>(ir)] = icn * (icn + 1) / 2;
    }
}

void set_refalpha_eeq(int z, double ga, double gc,
                      std::vector<std::array<double, kD4NFreq>>& alpha) {
    const int nref = nref_of(z);
    alpha.assign(static_cast<std::size_t>(nref + 1), {});
    for (int ir = 1; ir <= nref; ++ir) {
        const int isys = kD4RefSys[static_cast<std::size_t>(d4_idx2(ir, z))];
        if (isys == 0)
            continue;
        const double iz = kD4Zeff[static_cast<std::size_t>(isys)];
        const double hard = kD4Hardness[static_cast<std::size_t>(isys)];
        const double clh = kD4Clsh[static_cast<std::size_t>(d4_idx2(ir, z))];
        const double scale_h =
            kD4Hcount[static_cast<std::size_t>(d4_idx2(ir, z))];
        const double asc = kD4Ascale[static_cast<std::size_t>(d4_idx2(ir, z))];
        const double zf = zeta(ga, hard * gc, iz, clh + iz);
        const double ss = kD4Sscale[static_cast<std::size_t>(isys)];
        for (int k = 0; k < kD4NFreq; ++k) {
            const double sec =
                kD4SecAiw[static_cast<std::size_t>(d4_idx_sec(isys, k))];
            const double aiw = ss * sec * zf;
            const double raw =
                kD4Alphaiw[static_cast<std::size_t>(d4_idx_aiw(k, ir, z))];
            alpha[static_cast<std::size_t>(ir)][static_cast<std::size_t>(k)] =
                std::max(asc * (raw - scale_h * aiw), 0.0);
        }
    }
}

// Lattice translations within rthr for an orthorhombic cell (include origin).
void lattice_points(const double cell[3], double rthr, bool periodic,
                    std::vector<std::array<double, 3>>& trans) {
    trans.clear();
    if (!periodic) {
        trans.push_back({0.0, 0.0, 0.0});
        return;
    }
    const int rx = static_cast<int>(std::ceil(rthr / cell[0]));
    const int ry = static_cast<int>(std::ceil(rthr / cell[1]));
    const int rz = static_cast<int>(std::ceil(rthr / cell[2]));
    for (int ix = -rx; ix <= rx; ++ix)
        for (int iy = -ry; iy <= ry; ++iy)
            for (int iz = -rz; iz <= rz; ++iz)
                trans.push_back({ix * cell[0], iy * cell[1], iz * cell[2]});
}

// Fixed small real/reciprocal shells for Ewald A-matrix (multicharge rep=2).
void ewald_dir_trans(const double cell[3],
                     std::vector<std::array<double, 3>>& trans) {
    trans.clear();
    for (int ix = -2; ix <= 2; ++ix)
        for (int iy = -2; iy <= 2; ++iy)
            for (int iz = -2; iz <= 2; ++iz)
                trans.push_back({ix * cell[0], iy * cell[1], iz * cell[2]});
}

void ewald_rec_trans(const double cell[3],
                     std::vector<std::array<double, 3>>& gtrans) {
    gtrans.clear();
    const double g0 = 2.0 * units::pi / cell[0];
    const double g1 = 2.0 * units::pi / cell[1];
    const double g2 = 2.0 * units::pi / cell[2];
    for (int ix = -2; ix <= 2; ++ix)
        for (int iy = -2; iy <= 2; ++iy)
            for (int iz = -2; iz <= 2; ++iz) {
                if (ix == 0 && iy == 0 && iz == 0)
                    continue;
                gtrans.push_back({ix * g0, iy * g1, iz * g2});
            }
}

// Multicharge ewald.f90 helpers for optimal α.
double get_rec_term_3d(double gg, double alpha, double vol) {
    return 4.0 * units::pi
           * (std::exp(-0.25 * gg * gg / (alpha * alpha)) / (vol * gg * gg));
}
double get_dir_term(double rr, double alpha) {
    return std::erfc(alpha * rr) / rr;
}
double rec_dir_diff(double alpha, double rlen, double dlen, double vol) {
    return (get_rec_term_3d(4.0 * rlen, alpha, vol)
            - get_rec_term_3d(5.0 * rlen, alpha, vol))
           - (get_dir_term(2.0 * dlen, alpha)
              - get_dir_term(3.0 * dlen, alpha));
}

// Multicharge search_alpha (orthorhombic).
double ewald_alpha(const double cell[3]) {
    const double vol = cell[0] * cell[1] * cell[2];
    const double dlen = std::min({cell[0], cell[1], cell[2]});
    // Shortest reciprocal vector length: 2π / longest real edge.
    const double rlen = 2.0 * units::pi / std::max({cell[0], cell[1], cell[2]});
    const double tol = std::sqrt(std::numeric_limits<double>::epsilon());
    double alpha = 1.0e-8;
    double diff = rec_dir_diff(alpha, rlen, dlen, vol);
    while (diff < -tol && alpha < 1.0e300) {
        alpha *= 2.0;
        diff = rec_dir_diff(alpha, rlen, dlen, vol);
    }
    if (!(alpha > 1.0e-8) || alpha >= 1.0e300)
        return 0.25;
    double alpl = 0.5 * alpha;
    while (diff < tol && alpha < 1.0e300) {
        alpha *= 2.0;
        diff = rec_dir_diff(alpha, rlen, dlen, vol);
    }
    if (alpha >= 1.0e300)
        return 0.25;
    double alpr = alpha;
    alpha = 0.5 * (alpl + alpr);
    for (int it = 0; it < 30; ++it) {
        diff = rec_dir_diff(alpha, rlen, dlen, vol);
        if (std::abs(diff) <= tol)
            break;
        if (diff < 0.0)
            alpl = alpha;
        else
            alpr = alpha;
        alpha = 0.5 * (alpl + alpr);
    }
    return alpha;
}

void cn_d4(const std::vector<int>& Z,
           const std::vector<std::array<double, 3>>& pos, const double cell[3],
           bool periodic, double cn_cut, std::vector<double>& cn) {
    const int n = static_cast<int>(Z.size());
    cn.assign(static_cast<std::size_t>(n), 0.0);
    std::vector<std::array<double, 3>> trans;
    lattice_points(cell, cn_cut, periodic, trans);
    const double cut2 = cn_cut * cn_cut;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const int zi = Z[static_cast<std::size_t>(i)];
            const int zj = Z[static_cast<std::size_t>(j)];
            const double rc = kD4Rcov[static_cast<std::size_t>(zi)]
                              + kD4Rcov[static_cast<std::size_t>(zj)];
            const double den = en_factor(kD4En[static_cast<std::size_t>(zi)],
                                         kD4En[static_cast<std::size_t>(zj)]);
            for (const auto& t : trans) {
                const double dx =
                    pos[static_cast<std::size_t>(i)][0]
                    - (pos[static_cast<std::size_t>(j)][0] + t[0]);
                const double dy =
                    pos[static_cast<std::size_t>(i)][1]
                    - (pos[static_cast<std::size_t>(j)][1] + t[1]);
                const double dz =
                    pos[static_cast<std::size_t>(i)][2]
                    - (pos[static_cast<std::size_t>(j)][2] + t[2]);
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > cut2 || r2 < 1.0e-12)
                    continue;
                const double r = std::sqrt(r2);
                const double c = den * erf_count(r, rc, kKcn);
                cn[static_cast<std::size_t>(i)] += c;
                if (i != j)
                    cn[static_cast<std::size_t>(j)] += c;
            }
        }
    }
}

void cn_eeq(const std::vector<int>& Z,
            const std::vector<std::array<double, 3>>& pos, const double cell[3],
            bool periodic, double cn_cut, std::vector<double>& cn) {
    const int n = static_cast<int>(Z.size());
    cn.assign(static_cast<std::size_t>(n), 0.0);
    std::vector<std::array<double, 3>> trans;
    lattice_points(cell, cn_cut, periodic, trans);
    const double cut2 = cn_cut * cn_cut;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const int zi = Z[static_cast<std::size_t>(i)];
            const int zj = Z[static_cast<std::size_t>(j)];
            const double rc = kD4Rcov[static_cast<std::size_t>(zi)]
                              + kD4Rcov[static_cast<std::size_t>(zj)];
            for (const auto& t : trans) {
                const double dx =
                    pos[static_cast<std::size_t>(i)][0]
                    - (pos[static_cast<std::size_t>(j)][0] + t[0]);
                const double dy =
                    pos[static_cast<std::size_t>(i)][1]
                    - (pos[static_cast<std::size_t>(j)][1] + t[1]);
                const double dz =
                    pos[static_cast<std::size_t>(i)][2]
                    - (pos[static_cast<std::size_t>(j)][2] + t[2]);
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > cut2 || r2 < 1.0e-12)
                    continue;
                const double r = std::sqrt(r2);
                const double c = erf_count(r, rc, kKcn);
                cn[static_cast<std::size_t>(i)] += c;
                if (i != j)
                    cn[static_cast<std::size_t>(j)] += c;
            }
        }
    }
    for (int i = 0; i < n; ++i)
        cn[static_cast<std::size_t>(i)] =
            log_cn_cut(cn[static_cast<std::size_t>(i)], kCnMaxEeq);
}

double amat_dir(const double rij[3], double gam, double alp,
                const std::vector<std::array<double, 3>>& dtrans) {
    double s = 0.0;
    constexpr double eps = 1.0e-12;
    for (const auto& t : dtrans) {
        const double vx = rij[0] + t[0];
        const double vy = rij[1] + t[1];
        const double vz = rij[2] + t[2];
        const double r1 = std::sqrt(vx * vx + vy * vy + vz * vz);
        if (r1 < eps)
            continue;
        s += std::erf(gam * r1) / r1 - std::erf(alp * r1) / r1;
    }
    return s;
}

double amat_rec(const double rij[3], double vol, double alp,
                const std::vector<std::array<double, 3>>& gtrans) {
    double s = 0.0;
    const double fac = 4.0 * units::pi / vol;
    constexpr double eps = 1.0e-12;
    for (const auto& g : gtrans) {
        const double g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        if (g2 < eps)
            continue;
        const double phase = g[0] * rij[0] + g[1] * rij[1] + g[2] * rij[2];
        s += std::cos(phase) * fac * std::exp(-0.25 * g2 / (alp * alp)) / g2;
    }
    return s;
}

void eeq_charges(const std::vector<int>& Z,
                 const std::vector<std::array<double, 3>>& pos,
                 const double cell[3], bool periodic, double charge,
                 std::vector<double>& q) {
    const int n = static_cast<int>(Z.size());
    // EEQ CN uses fixed 25 Bohr cutoff (dftd4 charge.f90).
    constexpr double kEeqCnCut = 25.0;
    std::vector<double> cn;
    cn_eeq(Z, pos, cell, periodic, kEeqCnCut, cn);

    const int ndim = n + 1;
    std::vector<double> A(static_cast<std::size_t>(ndim * ndim), 0.0);
    std::vector<double> x(static_cast<std::size_t>(ndim), 0.0);
    auto at = [&](int i, int j) -> double& {
        return A[static_cast<std::size_t>(i * ndim + j)];
    };

    if (!periodic) {
        for (int i = 0; i < n; ++i) {
            const int zi = Z[static_cast<std::size_t>(i)];
            const double tmp =
                kD4EeqKcn[static_cast<std::size_t>(zi)]
                / std::sqrt(cn[static_cast<std::size_t>(i)] + 1.0e-14);
            x[static_cast<std::size_t>(i)] =
                -kD4EeqChi[static_cast<std::size_t>(zi)]
                + tmp * cn[static_cast<std::size_t>(i)];
            for (int j = 0; j < i; ++j) {
                const double dx = pos[static_cast<std::size_t>(j)][0]
                                  - pos[static_cast<std::size_t>(i)][0];
                const double dy = pos[static_cast<std::size_t>(j)][1]
                                  - pos[static_cast<std::size_t>(i)][1];
                const double dz = pos[static_cast<std::size_t>(j)][2]
                                  - pos[static_cast<std::size_t>(i)][2];
                const double r2 = dx * dx + dy * dy + dz * dz;
                const int zj = Z[static_cast<std::size_t>(j)];
                const double gam =
                    1.0
                    / (kD4EeqRad[static_cast<std::size_t>(zi)]
                           * kD4EeqRad[static_cast<std::size_t>(zi)]
                       + kD4EeqRad[static_cast<std::size_t>(zj)]
                             * kD4EeqRad[static_cast<std::size_t>(zj)]);
                const double r = std::sqrt(r2);
                const double tmpv = std::erf(std::sqrt(r2 * gam)) / r;
                at(i, j) = tmpv;
                at(j, i) = tmpv;
            }
            at(i, i) = kD4EeqEta[static_cast<std::size_t>(zi)]
                       + kSqrt2Pi / kD4EeqRad[static_cast<std::size_t>(zi)];
        }
    } else {
        std::vector<std::array<double, 3>> dtrans, gtrans;
        ewald_dir_trans(cell, dtrans);
        ewald_rec_trans(cell, gtrans);
        const double alp = ewald_alpha(cell);
        const double vol = cell[0] * cell[1] * cell[2];

        for (int i = 0; i < n; ++i) {
            const int zi = Z[static_cast<std::size_t>(i)];
            const double tmp =
                kD4EeqKcn[static_cast<std::size_t>(zi)]
                / std::sqrt(cn[static_cast<std::size_t>(i)] + 1.0e-14);
            x[static_cast<std::size_t>(i)] =
                -kD4EeqChi[static_cast<std::size_t>(zi)]
                + tmp * cn[static_cast<std::size_t>(i)];

            for (int j = 0; j < i; ++j) {
                const int zj = Z[static_cast<std::size_t>(j)];
                // multicharge: gam = 1/sqrt(rad_i^2 + rad_j^2)
                const double gam =
                    1.0
                    / std::sqrt(
                        kD4EeqRad[static_cast<std::size_t>(zi)]
                            * kD4EeqRad[static_cast<std::size_t>(zi)]
                        + kD4EeqRad[static_cast<std::size_t>(zj)]
                              * kD4EeqRad[static_cast<std::size_t>(zj)]);
                const double rij[3] = {
                    pos[static_cast<std::size_t>(i)][0]
                        - pos[static_cast<std::size_t>(j)][0],
                    pos[static_cast<std::size_t>(i)][1]
                        - pos[static_cast<std::size_t>(j)][1],
                    pos[static_cast<std::size_t>(i)][2]
                        - pos[static_cast<std::size_t>(j)][2]};
                const double aij = amat_dir(rij, gam, alp, dtrans)
                                   + amat_rec(rij, vol, alp, gtrans);
                at(i, j) = aij;
                at(j, i) = aij;
            }
            // Self term: gam = 1/sqrt(2 rad^2)
            const double gam_ii =
                1.0
                / std::sqrt(2.0 * kD4EeqRad[static_cast<std::size_t>(zi)]
                            * kD4EeqRad[static_cast<std::size_t>(zi)]);
            const double rij0[3] = {0.0, 0.0, 0.0};
            const double a_self = amat_dir(rij0, gam_ii, alp, dtrans)
                                  + amat_rec(rij0, vol, alp, gtrans);
            at(i, i) = a_self + kD4EeqEta[static_cast<std::size_t>(zi)]
                       + kSqrt2Pi / kD4EeqRad[static_cast<std::size_t>(zi)]
                       - 2.0 * alp / kSqrtPi;
        }
    }

    for (int i = 0; i < n; ++i) {
        at(n, i) = 1.0;
        at(i, n) = 1.0;
    }
    at(n, n) = 0.0;
    x[static_cast<std::size_t>(n)] = charge;

    // Gaussian elimination with partial pivoting.
    std::vector<double> aug(static_cast<std::size_t>(ndim * (ndim + 1)));
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j)
            aug[static_cast<std::size_t>(i * (ndim + 1) + j)] = at(i, j);
        aug[static_cast<std::size_t>(i * (ndim + 1) + ndim)] =
            x[static_cast<std::size_t>(i)];
    }
    for (int col = 0; col < ndim; ++col) {
        int piv = col;
        double best =
            std::abs(aug[static_cast<std::size_t>(col * (ndim + 1) + col)]);
        for (int r = col + 1; r < ndim; ++r) {
            const double v =
                std::abs(aug[static_cast<std::size_t>(r * (ndim + 1) + col)]);
            if (v > best) {
                best = v;
                piv = r;
            }
        }
        if (best < 1.0e-18)
            throw std::runtime_error("D4: singular EEQ matrix");
        if (piv != col) {
            for (int c = col; c <= ndim; ++c)
                std::swap(aug[static_cast<std::size_t>(col * (ndim + 1) + c)],
                          aug[static_cast<std::size_t>(piv * (ndim + 1) + c)]);
        }
        const double diag =
            aug[static_cast<std::size_t>(col * (ndim + 1) + col)];
        for (int r = col + 1; r < ndim; ++r) {
            const double f =
                aug[static_cast<std::size_t>(r * (ndim + 1) + col)] / diag;
            for (int c = col; c <= ndim; ++c)
                aug[static_cast<std::size_t>(r * (ndim + 1) + c)] -=
                    f * aug[static_cast<std::size_t>(col * (ndim + 1) + c)];
        }
    }
    q.resize(static_cast<std::size_t>(n));
    std::vector<double> sol(static_cast<std::size_t>(ndim));
    for (int i = ndim - 1; i >= 0; --i) {
        double s = aug[static_cast<std::size_t>(i * (ndim + 1) + ndim)];
        for (int j = i + 1; j < ndim; ++j)
            s -= aug[static_cast<std::size_t>(i * (ndim + 1) + j)]
                 * sol[static_cast<std::size_t>(j)];
        sol[static_cast<std::size_t>(i)] =
            s / aug[static_cast<std::size_t>(i * (ndim + 1) + i)];
    }
    for (int i = 0; i < n; ++i)
        q[static_cast<std::size_t>(i)] = sol[static_cast<std::size_t>(i)];
}

void weight_refs(const std::vector<int>& Z, const std::vector<double>& cn,
                 const std::vector<double>& q, double ga, double gc, double wf,
                 std::vector<std::vector<double>>& gw) {
    const int n = static_cast<int>(Z.size());
    int mref = 0;
    for (int zi : Z)
        mref = std::max(mref, nref_of(zi));
    gw.assign(static_cast<std::size_t>(n),
              std::vector<double>(static_cast<std::size_t>(mref + 1), 0.0));
    constexpr double eps_norm = 1.0e-150;

    for (int iat = 0; iat < n; ++iat) {
        const int zi = Z[static_cast<std::size_t>(iat)];
        const int nref = nref_of(zi);
        std::vector<int> ngw;
        get_ngw(zi, ngw);
        const double zi_eff = kD4Zeff[static_cast<std::size_t>(zi)];
        const double gi = kD4Hardness[static_cast<std::size_t>(zi)] * gc;
        const double cni = cn[static_cast<std::size_t>(iat)];
        const double qi = q[static_cast<std::size_t>(iat)];

        double norm = 0.0;
        for (int ir = 1; ir <= nref; ++ir) {
            const double cnr =
                kD4RefCovCn[static_cast<std::size_t>(d4_idx2(ir, zi))];
            for (int igw = 1; igw <= ngw[static_cast<std::size_t>(ir)]; ++igw)
                norm += weight_cn(static_cast<double>(igw) * wf, cni, cnr);
        }
        if (std::abs(norm) > eps_norm)
            norm = 1.0 / norm;
        else
            norm = 0.0;

        double maxcn = 0.0;
        for (int ir = 1; ir <= nref; ++ir)
            maxcn = std::max(
                maxcn, kD4RefCovCn[static_cast<std::size_t>(d4_idx2(ir, zi))]);

        for (int ir = 1; ir <= nref; ++ir) {
            const double cnr =
                kD4RefCovCn[static_cast<std::size_t>(d4_idx2(ir, zi))];
            double expw = 0.0;
            for (int igw = 1; igw <= ngw[static_cast<std::size_t>(ir)]; ++igw)
                expw += weight_cn(static_cast<double>(igw) * wf, cni, cnr);
            double gwk = expw * norm;
            if (!std::isfinite(gwk) || norm == 0.0) {
                gwk = (std::abs(maxcn - cnr) < 1.0e-12) ? 1.0 : 0.0;
            }
            const double qref =
                kD4Clsq[static_cast<std::size_t>(d4_idx2(ir, zi))];
            gw[static_cast<std::size_t>(iat)][static_cast<std::size_t>(ir)] =
                gwk * zeta(ga, gi, qref + zi_eff, qi + zi_eff);
        }
    }
}

void atomic_c6(const std::vector<int>& Z,
               const std::vector<std::vector<double>>& gw, double ga, double gc,
               std::vector<std::vector<double>>& c6) {
    const int n = static_cast<int>(Z.size());
    std::vector<int> uniq = Z;
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

    std::vector<std::vector<std::array<double, kD4NFreq>>> ref_alpha(
        uniq.size());
    for (std::size_t u = 0; u < uniq.size(); ++u)
        set_refalpha_eeq(uniq[u], ga, gc, ref_alpha[u]);

    auto slot = [&](int z) -> std::size_t {
        return static_cast<std::size_t>(std::find(uniq.begin(), uniq.end(), z)
                                        - uniq.begin());
    };

    c6.assign(static_cast<std::size_t>(n),
              std::vector<double>(static_cast<std::size_t>(n), 0.0));
    for (int i = 0; i < n; ++i) {
        const int zi = Z[static_cast<std::size_t>(i)];
        const int nrefi = nref_of(zi);
        const auto& ai = ref_alpha[slot(zi)];
        for (int j = 0; j < n; ++j) {
            const int zj = Z[static_cast<std::size_t>(j)];
            const int nrefj = nref_of(zj);
            const auto& aj = ref_alpha[slot(zj)];
            double s = 0.0;
            for (int ir = 1; ir <= nrefi; ++ir) {
                for (int jr = 1; jr <= nrefj; ++jr) {
                    std::array<double, kD4NFreq> pol{};
                    for (int k = 0; k < kD4NFreq; ++k)
                        pol[static_cast<std::size_t>(k)] =
                            ai[static_cast<std::size_t>(ir)]
                              [static_cast<std::size_t>(k)]
                            * aj[static_cast<std::size_t>(jr)]
                                [static_cast<std::size_t>(k)];
                    const double rc6 = kThoPi * trapzd(pol);
                    s += gw[static_cast<std::size_t>(i)]
                           [static_cast<std::size_t>(ir)]
                         * gw[static_cast<std::size_t>(j)]
                             [static_cast<std::size_t>(jr)]
                         * rc6;
                }
            }
            c6[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = s;
        }
    }
}

double disp2(const std::vector<int>& Z,
             const std::vector<std::array<double, 3>>& pos,
             const double cell[3], bool periodic,
             const std::vector<std::vector<double>>& c6, double s6, double s8,
             double a1, double a2, double disp2_cut) {
    const int n = static_cast<int>(Z.size());
    std::vector<std::array<double, 3>> trans;
    lattice_points(cell, disp2_cut, periodic, trans);
    const double cut2 = disp2_cut * disp2_cut;
    double e = 0.0;
    // Matches simple-dftd4 rational damping: loop iat, jat<=iat, all T;
    // partition 0.5 to each endpoint so self-image pairs get half weight.
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const int zi = Z[static_cast<std::size_t>(i)];
            const int zj = Z[static_cast<std::size_t>(j)];
            const double rrij = 3.0 * kD4R4r2[static_cast<std::size_t>(zi)]
                                * kD4R4r2[static_cast<std::size_t>(zj)];
            const double r0ij = a1 * std::sqrt(rrij) + a2;
            const double c6ij =
                c6[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
            for (const auto& t : trans) {
                const double dx =
                    pos[static_cast<std::size_t>(i)][0]
                    - (pos[static_cast<std::size_t>(j)][0] + t[0]);
                const double dy =
                    pos[static_cast<std::size_t>(i)][1]
                    - (pos[static_cast<std::size_t>(j)][1] + t[1]);
                const double dz =
                    pos[static_cast<std::size_t>(i)][2]
                    - (pos[static_cast<std::size_t>(j)][2] + t[2]);
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > cut2 || r2 < 1.0e-24)
                    continue;
                const double t6 = 1.0 / (r2 * r2 * r2 + std::pow(r0ij, 6));
                const double t8 = 1.0 / (r2 * r2 * r2 * r2 + std::pow(r0ij, 8));
                const double edisp = s6 * t6 + s8 * rrij * t8;
                // Total energy contribution for this (pair, T):
                // i!=j: full -c6*edisp; i==j: half (dftd4 0.5 partition).
                if (i != j)
                    e += -c6ij * edisp;
                else
                    e += -0.5 * c6ij * edisp;
            }
        }
    }
    return e;
}

// dftd4 damping/atm.f90 triple_scale — distribute image-self triples.
double triple_scale(int ii, int jj, int kk) {
    if (ii == jj) {
        if (ii == kk)
            return 1.0 / 6.0; // ii'i"
        return 0.5;           // ii'j
    }
    if (ii != kk && jj != kk)
        return 1.0; // ijk
    return 0.5;     // ijj' / iji'
}

/// @brief Axilrod–Teller–Muto three-body (dftd4 get_atm_dispersion_energy).
///
/// C6 matrix must be the **neutral-charge** reweight (q=0), matching
/// simple-dftd4 `disp.f90`. Zero-damping alp=16 with BJ critical radii.
double disp3_atm(const std::vector<int>& Z,
                 const std::vector<std::array<double, 3>>& pos,
                 const double cell[3], bool periodic,
                 const std::vector<std::vector<double>>& c6, double s9,
                 double a1, double a2, double alp, double disp3_cut) {
    if (std::abs(s9) < 1.0e-30)
        return 0.0;
    const int n = static_cast<int>(Z.size());
    if (n < 1)
        return 0.0;

    std::vector<std::array<double, 3>> trans;
    lattice_points(cell, disp3_cut, periodic, trans);
    const double cutoff2 = disp3_cut * disp3_cut;
    constexpr double eps = 1.0e-24;
    double energy = 0.0;

    for (int iat = 0; iat < n; ++iat) {
        const int iz = Z[static_cast<std::size_t>(iat)];
        for (int jat = 0; jat <= iat; ++jat) {
            const int jz = Z[static_cast<std::size_t>(jat)];
            const double c6ij =
                c6[static_cast<std::size_t>(jat)][static_cast<std::size_t>(iat)];
            const double r0ij =
                a1
                    * std::sqrt(3.0 * kD4R4r2[static_cast<std::size_t>(jz)]
                                * kD4R4r2[static_cast<std::size_t>(iz)])
                + a2;
            for (const auto& tj : trans) {
                const double vij0 =
                    pos[static_cast<std::size_t>(jat)][0] + tj[0]
                    - pos[static_cast<std::size_t>(iat)][0];
                const double vij1 =
                    pos[static_cast<std::size_t>(jat)][1] + tj[1]
                    - pos[static_cast<std::size_t>(iat)][1];
                const double vij2 =
                    pos[static_cast<std::size_t>(jat)][2] + tj[2]
                    - pos[static_cast<std::size_t>(iat)][2];
                const double r2ij = vij0 * vij0 + vij1 * vij1 + vij2 * vij2;
                if (r2ij > cutoff2 || r2ij < eps)
                    continue;
                for (int kat = 0; kat <= jat; ++kat) {
                    const int kz = Z[static_cast<std::size_t>(kat)];
                    const double c6ik =
                        c6[static_cast<std::size_t>(kat)]
                          [static_cast<std::size_t>(iat)];
                    const double c6jk =
                        c6[static_cast<std::size_t>(kat)]
                          [static_cast<std::size_t>(jat)];
                    const double c9 =
                        -s9 * std::sqrt(std::abs(c6ij * c6ik * c6jk));
                    const double r0ik =
                        a1
                            * std::sqrt(3.0 * kD4R4r2[static_cast<std::size_t>(kz)]
                                        * kD4R4r2[static_cast<std::size_t>(iz)])
                        + a2;
                    const double r0jk =
                        a1
                            * std::sqrt(3.0 * kD4R4r2[static_cast<std::size_t>(kz)]
                                        * kD4R4r2[static_cast<std::size_t>(jz)])
                        + a2;
                    const double r0 = r0ij * r0ik * r0jk;
                    const double triple = triple_scale(iat, jat, kat);
                    for (const auto& tk : trans) {
                        const double vik0 =
                            pos[static_cast<std::size_t>(kat)][0] + tk[0]
                            - pos[static_cast<std::size_t>(iat)][0];
                        const double vik1 =
                            pos[static_cast<std::size_t>(kat)][1] + tk[1]
                            - pos[static_cast<std::size_t>(iat)][1];
                        const double vik2 =
                            pos[static_cast<std::size_t>(kat)][2] + tk[2]
                            - pos[static_cast<std::size_t>(iat)][2];
                        const double r2ik =
                            vik0 * vik0 + vik1 * vik1 + vik2 * vik2;
                        if (r2ik > cutoff2 || r2ik < eps)
                            continue;
                        // vjk = (kat + tk) - (jat + tj)
                        const double vjk0 =
                            pos[static_cast<std::size_t>(kat)][0] + tk[0]
                            - (pos[static_cast<std::size_t>(jat)][0] + tj[0]);
                        const double vjk1 =
                            pos[static_cast<std::size_t>(kat)][1] + tk[1]
                            - (pos[static_cast<std::size_t>(jat)][1] + tj[1]);
                        const double vjk2 =
                            pos[static_cast<std::size_t>(kat)][2] + tk[2]
                            - (pos[static_cast<std::size_t>(jat)][2] + tj[2]);
                        const double r2jk =
                            vjk0 * vjk0 + vjk1 * vjk1 + vjk2 * vjk2;
                        if (r2jk > cutoff2 || r2jk < eps)
                            continue;

                        const double r2 = r2ij * r2ik * r2jk;
                        const double r1 = std::sqrt(r2);
                        const double r3 = r2 * r1;
                        const double r5 = r3 * r2;

                        const double fdmp =
                            1.0
                            / (1.0
                               + 6.0 * std::pow(r0 / r1, alp / 3.0));
                        const double ang =
                            0.375
                                * (r2ij + r2jk - r2ik) * (r2ij - r2jk + r2ik)
                                * (-r2ij + r2jk + r2ik) / r5
                            + 1.0 / r3;
                        const double rr = ang * fdmp;
                        // Per-atom: energy_i -= dE/3 with dE = rr*c9*triple
                        // → total E += -rr*c9*triple
                        energy += -rr * c9 * triple;
                    }
                }
            }
        }
    }
    return energy;
}

} // namespace

double dftd4_energy(const std::vector<int>& Z,
                    const std::vector<std::array<double, 3>>& pos_bohr,
                    const double cell_bohr[3], bool periodic, double s6,
                    double s8, double a1, double a2, double ga, double gc,
                    double wf, double charge, double cn_cut, double disp2_cut,
                    double s9, double disp3_cut) {
    const int n = static_cast<int>(Z.size());
    if (n <= 0)
        return 0.0;
    if (static_cast<int>(pos_bohr.size()) != n)
        throw std::runtime_error("D4: Z and pos size mismatch");
    for (int z : Z)
        require_z(z);
    if (periodic
        && (cell_bohr[0] <= 0.0 || cell_bohr[1] <= 0.0 || cell_bohr[2] <= 0.0))
        throw std::runtime_error("D4: periodic cell edges must be positive");

    const double cell[3] = {cell_bohr[0], cell_bohr[1], cell_bohr[2]};

    std::vector<double> cn_d, q;
    cn_d4(Z, pos_bohr, cell, periodic, cn_cut, cn_d);
    eeq_charges(Z, pos_bohr, cell, periodic, charge, q);

    // Two-body: EEQ-charged C6 (dftd4 default).
    std::vector<std::vector<double>> gw;
    weight_refs(Z, cn_d, q, ga, gc, wf, gw);
    std::vector<std::vector<double>> c6;
    atomic_c6(Z, gw, ga, gc, c6);
    double e = disp2(Z, pos_bohr, cell, periodic, c6, s6, s8, a1, a2, disp2_cut);

    // ATM three-body: reweight C6 at q = 0 (dftd4 disp.f90).
    if (std::abs(s9) > 1.0e-30) {
        std::vector<double> q0(static_cast<std::size_t>(n), 0.0);
        std::vector<std::vector<double>> gw0;
        weight_refs(Z, cn_d, q0, ga, gc, wf, gw0);
        std::vector<std::vector<double>> c6_atm;
        atomic_c6(Z, gw0, ga, gc, c6_atm);
        constexpr double kAlpAtm = 16.0;
        e += disp3_atm(Z, pos_bohr, cell, periodic, c6_atm, s9, a1, a2, kAlpAtm,
                       disp3_cut);
    }
    return e;
}

D4Result dftd4(const std::vector<int>& Z,
               const std::vector<std::array<double, 3>>& pos_bohr,
               const double cell_bohr[3], bool periodic, double s6, double s8,
               double a1, double a2, double ga, double gc, double wf,
               double charge, double cn_cut, double disp2_cut, double s9,
               double disp3_cut, double grad_h) {
    const int n = static_cast<int>(Z.size());
    D4Result res;
    res.forces.assign(static_cast<std::size_t>(n), {0.0, 0.0, 0.0});
    if (n <= 0)
        return res;
    if (!(grad_h > 0.0) || !std::isfinite(grad_h))
        throw std::runtime_error("D4: grad_h must be positive and finite");

    res.energy = dftd4_energy(Z, pos_bohr, cell_bohr, periodic, s6, s8, a1, a2,
                             ga, gc, wf, charge, cn_cut, disp2_cut, s9,
                             disp3_cut);

    // Exact ∇E of the geometry-only energy functional via central differences
    // of dftd4_energy (CN + EEQ + BJ pair + ATM all move with nuclei).
    std::vector<std::array<double, 3>> pos = pos_bohr;
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < 3; ++k) {
            const double x0 = pos[static_cast<std::size_t>(i)]
                                  [static_cast<std::size_t>(k)];
            pos[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                x0 + grad_h;
            const double ep =
                dftd4_energy(Z, pos, cell_bohr, periodic, s6, s8, a1, a2, ga,
                             gc, wf, charge, cn_cut, disp2_cut, s9, disp3_cut);
            pos[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                x0 - grad_h;
            const double em =
                dftd4_energy(Z, pos, cell_bohr, periodic, s6, s8, a1, a2, ga,
                             gc, wf, charge, cn_cut, disp2_cut, s9, disp3_cut);
            pos[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = x0;
            // F = −dE/dR
            res.forces[static_cast<std::size_t>(i)]
                      [static_cast<std::size_t>(k)] = -(ep - em) / (2.0 * grad_h);
        }
    }
    return res;
}

double D4::energy(const Frame& frame) const {
    const int n = static_cast<int>(frame.nlocal());
    if (n <= 0)
        return 0.0;

    const auto type = frame.span(kAtom_type);
    const auto x = frame.span(kAtom_x);
    const auto y = frame.span(kAtom_y);
    const auto z = frame.span(kAtom_z);

    std::vector<int> Z(static_cast<std::size_t>(n));
    std::vector<std::array<double, 3>> pos(static_cast<std::size_t>(n));
    const double a2b = units::kAngstromInBohr;
    for (int i = 0; i < n; ++i) {
        const int zi = static_cast<int>(type[i]);
        require_z(zi);
        Z[static_cast<std::size_t>(i)] = zi;
        pos[static_cast<std::size_t>(i)] = {x[i] * a2b, y[i] * a2b, z[i] * a2b};
    }
    const double cell[3] = {0.0, 0.0, 0.0};
    return dftd4_energy(Z, pos, cell, /*periodic=*/false, s6, s8, a1, a2, ga,
                        gc, wf, charge, /*cn_cut=*/30.0, /*disp2_cut=*/60.0,
                        s9);
}

} // namespace scf
} // namespace cpu
} // namespace ATV
