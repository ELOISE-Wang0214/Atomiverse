/// @file cpu/scf/dftd4.cpp
/// @brief Molecular D4::energy — Frame (Å) → PBE-D4 energy [Hartree].
///
/// Two-body DFT-D4 with EEQ2019 charges, EN-weighted erf coordination numbers,
/// and Becke–Johnson damping. Algorithm and tables match dftd4 (Caldeweyher
/// et al., J. Chem. Phys. 150, 154122 (2019)); three-body ATM is off.

#include "atomiverse/cpu/scf/dftd4.h"

#include <algorithm>
#include <array>
#include <cmath>
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

constexpr double kPi = 3.14159265358979323846;
constexpr double kSqrt2Pi = 0.79788456080286535588; // sqrt(2/pi)
constexpr double kThoPi = 3.0 / kPi;
constexpr double kKcn = 7.5;
constexpr double kCnMaxEeq = 8.0;
constexpr double kK4 = 4.10451;
constexpr double kK5 = 19.08857;
constexpr double kK6 = 2.0 * 11.28174 * 11.28174;
constexpr double kCnCutoff = 25.0;

// Casimir–Polder quadrature (dftd4 trapzd weights).
constexpr std::array<double, kD4NFreq> kFreq = {
    0.000001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.2,      1.4,  1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0};

std::array<double, kD4NFreq> trapzd_weights() {
    std::array<double, kD4NFreq> w{};
    w[0] = 0.5 * (kFreq[1] - kFreq[0]);
    for (int i = 1; i < kD4NFreq - 1; ++i)
        w[static_cast<std::size_t>(i)] =
            0.5 * ((kFreq[static_cast<std::size_t>(i)] -
                    kFreq[static_cast<std::size_t>(i - 1)]) +
                   (kFreq[static_cast<std::size_t>(i + 1)] -
                    kFreq[static_cast<std::size_t>(i)]));
    w[kD4NFreq - 1] = 0.5 * (kFreq[kD4NFreq - 1] - kFreq[kD4NFreq - 2]);
    return w;
}

const std::array<double, kD4NFreq> kTrapW = trapzd_weights();

double trapzd(const std::array<double, kD4NFreq>& pol) {
    double s = 0.0;
    for (int k = 0; k < kD4NFreq; ++k)
        s += pol[static_cast<std::size_t>(k)] * kTrapW[static_cast<std::size_t>(k)];
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
    return std::log(1.0 + std::exp(cnmax)) -
           std::log(1.0 + std::exp(cnmax - cn));
}

void require_z(int z) {
    if (z < 1 || z > kD4MaxElem || kD4NRef[static_cast<std::size_t>(z)] <= 0)
        throw std::runtime_error(
            "D4: element Z=" + std::to_string(z) +
            " is not in the in-tree D4 table (no fallback)");
}

int nref_of(int z) { return kD4NRef[static_cast<std::size_t>(z)]; }

void get_ngw(int z, std::vector<int>& ngw) {
    const int nref = nref_of(z);
    constexpr int max_cn = 19;
    int cnc[max_cn + 1] = {1};
    for (int ir = 1; ir <= nref; ++ir) {
        const int icn = std::min(
            static_cast<int>(std::lround(
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

// Reference polarizabilities for element z: alpha[ir][k], ir=1..nref.
void set_refalpha_eeq(int z, double ga, double gc,
                      std::vector<std::array<double, kD4NFreq>>& alpha) {
    const int nref = nref_of(z);
    alpha.assign(static_cast<std::size_t>(nref + 1), {});
    for (int ir = 1; ir <= nref; ++ir) {
        const int isys =
            kD4RefSys[static_cast<std::size_t>(d4_idx2(ir, z))];
        if (isys == 0)
            continue;
        const double iz = kD4Zeff[static_cast<std::size_t>(isys)];
        const double hard = kD4Hardness[static_cast<std::size_t>(isys)];
        const double clh =
            kD4Clsh[static_cast<std::size_t>(d4_idx2(ir, z))];
        const double scale_h =
            kD4Hcount[static_cast<std::size_t>(d4_idx2(ir, z))];
        const double asc =
            kD4Ascale[static_cast<std::size_t>(d4_idx2(ir, z))];
        const double zf = zeta(ga, hard * gc, iz, clh + iz);
        const double ss = kD4Sscale[static_cast<std::size_t>(isys)];
        for (int k = 0; k < kD4NFreq; ++k) {
            const double sec =
                kD4SecAiw[static_cast<std::size_t>(d4_idx_sec(isys, k))];
            const double aiw = ss * sec * zf;
            const double raw = kD4Alphaiw[static_cast<std::size_t>(
                d4_idx_aiw(k, ir, z))];
            alpha[static_cast<std::size_t>(ir)][static_cast<std::size_t>(k)] =
                std::max(asc * (raw - scale_h * aiw), 0.0);
        }
    }
}

void cn_d4(const std::vector<int>& Z,
           const std::vector<std::array<double, 3>>& pos,
           std::vector<double>& cn) {
    const int n = static_cast<int>(Z.size());
    cn.assign(static_cast<std::size_t>(n), 0.0);
    const double cut2 = kCnCutoff * kCnCutoff;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double dx = pos[static_cast<std::size_t>(i)][0] -
                              pos[static_cast<std::size_t>(j)][0];
            const double dy = pos[static_cast<std::size_t>(i)][1] -
                              pos[static_cast<std::size_t>(j)][1];
            const double dz = pos[static_cast<std::size_t>(i)][2] -
                              pos[static_cast<std::size_t>(j)][2];
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > cut2 || r2 < 1.0e-24)
                continue;
            const double r = std::sqrt(r2);
            const int zi = Z[static_cast<std::size_t>(i)];
            const int zj = Z[static_cast<std::size_t>(j)];
            const double rc = kD4Rcov[static_cast<std::size_t>(zi)] +
                              kD4Rcov[static_cast<std::size_t>(zj)];
            const double den =
                en_factor(kD4En[static_cast<std::size_t>(zi)],
                         kD4En[static_cast<std::size_t>(zj)]);
            const double c = den * erf_count(r, rc, kKcn);
            cn[static_cast<std::size_t>(i)] += c;
            if (i != j)
                cn[static_cast<std::size_t>(j)] += c;
        }
    }
}

void cn_eeq(const std::vector<int>& Z,
            const std::vector<std::array<double, 3>>& pos,
            std::vector<double>& cn) {
    const int n = static_cast<int>(Z.size());
    cn.assign(static_cast<std::size_t>(n), 0.0);
    const double cut2 = kCnCutoff * kCnCutoff;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double dx = pos[static_cast<std::size_t>(i)][0] -
                              pos[static_cast<std::size_t>(j)][0];
            const double dy = pos[static_cast<std::size_t>(i)][1] -
                              pos[static_cast<std::size_t>(j)][1];
            const double dz = pos[static_cast<std::size_t>(i)][2] -
                              pos[static_cast<std::size_t>(j)][2];
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > cut2 || r2 < 1.0e-24)
                continue;
            const double r = std::sqrt(r2);
            const int zi = Z[static_cast<std::size_t>(i)];
            const int zj = Z[static_cast<std::size_t>(j)];
            const double rc = kD4Rcov[static_cast<std::size_t>(zi)] +
                              kD4Rcov[static_cast<std::size_t>(zj)];
            const double c = erf_count(r, rc, kKcn);
            cn[static_cast<std::size_t>(i)] += c;
            if (i != j)
                cn[static_cast<std::size_t>(j)] += c;
        }
    }
    for (int i = 0; i < n; ++i)
        cn[static_cast<std::size_t>(i)] =
            log_cn_cut(cn[static_cast<std::size_t>(i)], kCnMaxEeq);
}

void eeq_charges(const std::vector<int>& Z,
                 const std::vector<std::array<double, 3>>& pos, double charge,
                 std::vector<double>& q) {
    const int n = static_cast<int>(Z.size());
    std::vector<double> cn;
    cn_eeq(Z, pos, cn);
    const int ndim = n + 1;
    std::vector<double> A(static_cast<std::size_t>(ndim * ndim), 0.0);
    std::vector<double> x(static_cast<std::size_t>(ndim), 0.0);
    auto at = [&](int i, int j) -> double& {
        return A[static_cast<std::size_t>(i * ndim + j)];
    };
    for (int i = 0; i < n; ++i) {
        const int zi = Z[static_cast<std::size_t>(i)];
        const double tmp =
            kD4EeqKcn[static_cast<std::size_t>(zi)] /
            std::sqrt(cn[static_cast<std::size_t>(i)] + 1.0e-14);
        x[static_cast<std::size_t>(i)] =
            -kD4EeqChi[static_cast<std::size_t>(zi)] +
            tmp * cn[static_cast<std::size_t>(i)];
        for (int j = 0; j < i; ++j) {
            const double dx = pos[static_cast<std::size_t>(j)][0] -
                              pos[static_cast<std::size_t>(i)][0];
            const double dy = pos[static_cast<std::size_t>(j)][1] -
                              pos[static_cast<std::size_t>(i)][1];
            const double dz = pos[static_cast<std::size_t>(j)][2] -
                              pos[static_cast<std::size_t>(i)][2];
            const double r2 = dx * dx + dy * dy + dz * dz;
            const int zj = Z[static_cast<std::size_t>(j)];
            const double gam =
                1.0 / (kD4EeqRad[static_cast<std::size_t>(zi)] *
                           kD4EeqRad[static_cast<std::size_t>(zi)] +
                       kD4EeqRad[static_cast<std::size_t>(zj)] *
                           kD4EeqRad[static_cast<std::size_t>(zj)]);
            const double r = std::sqrt(r2);
            const double tmpv = std::erf(std::sqrt(r2 * gam)) / r;
            at(i, j) = tmpv;
            at(j, i) = tmpv;
        }
        at(i, i) = kD4EeqEta[static_cast<std::size_t>(zi)] +
                   kSqrt2Pi / kD4EeqRad[static_cast<std::size_t>(zi)];
    }
    for (int i = 0; i < n; ++i) {
        at(n, i) = 1.0;
        at(i, n) = 1.0;
    }
    at(n, n) = 0.0;
    x[static_cast<std::size_t>(n)] = charge;

    // Gaussian elimination with partial pivoting (n+1 is tiny for molecules).
    std::vector<double> aug(static_cast<std::size_t>(ndim * (ndim + 1)));
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j)
            aug[static_cast<std::size_t>(i * (ndim + 1) + j)] = at(i, j);
        aug[static_cast<std::size_t>(i * (ndim + 1) + ndim)] =
            x[static_cast<std::size_t>(i)];
    }
    for (int col = 0; col < ndim; ++col) {
        int piv = col;
        double best = std::abs(
            aug[static_cast<std::size_t>(col * (ndim + 1) + col)]);
        for (int r = col + 1; r < ndim; ++r) {
            const double v = std::abs(
                aug[static_cast<std::size_t>(r * (ndim + 1) + col)]);
            if (v > best) {
                best = v;
                piv = r;
            }
        }
        if (best < 1.0e-18)
            throw std::runtime_error("D4: singular EEQ matrix");
        if (piv != col) {
            for (int c = col; c <= ndim; ++c)
                std::swap(
                    aug[static_cast<std::size_t>(col * (ndim + 1) + c)],
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
            s -= aug[static_cast<std::size_t>(i * (ndim + 1) + j)] *
                 sol[static_cast<std::size_t>(j)];
        sol[static_cast<std::size_t>(i)] =
            s / aug[static_cast<std::size_t>(i * (ndim + 1) + i)];
    }
    for (int i = 0; i < n; ++i)
        q[static_cast<std::size_t>(i)] = sol[static_cast<std::size_t>(i)];
}

void weight_refs(const std::vector<int>& Z, const std::vector<double>& cn,
                 const std::vector<double>& q, double ga, double gc,
                 double wf, std::vector<std::vector<double>>& gw) {
    const int n = static_cast<int>(Z.size());
    int mref = 0;
    for (int zi : Z)
        mref = std::max(mref, nref_of(zi));
    gw.assign(static_cast<std::size_t>(n),
              std::vector<double>(static_cast<std::size_t>(mref + 1), 0.0));
    constexpr double eps_norm = 1.0e-150; // ~ sqrt(tiny)

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
            for (int igw = 1; igw <= ngw[static_cast<std::size_t>(ir)];
                 ++igw)
                norm += weight_cn(static_cast<double>(igw) * wf, cni, cnr);
        }
        if (std::abs(norm) > eps_norm)
            norm = 1.0 / norm;
        else
            norm = 0.0;

        double maxcn = 0.0;
        for (int ir = 1; ir <= nref; ++ir)
            maxcn = std::max(
                maxcn,
                kD4RefCovCn[static_cast<std::size_t>(d4_idx2(ir, zi))]);

        for (int ir = 1; ir <= nref; ++ir) {
            const double cnr =
                kD4RefCovCn[static_cast<std::size_t>(d4_idx2(ir, zi))];
            double expw = 0.0;
            for (int igw = 1; igw <= ngw[static_cast<std::size_t>(ir)];
                 ++igw)
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
               const std::vector<std::vector<double>>& gw, double ga,
               double gc, std::vector<std::vector<double>>& c6) {
    const int n = static_cast<int>(Z.size());
    // Unique elements → reference polarizabilities.
    std::vector<int> uniq = Z;
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

    // ref_alpha[zi_slot][ir][k]
    std::vector<std::vector<std::array<double, kD4NFreq>>> ref_alpha(
        uniq.size());
    for (std::size_t u = 0; u < uniq.size(); ++u)
        set_refalpha_eeq(uniq[u], ga, gc, ref_alpha[u]);

    auto slot = [&](int z) -> std::size_t {
        return static_cast<std::size_t>(
            std::find(uniq.begin(), uniq.end(), z) - uniq.begin());
    };

    // refc6[(zi,ir,zj,jr)] via trapzd of product polarizabilities
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
                              [static_cast<std::size_t>(k)] *
                            aj[static_cast<std::size_t>(jr)]
                              [static_cast<std::size_t>(k)];
                    const double rc6 = kThoPi * trapzd(pol);
                    s += gw[static_cast<std::size_t>(i)]
                           [static_cast<std::size_t>(ir)] *
                         gw[static_cast<std::size_t>(j)]
                           [static_cast<std::size_t>(jr)] *
                         rc6;
                }
            }
            c6[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = s;
        }
    }
}

double disp2(const std::vector<int>& Z,
             const std::vector<std::array<double, 3>>& pos,
             const std::vector<std::vector<double>>& c6, double s6, double s8,
             double a1, double a2) {
    const int n = static_cast<int>(Z.size());
    double e = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            const double dx = pos[static_cast<std::size_t>(i)][0] -
                              pos[static_cast<std::size_t>(j)][0];
            const double dy = pos[static_cast<std::size_t>(i)][1] -
                              pos[static_cast<std::size_t>(j)][1];
            const double dz = pos[static_cast<std::size_t>(i)][2] -
                              pos[static_cast<std::size_t>(j)][2];
            const double r2 = dx * dx + dy * dy + dz * dz;
            const int zi = Z[static_cast<std::size_t>(i)];
            const int zj = Z[static_cast<std::size_t>(j)];
            const double rrij = 3.0 * kD4R4r2[static_cast<std::size_t>(zi)] *
                                kD4R4r2[static_cast<std::size_t>(zj)];
            const double r0ij = a1 * std::sqrt(rrij) + a2;
            const double t6 = 1.0 / (r2 * r2 * r2 + std::pow(r0ij, 6));
            const double t8 = 1.0 / (r2 * r2 * r2 * r2 + std::pow(r0ij, 8));
            const double edisp = s6 * t6 + s8 * rrij * t8;
            e += -c6[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
                 edisp;
        }
    }
    return e;
}

} // namespace

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
        pos[static_cast<std::size_t>(i)] = {x[i] * a2b, y[i] * a2b,
                                            z[i] * a2b};
    }

    std::vector<double> cn_d, q;
    cn_d4(Z, pos, cn_d);
    eeq_charges(Z, pos, charge, q);

    std::vector<std::vector<double>> gw;
    weight_refs(Z, cn_d, q, ga, gc, wf, gw);

    std::vector<std::vector<double>> c6;
    atomic_c6(Z, gw, ga, gc, c6);

    return disp2(Z, pos, c6, s6, s8, a1, a2);
}

} // namespace scf
} // namespace cpu
} // namespace ATV
