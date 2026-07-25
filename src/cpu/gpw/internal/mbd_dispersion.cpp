/// @file cpu/gpw/internal/mbd_dispersion.cpp
/// @brief MBD@rsSCS free-atom energy (isolated), matching libMBD pure-Python.

#include "atomiverse/cpu/gpw/internal/mbd_dispersion.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace gpw {
namespace {

using linalg::RealMatrix;
using linalg::RealVector;

// Gauss–Legendre n = 15 on [-1,1] (numpy.polynomial.legendre.leggauss).
constexpr int kNLeg = 15;
constexpr double kLegX[kNLeg] = {
    -0.9879925180204854, -0.9372733924007058, -0.8482065834104272,
    -0.7244177313601701, -0.5709721726085388, -0.3941513470775634,
    -0.20119409399743451, 0.0,                 0.20119409399743451,
    0.3941513470775634,  0.5709721726085388,  0.7244177313601701,
    0.8482065834104272,  0.9372733924007058,  0.9879925180204854};
constexpr double kLegW[kNLeg] = {
    0.030753241996117203, 0.0703660474881084,  0.10715922046717141,
    0.13957067792615444,  0.16626920581699398, 0.1861610000155622,
    0.1984314853271116,   0.2025782419255613,  0.1984314853271116,
    0.1861610000155622,   0.16626920581699398, 0.13957067792615444,
    0.10715922046717141,  0.0703660474881084,  0.030753241996117203};

// libMBD freq_grid(n=15, L=0.6): prepend u=0 with weight 0, reverse-mapped
// Gauss–Legendre half-line.
constexpr int kNFreq = 16; // 1 + n

void free_atom_params(int Z, double& alpha, double& c6, double& rvdw) {
    if (!ts_supports_element(Z))
        throw std::runtime_error(
            "GPW MBD: element Z=" + std::to_string(Z)
            + " is not in the in-tree free-atom table (no fallback)");
    const double* row = kTsFreeAtom[Z - 1];
    alpha = row[0];
    c6 = row[1];
    rvdw = row[2];
}

/// @brief Fermi damping f = 1 / (1 + exp(−a (R/S − 1))).
double damping_fermi(double r, double s_vdw, double a) {
    if (!(s_vdw > 0.0))
        return 0.0;
    const double x = -a * (r / s_vdw - 1.0);
    if (x > 60.0)
        return 0.0;
    if (x < -60.0)
        return 1.0;
    return 1.0 / (1.0 + std::exp(x));
}

/// @brief Bare dipole tensor T_ab = (3 r_a r_b − r² δ_ab) / r⁵  (note sign:
///        libMBD returns the negative of the electrostatic ∇∇(1/r) convention
///        used in some texts — we match libMBD T_bare).
void T_bare(double Rx, double Ry, double Rz, double T[3][3]) {
    const double r2 = Rx * Rx + Ry * Ry + Rz * Rz;
    if (r2 < 1e-30) {
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                T[a][b] = 0.0;
        return;
    }
    const double r = std::sqrt(r2);
    const double r5 = r2 * r2 * r; // r^5
    const double inv_r5 = 1.0 / r5;
    const double R[3] = {Rx, Ry, Rz};
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            // libMBD: (-3 R_a R_b + r² δ_ab) / r⁵
            T[a][b] = (-3.0 * R[a] * R[b] + (a == b ? r2 : 0.0)) * inv_r5;
        }
    }
}

/// @brief Gaussian-damped Coulomb dipole tensor (libMBD T_erf_coulomb).
void T_erf_coulomb(double Rx, double Ry, double Rz, double sigma, double T[3][3]) {
    double bare[3][3];
    T_bare(Rx, Ry, Rz, bare);
    const double r2 = Rx * Rx + Ry * Ry + Rz * Rz;
    if (r2 < 1e-30 || !(sigma > 0.0)) {
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                T[a][b] = 0.0;
        return;
    }
    const double r = std::sqrt(r2);
    const double r5 = r2 * r2 * r;
    const double inv_r5 = 1.0 / r5;
    const double zeta = r / sigma;
    const double theta =
        2.0 * zeta / std::sqrt(M_PI) * std::exp(-zeta * zeta);
    const double erf_theta = std::erf(zeta) - theta;
    const double R[3] = {Rx, Ry, Rz};
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const double rr_r5 = R[a] * R[b] * inv_r5;
            T[a][b] =
                erf_theta * bare[a][b] + (2.0 * zeta * zeta * theta) * rr_r5;
        }
    }
}

/// @brief Build 3N×3N dipole interaction matrix (isolated, no lattice).
///
/// damping:
///  - "fermi,dip,gg": short-range Gaussian (1−f)×T_erf  for rsSCS
///  - "fermi,dip":    long-range Fermi f×T_bare for MBD Hamiltonian
RealMatrix dipole_matrix_isolated(
    const std::vector<std::array<double, 3>>& pos, const std::vector<double>& R_vdw,
    const std::vector<double>& sigma, double beta, double a_fermi, bool gg) {
    const int n = static_cast<int>(pos.size());
    const int n3 = 3 * n;
    RealMatrix dip = xt::zeros<double>({n3, n3});
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j)
                continue;
            const double Rx = pos[static_cast<std::size_t>(i)][0]
                              - pos[static_cast<std::size_t>(j)][0];
            const double Ry = pos[static_cast<std::size_t>(i)][1]
                              - pos[static_cast<std::size_t>(j)][1];
            const double Rz = pos[static_cast<std::size_t>(i)][2]
                              - pos[static_cast<std::size_t>(j)][2];
            const double r =
                std::sqrt(Rx * Rx + Ry * Ry + Rz * Rz);
            if (r < 1e-15)
                continue;
            const double s_vdw =
                beta
                * (R_vdw[static_cast<std::size_t>(i)]
                   + R_vdw[static_cast<std::size_t>(j)]);
            const double f = damping_fermi(r, s_vdw, a_fermi);
            double T[3][3];
            if (gg) {
                const double sig = std::sqrt(
                    sigma[static_cast<std::size_t>(i)]
                        * sigma[static_cast<std::size_t>(i)]
                    + sigma[static_cast<std::size_t>(j)]
                          * sigma[static_cast<std::size_t>(j)]);
                T_erf_coulomb(Rx, Ry, Rz, sig, T);
                const double pref = 1.0 - f;
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        T[a][b] *= pref;
            } else {
                T_bare(Rx, Ry, Rz, T);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        T[a][b] *= f;
            }
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    // Column-major: row = 3*i+a, col = 3*j+b
                    dip(3 * i + a, 3 * j + b) += T[a][b];
                }
            }
        }
    }
    return dip;
}

/// @brief α_nlc = (α0^{-1} I + T)^{-1}  (3N×3N).
RealMatrix screening_matrix(const std::vector<double>& alpha_0,
                            const RealMatrix& dip) {
    const int n = static_cast<int>(alpha_0.size());
    const int n3 = 3 * n;
    RealMatrix A = dip;
    for (int i = 0; i < n; ++i) {
        const double inv_a = 1.0 / alpha_0[static_cast<std::size_t>(i)];
        for (int a = 0; a < 3; ++a)
            A(3 * i + a, 3 * i + a) += inv_a;
    }
    return xt::linalg::inv(A);
}

/// @brief Contract 3N×3N screened matrix → isotropic per-atom α.
std::vector<double> contract_atomic_alpha(const RealMatrix& a_nlc, int n) {
    std::vector<double> a(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        double tr = 0.0;
        for (int a = 0; a < 3; ++a)
            tr += a_nlc(3 * i + a, 3 * i + a);
        // sum over second atom of diagonal blocks: einsum aii->a after reshape
        // Full contraction: sum_j sum_a a_nlc[3i+a, 3j+a] / 3
        double sum = 0.0;
        for (int j = 0; j < n; ++j)
            for (int a = 0; a < 3; ++a)
                sum += a_nlc(3 * i + a, 3 * j + a);
        a[static_cast<std::size_t>(i)] = sum / 3.0;
        (void)tr;
    }
    return a;
}

void freq_grid(std::vector<double>& freq, std::vector<double>& freq_w) {
    constexpr double L = 0.6;
    freq.resize(kNFreq);
    freq_w.resize(kNFreq);
    freq[0] = 0.0;
    freq_w[0] = 0.0;
    // x,w from leggauss; map to half-line then reverse (libMBD).
    std::array<double, kNLeg> xmap{}, wmap{};
    for (int i = 0; i < kNLeg; ++i) {
        const double x = kLegX[i];
        const double w = kLegW[i];
        wmap[static_cast<std::size_t>(i)] = 2.0 * L / ((1.0 - x) * (1.0 - x)) * w;
        xmap[static_cast<std::size_t>(i)] = L * (1.0 + x) / (1.0 - x);
    }
    for (int i = 0; i < kNLeg; ++i) {
        // reverse
        freq[static_cast<std::size_t>(1 + i)] =
            xmap[static_cast<std::size_t>(kNLeg - 1 - i)];
        freq_w[static_cast<std::size_t>(1 + i)] =
            wmap[static_cast<std::size_t>(kNLeg - 1 - i)];
    }
}

} // namespace

MbdResult mbd_rsscs_pbe(const std::vector<int>& Z,
                        const std::vector<std::array<double, 3>>& pos_bohr,
                        const double cell_bohr[3], bool periodic, double beta,
                        double a_fermi, const std::vector<double>& vol_ratio) {
    (void)cell_bohr;
    const std::size_t n = Z.size();
    if (pos_bohr.size() != n)
        throw std::runtime_error("GPW MBD: Z and pos size mismatch");
    if (!vol_ratio.empty() && vol_ratio.size() != n)
        throw std::runtime_error("GPW MBD: vol_ratio size mismatch");
    if (periodic)
        throw std::runtime_error(
            "GPW MBD: periodic Ewald MBD is not in this product path — "
            "use isolated free-atom MBD@rsSCS (set nonperiodic or call "
            "mbd_rsscs_pbe with periodic=false)");
    if (!(beta > 0.0) || !(a_fermi > 0.0))
        throw std::runtime_error("GPW MBD: beta and a_fermi must be positive");

    MbdResult res;
    if (n == 0)
        return res;

    std::vector<double> alpha0(n), c6(n), rvdw(n);
    for (std::size_t i = 0; i < n; ++i) {
        free_atom_params(Z[i], alpha0[i], c6[i], rvdw[i]);
        const double v = vol_ratio.empty() ? 1.0 : vol_ratio[i];
        if (!(v > 0.0) || !std::isfinite(v))
            throw std::runtime_error(
                "GPW MBD: volume ratio must be positive finite");
        alpha0[i] *= v;
        c6[i] *= v * v;
        rvdw[i] *= std::cbrt(v);
    }

    // Characteristic frequencies ω = (4/3) C6 / α²
    std::vector<double> omega(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double a2 = alpha0[i] * alpha0[i];
        if (!(a2 > 0.0))
            throw std::runtime_error("GPW MBD: alpha_0 must be positive");
        omega[i] = (4.0 / 3.0) * c6[i] / a2;
    }

    std::vector<double> freq, freq_w;
    freq_grid(freq, freq_w);

    // Dynamical polarizabilities on the frequency grid → rsSCS screening.
    // alpha_dyn_rsscs[ifreq][iatom]
    const int nf = static_cast<int>(freq.size());
    std::vector<std::vector<double>> alpha_dyn_rsscs(
        static_cast<std::size_t>(nf), std::vector<double>(n, 0.0));

    for (int iu = 0; iu < nf; ++iu) {
        const double u = freq[static_cast<std::size_t>(iu)];
        std::vector<double> a_dyn(n);
        std::vector<double> sigma(n);
        for (std::size_t i = 0; i < n; ++i) {
            a_dyn[i] = alpha0[i] / (1.0 + (u / omega[i]) * (u / omega[i]));
            // σ = (√(2/π) α / 3)^{1/3}
            sigma[i] = std::cbrt(std::sqrt(2.0 / M_PI) * a_dyn[i] / 3.0);
        }
        const RealMatrix dip =
            dipole_matrix_isolated(pos_bohr, rvdw, sigma, beta, a_fermi,
                                   /*gg=*/true);
        const RealMatrix a_nlc = screening_matrix(a_dyn, dip);
        alpha_dyn_rsscs[static_cast<std::size_t>(iu)] =
            contract_atomic_alpha(a_nlc, static_cast<int>(n));
    }

    // C6_rsscs = 3/π Σ_w w(u) α(u)² ; R_vdw_rsscs = R * (α(0)/α0)^{1/3}
    std::vector<double> alpha_rsscs(n), c6_rsscs(n), rvdw_rsscs(n);
    for (std::size_t i = 0; i < n; ++i) {
        alpha_rsscs[i] = alpha_dyn_rsscs[0][i];
        double s = 0.0;
        for (int iu = 0; iu < nf; ++iu) {
            const double au = alpha_dyn_rsscs[static_cast<std::size_t>(iu)][i];
            s += freq_w[static_cast<std::size_t>(iu)] * au * au;
        }
        c6_rsscs[i] = 3.0 / M_PI * s;
        const double ratio =
            (alpha0[i] > 0.0) ? (alpha_rsscs[i] / alpha0[i]) : 1.0;
        if (!(ratio > 0.0))
            throw std::runtime_error("GPW MBD: screened alpha non-positive");
        rvdw_rsscs[i] = rvdw[i] * std::cbrt(ratio);
    }

    // MBD Hamiltonian eigenvalues
    std::vector<double> omega_rsscs(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double a2 = alpha_rsscs[i] * alpha_rsscs[i];
        if (!(a2 > 0.0))
            throw std::runtime_error("GPW MBD: screened alpha_0 must be positive");
        omega_rsscs[i] = (4.0 / 3.0) * c6_rsscs[i] / a2;
    }

    std::vector<double> sigma_dummy(n, 0.0); // unused for fermi,dip
    const RealMatrix Tlr = dipole_matrix_isolated(
        pos_bohr, rvdw_rsscs, sigma_dummy, beta, a_fermi, /*gg=*/false);

    const int n3 = 3 * static_cast<int>(n);
    RealMatrix C = xt::zeros<double>({n3, n3});
    // pre_i = ω_i √α_i  (repeated for x,y,z)
    std::vector<double> pre(static_cast<std::size_t>(n3));
    for (std::size_t i = 0; i < n; ++i) {
        const double p = omega_rsscs[i] * std::sqrt(alpha_rsscs[i]);
        for (int a = 0; a < 3; ++a)
            pre[3 * i + static_cast<std::size_t>(a)] = p;
    }
    for (int p = 0; p < n3; ++p) {
        C(p, p) +=
            omega_rsscs[static_cast<std::size_t>(p / 3)]
            * omega_rsscs[static_cast<std::size_t>(p / 3)];
        for (int q = 0; q < n3; ++q)
            C(p, q) += pre[static_cast<std::size_t>(p)]
                       * pre[static_cast<std::size_t>(q)] * Tlr(p, q);
    }

    auto eigs = xt::linalg::eigvalsh(C); // ascending
    double sum_sqrt = 0.0;
    for (int k = 0; k < n3; ++k) {
        const double lam = eigs(k);
        if (lam < -1e-10)
            throw std::runtime_error(
                "GPW MBD: negative eigenvalue of dynamical matrix ("
                + std::to_string(lam) + ")");
        sum_sqrt += std::sqrt(std::max(lam, 0.0));
    }
    double sum_omega = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        sum_omega += omega_rsscs[i];

    res.energy = 0.5 * sum_sqrt - 1.5 * sum_omega;
    return res;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
