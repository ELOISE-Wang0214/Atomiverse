/// @file cpu/gpw/internal/vdw_df.cpp
/// @brief RPS evaluation of vdW-DF (DRSLL / DF1) and vdW-DF2 (LMKLL / DF2).

#include "atomiverse/cpu/gpw/internal/vdw_df.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include <fftw3.h>

#include "atomiverse/cpu/gpw/internal/vdw_df_kernel_data.h"

namespace ATV {
namespace cpu {
namespace gpw {
namespace {

// Perdew–Wang LDA constants used inside q₀ (libvdwxc / Dion).
constexpr double kPi = 3.14159265358979323846;
constexpr double kKfPref = 3.0936677262801355; // (3 π²)^{1/3}
constexpr double kRsPref = 0.6203504908994001; // (3/(4π))^{1/3}
constexpr double kLdaA = 0.031091;
constexpr double kLdaA1 = 0.2137;
constexpr double kLdaB1 = 7.5957;
constexpr double kLdaB2 = 3.5876;
constexpr double kLdaB3 = 1.6382;
constexpr double kLdaB4 = 0.49294;
constexpr double kRhoCut = 1e-12;

void interpolate_kernel(double k, double* kernel_aa) {
    using namespace vdw_df_data;
    const float* tab = kernel_kaa();
    const int asqr = kNAlpha * kNAlpha;
    const double index = k / kDeltaK;
    int istart = static_cast<int>(std::floor(index));
    if (k < 0.0 || istart > kNKpoints - 2) {
        for (int i = 0; i < asqr; ++i)
            kernel_aa[i] = 0.0;
        return;
    }
    const double t = index - static_cast<double>(istart);
    const double t2 = 1.0 - t;
    const float* left = tab + static_cast<std::size_t>(istart) * asqr;
    const float* right = tab + static_cast<std::size_t>(istart + 1) * asqr;
    for (int i = 0; i < asqr; ++i)
        kernel_aa[i] = static_cast<double>(left[i]) * t2
                       + static_cast<double>(right[i]) * t;
}

/// @brief Evaluate p_α(q) cubic splines; write nalpha values into out.
void palpha_splines(double q, double prefactor, double* out) {
    using namespace vdw_df_data;
    const double* qm = qmesh();
    const float* spl = splines_arc();
    int alpha = 0;
    if (q > qm[kNAlpha - 1]) {
        // Soft-cutoff already keeps q ≤ q_cut; clamp.
        q = qm[kNAlpha - 1];
        alpha = kNAlpha - 2;
    } else if (q <= qm[0]) {
        q = qm[0];
        alpha = 0;
    } else {
        for (alpha = 1; alpha < kNAlpha; ++alpha) {
            if (qm[alpha] >= q) {
                --alpha;
                break;
            }
        }
    }
    const double dq = q - qm[alpha];
    for (int beta = 0; beta < kNAlpha; ++beta) {
        const float* poly =
            spl + (beta * kNRegions + alpha) * kNPCoef;
        out[beta] = prefactor
                    * (static_cast<double>(poly[0])
                       + dq
                             * (static_cast<double>(poly[1])
                                + dq
                                      * (static_cast<double>(poly[2])
                                         + dq * static_cast<double>(poly[3]))));
    }
}

/// @brief Shared RPS energy with explicit Z_ab (DF1 or DF2).
VdwDfResult vdw_df_rps(const UniformGrid& grid, const std::vector<double>& rho,
                       const std::vector<double>& sigma, double zab) {
    using namespace vdw_df_data;
    const int n0 = grid.n[0];
    const int n1 = grid.n[1];
    const int n2 = grid.n[2];
    const long ntot = static_cast<long>(n0) * n1 * n2;
    if (static_cast<long>(rho.size()) != ntot
        || static_cast<long>(sigma.size()) != ntot)
        throw std::runtime_error("vdW-DF: density size mismatch vs grid");

    std::vector<double> q0;
    vdw_df_q0(rho, sigma, q0, zab);

    const int n2c = n2 / 2 + 1;
    const long ncplx = static_cast<long>(n0) * n1 * n2c;
    const int nalpha = kNAlpha;

    // θ_α(r) then FFT → θ_α(k). Layout real: [g][α]; complex: [k][α].
    double* rbuf = static_cast<double*>(
        fftw_malloc(sizeof(double) * static_cast<std::size_t>(ntot * nalpha)));
    fftw_complex* cbuf = static_cast<fftw_complex*>(fftw_malloc(
        sizeof(fftw_complex) * static_cast<std::size_t>(ncplx * nalpha)));
    if (!rbuf || !cbuf)
        throw std::runtime_error("vdW-DF: FFTW alloc failed");

    // Build θ_α = ρ p_α(q0) as α-major real planes.
    std::vector<double> ptmp(static_cast<std::size_t>(nalpha));
    for (long g = 0; g < ntot; ++g) {
        palpha_splines(q0[static_cast<std::size_t>(g)],
                       rho[static_cast<std::size_t>(g)], ptmp.data());
        for (int a = 0; a < nalpha; ++a)
            rbuf[a * ntot + g] = ptmp[static_cast<std::size_t>(a)];
    }

    // FFT each α plane (r2c, same convention as poisson_fft).
    for (int a = 0; a < nalpha; ++a) {
        double* rin = rbuf + a * ntot;
        fftw_complex* cout = cbuf + a * ncplx;
        fftw_plan fwd =
            fftw_plan_dft_r2c_3d(n0, n1, n2, rin, cout, FFTW_ESTIMATE);
        fftw_execute(fwd);
        fftw_destroy_plan(fwd);
    }

    // Reciprocal lattice (orthorhombic): G_d = 2π m_d / L_d.
    const double two_pi = 2.0 * kPi;
    const double gx_unit = two_pi / grid.L[0];
    const double gy_unit = two_pi / grid.L[1];
    const double gz_unit = two_pi / grid.L[2];

    std::vector<double> kernel_aa(static_cast<std::size_t>(nalpha * nalpha));
    double energy = 0.0;

    for (int i0 = 0; i0 < n0; ++i0) {
        const int m0 = (i0 <= n0 / 2) ? i0 : i0 - n0;
        const double gx = m0 * gx_unit;
        for (int i1 = 0; i1 < n1; ++i1) {
            const int m1 = (i1 <= n1 / 2) ? i1 : i1 - n1;
            const double gy = m1 * gy_unit;
            for (int i2 = 0; i2 < n2c; ++i2) {
                const int m2 = i2;
                const double gz = m2 * gz_unit;
                const long kindex =
                    (static_cast<long>(i0) * n1 + i1) * n2c + i2;
                // Hermitian weight: endpoints of half-range once, interior ×2.
                const int weight =
                    (i2 == 0 || (n2 % 2 == 0 && i2 == n2c - 1)) ? 1 : 2;
                const double k = std::sqrt(gx * gx + gy * gy + gz * gz);
                interpolate_kernel(k, kernel_aa.data());

                // F_α = Σ_β φ_αβ θ_β ; energy += w Re[θ_α* F_α]
                for (int a1 = 0; a1 < nalpha; ++a1) {
                    std::complex<double> F(0.0, 0.0);
                    for (int a2 = 0; a2 < nalpha; ++a2) {
                        const fftw_complex& th =
                            cbuf[a2 * ncplx + kindex];
                        const double phi =
                            kernel_aa[static_cast<std::size_t>(a1 * nalpha
                                                               + a2)];
                        F += std::complex<double>(th[0], th[1]) * phi;
                    }
                    const fftw_complex& th1 = cbuf[a1 * ncplx + kindex];
                    energy += weight
                              * (th1[0] * F.real() + th1[1] * F.imag());
                }
            }
        }
    }

    const double dV = grid.weight(); // V / N
    energy *= 0.5 * dV / static_cast<double>(ntot);

    fftw_free(rbuf);
    fftw_free(cbuf);

    VdwDfResult res;
    res.energy = energy;
    return res;
}

} // namespace

void vdw_df_q0(const std::vector<double>& rho, const std::vector<double>& sigma,
               std::vector<double>& q0) {
    vdw_df_q0(rho, sigma, q0, kVdwDf1Zab);
}

void vdw_df_q0(const std::vector<double>& rho, const std::vector<double>& sigma,
               std::vector<double>& q0, double zab) {
    const std::size_t n = rho.size();
    if (sigma.size() != n)
        throw std::runtime_error("vdW-DF q0: rho/sigma size mismatch");
    q0.resize(n);
    const double qcut = kVdwDfQCut;
    for (std::size_t i = 0; i < n; ++i) {
        double mrho = rho[i];
        double msigma = sigma[i];
        if (mrho < kRhoCut) {
            q0[i] = qcut;
            continue;
        }
        const double rho1_root3 = std::cbrt(mrho);
        const double kF = kKfPref * rho1_root3;
        const double rs = kRsPref / rho1_root3;
        const double rs_sqr = rs * rs;
        const double rs_root2 = std::sqrt(rs);
        const double rho_sqr = mrho * mrho;
        const double gradient_correction =
            -zab / (36.0 * kF * rho_sqr) * msigma;
        const double LDA_1 =
            8.0 * kPi / 3.0 * (kLdaA * (1.0 + kLdaA1 * rs));
        const double LDA_2 =
            2.0 * kLdaA
            * (kLdaB1 * rs_root2 + kLdaB2 * rs + kLdaB3 * rs * rs_root2
               + kLdaB4 * rs_sqr);
        const double q =
            kF + LDA_1 * std::log(1.0 + 1.0 / LDA_2) + gradient_correction;
        double exponent = 0.0;
        double qscaled0 = q / qcut;
        double qscaled = 1.0;
        for (int m = 1; m <= 12; ++m) {
            qscaled *= qscaled0;
            exponent += qscaled / static_cast<double>(m);
        }
        q0[i] = qcut * (1.0 - std::exp(-exponent));
    }
}

VdwDfResult vdw_df_drsll(const UniformGrid& grid, const std::vector<double>& rho,
                         const std::vector<double>& sigma) {
    return vdw_df_rps(grid, rho, sigma, kVdwDf1Zab);
}

VdwDfResult vdw_df_lmkll(const UniformGrid& grid, const std::vector<double>& rho,
                         const std::vector<double>& sigma) {
    return vdw_df_rps(grid, rho, sigma, kVdwDf2Zab);
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
