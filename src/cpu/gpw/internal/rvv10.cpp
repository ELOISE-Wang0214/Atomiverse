/// @file cpu/gpw/internal/rvv10.cpp
/// @brief RPS evaluation of rVV10 nonlocal correlation (CP2K TYPE RVV10).

#include "atomiverse/cpu/gpw/internal/rvv10.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include <fftw3.h>

#include "atomiverse/cpu/gpw/internal/rvv10_kernel_data.h"

namespace ATV {
namespace cpu {
namespace gpw {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kRhoCut = 1e-10;

/// @brief Natural cubic spline second derivatives of Kronecker-δ basis on q-mesh
///        (CP2K initialize_spline_interpolation).
void init_p_spline_d2(const double* x, int n, std::vector<double>& d2) {
    // d2 layout [P_i][index], n × n
    d2.assign(static_cast<std::size_t>(n * n), 0.0);
    std::vector<double> y(static_cast<std::size_t>(n));
    std::vector<double> temp(static_cast<std::size_t>(n));
    for (int P_i = 0; P_i < n; ++P_i) {
        for (int j = 0; j < n; ++j)
            y[static_cast<std::size_t>(j)] = 0.0;
        y[static_cast<std::size_t>(P_i)] = 1.0;
        d2[static_cast<std::size_t>(P_i * n + 0)] = 0.0;
        temp[0] = 0.0;
        for (int index = 1; index < n - 1; ++index) {
            const double temp1 =
                (x[index] - x[index - 1]) / (x[index + 1] - x[index - 1]);
            const double temp2 =
                temp1 * d2[static_cast<std::size_t>(P_i * n + index - 1)]
                + 2.0;
            d2[static_cast<std::size_t>(P_i * n + index)] =
                (temp1 - 1.0) / temp2;
            temp[static_cast<std::size_t>(index)] =
                (y[static_cast<std::size_t>(index + 1)]
                 - y[static_cast<std::size_t>(index)])
                    / (x[index + 1] - x[index])
                - (y[static_cast<std::size_t>(index)]
                   - y[static_cast<std::size_t>(index - 1)])
                      / (x[index] - x[index - 1]);
            temp[static_cast<std::size_t>(index)] =
                (6.0 * temp[static_cast<std::size_t>(index)]
                     / (x[index + 1] - x[index - 1])
                 - temp1 * temp[static_cast<std::size_t>(index - 1)])
                / temp2;
        }
        d2[static_cast<std::size_t>(P_i * n + (n - 1))] = 0.0;
        for (int index = n - 2; index >= 0; --index) {
            d2[static_cast<std::size_t>(P_i * n + index)] =
                d2[static_cast<std::size_t>(P_i * n + index)]
                    * d2[static_cast<std::size_t>(P_i * n + index + 1)]
                + temp[static_cast<std::size_t>(index)];
        }
    }
}

/// @brief Evaluate all p_α(q) via cubic spline on Kronecker basis (CP2K).
void palpha_at_q(const double* x, const std::vector<double>& d2, int n,
                 double q, double* out) {
    int lower = 0;
    int upper = n - 1;
    if (q <= x[0]) {
        for (int a = 0; a < n; ++a)
            out[a] = (a == 0) ? 1.0 : 0.0;
        return;
    }
    if (q >= x[n - 1]) {
        for (int a = 0; a < n; ++a)
            out[a] = (a == n - 1) ? 1.0 : 0.0;
        return;
    }
    while (upper - lower > 1) {
        const int mid = (upper + lower) / 2;
        if (q > x[mid])
            lower = mid;
        else
            upper = mid;
    }
    const double dx = x[upper] - x[lower];
    const double const_ = dx * dx / 6.0;
    const double a = (x[upper] - q) / dx;
    const double b = (q - x[lower]) / dx;
    const double c = (a * a * a - a) * const_;
    const double d = (b * b * b - b) * const_;
    for (int P_i = 0; P_i < n; ++P_i) {
        // y is Kronecker: y[lower]=δ, y[upper]=δ
        const double yl = (P_i == lower) ? 1.0 : 0.0;
        const double yu = (P_i == upper) ? 1.0 : 0.0;
        out[P_i] = a * yl + b * yu
                   + (c * d2[static_cast<std::size_t>(P_i * n + lower)]
                      + d * d2[static_cast<std::size_t>(P_i * n + upper)]);
    }
}

void interpolate_kernel_linear(double k, double* kernel_aa) {
    using namespace rvv10_data;
    const float* tab = kernel_kaa();
    const int asqr = kNAlpha * kNAlpha;
    const double index = k / kDeltaK;
    int istart = static_cast<int>(std::floor(index));
    if (k < 0.0 || istart > kNKpoints - 2) {
        for (int i = 0; i < asqr; ++i)
            kernel_aa[i] = 0.0;
        return;
    }
    if (istart < 0)
        istart = 0;
    const double t = index - static_cast<double>(istart);
    const double t2 = 1.0 - t;
    const float* left = tab + static_cast<std::size_t>(istart) * asqr;
    const float* right = tab + static_cast<std::size_t>(istart + 1) * asqr;
    for (int i = 0; i < asqr; ++i)
        kernel_aa[i] = static_cast<double>(left[i]) * t2
                       + static_cast<double>(right[i]) * t;
}

} // namespace

void rvv10_q0(const std::vector<double>& rho, const std::vector<double>& sigma,
              std::vector<double>& q0, double b, double C) {
    const std::size_t n = rho.size();
    if (sigma.size() != n)
        throw std::runtime_error("rVV10 q0: rho/sigma size mismatch");
    q0.resize(n);
    const double qcut = kRvv10QCut;
    const double qmin = kRvv10QMin;
    for (std::size_t i = 0; i < n; ++i) {
        const double mrho = rho[i];
        if (mrho < kRhoCut) {
            q0[i] = qcut;
            continue;
        }
        const double gmod2 = sigma[i]; // |∇ρ|²
        // wp² = 16 π ρ ; wg² = 4 C (|∇ρ|/ρ)⁴ = 4 C (|∇ρ|²)² / ρ⁴
        const double wp2 = 16.0 * kPi * mrho;
        const double rho2 = mrho * mrho;
        const double rho4 = rho2 * rho2;
        const double wg2 = 4.0 * C * (gmod2 * gmod2) / rho4;
        // κ = b · 3π · (ρ/(9π))^{1/6}
        const double k =
            b * 3.0 * kPi * std::pow(mrho / (9.0 * kPi), 1.0 / 6.0);
        const double w0 = std::sqrt(wg2 + wp2 / 3.0);
        double q = w0 / k;
        // Soft saturation (Soler eq. 5), m_cut = 12
        double exponent = 0.0;
        double qscaled0 = q / qcut;
        double qscaled = 1.0;
        for (int m = 1; m <= 12; ++m) {
            qscaled *= qscaled0;
            exponent += qscaled / static_cast<double>(m);
        }
        double qsat = qcut * (1.0 - std::exp(-exponent));
        if (qsat < qmin)
            qsat = qmin;
        q0[i] = qsat;
    }
}

Rvv10Result rvv10_energy(const UniformGrid& grid, const std::vector<double>& rho,
                         const std::vector<double>& sigma, double b, double C) {
    using namespace rvv10_data;
    const int n0 = grid.n[0];
    const int n1 = grid.n[1];
    const int n2 = grid.n[2];
    const long ntot = static_cast<long>(n0) * n1 * n2;
    if (static_cast<long>(rho.size()) != ntot
        || static_cast<long>(sigma.size()) != ntot)
        throw std::runtime_error("rVV10: density size mismatch vs grid");

    std::vector<double> q0;
    rvv10_q0(rho, sigma, q0, b, C);

    // θ_α = const * ρ^{3/4} * p_α(q0)
    // const = 1 / (3 √π b^{3/2} π^{3/4})   [CP2K qs_dispersion_nonloc.F]
    const double rootpi = std::sqrt(kPi);
    const double const_theta =
        1.0 / (3.0 * rootpi * std::pow(b, 1.5) * std::pow(kPi, 0.75));

    // β = (1/32) (3 / b²)^{3/4}
    const double beta = 0.03125 * std::pow(3.0 / (b * b), 0.75);

    const double* qm = qmesh();
    std::vector<double> d2y;
    init_p_spline_d2(qm, kNAlpha, d2y);

    const int n2c = n2 / 2 + 1;
    const long ncplx = static_cast<long>(n0) * n1 * n2c;
    const int nalpha = kNAlpha;

    double* rbuf = static_cast<double*>(
        fftw_malloc(sizeof(double) * static_cast<std::size_t>(ntot * nalpha)));
    fftw_complex* cbuf = static_cast<fftw_complex*>(fftw_malloc(
        sizeof(fftw_complex) * static_cast<std::size_t>(ncplx * nalpha)));
    if (!rbuf || !cbuf)
        throw std::runtime_error("rVV10: FFTW alloc failed");

    std::vector<double> ptmp(static_cast<std::size_t>(nalpha));
    double nelec = 0.0;
    const double dV = grid.weight();
    for (long g = 0; g < ntot; ++g) {
        const double rg = rho[static_cast<std::size_t>(g)];
        nelec += rg;
        if (rg < kRhoCut) {
            for (int a = 0; a < nalpha; ++a)
                rbuf[a * ntot + g] = 0.0;
            continue;
        }
        palpha_at_q(qm, d2y, nalpha, q0[static_cast<std::size_t>(g)],
                    ptmp.data());
        const double pref = const_theta * std::pow(rg, 0.75);
        for (int a = 0; a < nalpha; ++a)
            rbuf[a * ntot + g] = pref * ptmp[static_cast<std::size_t>(a)];
    }
    nelec *= dV;

    for (int a = 0; a < nalpha; ++a) {
        double* rin = rbuf + a * ntot;
        fftw_complex* cout = cbuf + a * ncplx;
        fftw_plan fwd =
            fftw_plan_dft_r2c_3d(n0, n1, n2, rin, cout, FFTW_ESTIMATE);
        fftw_execute(fwd);
        fftw_destroy_plan(fwd);
    }

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
                const double gz = i2 * gz_unit;
                const long kindex =
                    (static_cast<long>(i0) * n1 + i1) * n2c + i2;
                const int weight =
                    (i2 == 0 || (n2 % 2 == 0 && i2 == n2c - 1)) ? 1 : 2;
                const double k = std::sqrt(gx * gx + gy * gy + gz * gz);
                interpolate_kernel_linear(k, kernel_aa.data());

                for (int a1 = 0; a1 < nalpha; ++a1) {
                    std::complex<double> F(0.0, 0.0);
                    for (int a2 = 0; a2 < nalpha; ++a2) {
                        const fftw_complex& th =
                            cbuf[a2 * ncplx + kindex];
                        const double phi = kernel_aa[static_cast<std::size_t>(
                            a1 * nalpha + a2)];
                        F += std::complex<double>(th[0], th[1]) * phi;
                    }
                    const fftw_complex& th1 = cbuf[a1 * ncplx + kindex];
                    energy += weight
                              * (th1[0] * F.real() + th1[1] * F.imag());
                }
            }
        }
    }

    // Same FFT convention as DRSLL (matches unnormalized r2c + dV/N):
    //   E_kernel_raw = ½ * (dV/N) * sum_k w Re(θ* φ θ)
    // CP2K applies an extra ½ for RVV10 after its own ½·vol factor.
    energy *= 0.5 * dV / static_cast<double>(ntot);
    const double e_kernel = 0.5 * energy;
    const double e_beta = beta * nelec;

    fftw_free(rbuf);
    fftw_free(cbuf);

    Rvv10Result res;
    res.e_kernel = e_kernel;
    res.e_beta = e_beta;
    res.energy = e_kernel + e_beta;
    return res;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
