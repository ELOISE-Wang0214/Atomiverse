/// @file longrange/pme_influence.h
/// @brief Backend-agnostic PME influence function and force interpolation
///        (no Frame / Component / MPI / FFT engine).
///
/// Smooth PME of Essmann et al., *J. Chem. Phys.* **103**, 8577 (1995).
/// DOI: 10.1063/1.470117. Operates on a complex grid view after the forward
/// FFT (influence) and after the inverse FFT (force interpolation). FFT
/// execution stays in the backend wrapper.
///
/// Units: energy / virial [eV], force [eV/Å], volume [Å³], α [1/Å].
#pragma once

#include <cmath>
#include <complex>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/precision.h"
#include "atomiverse/core/units.h"

namespace ATV {
namespace longrange {

/// @brief Apply the PME influence function in place: energy, virial, H·Q(m).
///
/// Reciprocal energy (Essmann et al.):
///   E = (k_e/(2πV)) Σ_{m≠0} exp(−π²|m|²/α²)/|m|² · B(m) · |Q(m)|²
/// with m built from inv_cell. Per-mode virial:
///   V_ab = ug · (δ_ab + vt · m_a m_b),  vt = −2(1/|m|² + π²/α²).
/// Grid layout: idx = m2 + K2·(m1 + K1·m0). Zero mode is zeroed.
///
/// @param[in,out] grid Complex grid (post forward FFT), length K0·K1·K2;
///                     on exit each mode is multiplied by H(m).
/// @param K            Grid dimensions {K0, K1, K2}.
/// @param bsp_mod0,1,2 Euler moduli |b_d(m)|², length K_d each.
/// @param inv_cell     Inverse cell matrix [1/Å].
/// @param alpha        Ewald α [1/Å].
/// @param volume       Cell volume [Å³].
/// @param[out] e_recip Scaled reciprocal energy [eV].
/// @param[out] stress  Scaled virial Voigt xx,yy,zz,xy,xz,yz [eV].
inline void pme_apply_influence(std::complex<double>* grid, const int K[3],
                                const double* bsp_mod0, const double* bsp_mod1,
                                const double* bsp_mod2, const Mat3& inv_cell,
                                double alpha, double volume, double& e_recip,
                                double stress[6]) {
    const int K0 = K[0], K1 = K[1], K2 = K[2];
    const double inv_pi2_a2 =
        ATV::units::pi * ATV::units::pi / (alpha * alpha);
    const double vt_const = inv_pi2_a2; // π²/α² term inside vt
    e_recip = 0;
    double s_xx = 0, s_yy = 0, s_zz = 0, s_xy = 0, s_xz = 0, s_yz = 0;
    for (int m0 = 0; m0 < K0; ++m0) {
        const int h0 = (m0 <= K0 / 2) ? m0 : m0 - K0;
        const double bm0 = bsp_mod0[m0];
        for (int m1 = 0; m1 < K1; ++m1) {
            const int h1 = (m1 <= K1 / 2) ? m1 : m1 - K1;
            const double bm01 = bm0 * bsp_mod1[m1];
            for (int m2 = 0; m2 < K2; ++m2) {
                const int idx = m2 + K2 * (m1 + K1 * m0);
                if (m0 == 0 && m1 == 0 && m2 == 0) {
                    grid[idx] = 0;
                    continue;
                }
                const int h2 = (m2 <= K2 / 2) ? m2 : m2 - K2;
                // Reciprocal vector m = Σ_d h_d · (column d of cell⁻¹).
                double mx = 0, my = 0, mz = 0;
                const int hh[3] = {h0, h1, h2};
                for (int d = 0; d < 3; ++d) {
                    mx += hh[d] * inv_cell(0, d);
                    my += hh[d] * inv_cell(1, d);
                    mz += hh[d] * inv_cell(2, d);
                }
                const double msq = mx * mx + my * my + mz * mz;
                const double B = bm01 * bsp_mod2[m2];
                const double H = std::exp(-inv_pi2_a2 * msq) / msq * B;
                const std::complex<double> Qm = grid[idx];
                const double q2 = Qm.real() * Qm.real() + Qm.imag() * Qm.imag();
                e_recip += H * q2;
                const double vt = -2.0 * (1.0 / msq + vt_const);
                const double ug = H * q2;
                s_xx += ug * (1.0 + vt * mx * mx);
                s_yy += ug * (1.0 + vt * my * my);
                s_zz += ug * (1.0 + vt * mz * mz);
                s_xy += ug * (vt * mx * my);
                s_xz += ug * (vt * mx * mz);
                s_yz += ug * (vt * my * mz);
                grid[idx] = H * Qm;
            }
        }
    }
    const double escale = ATV::units::coulomb / (2.0 * ATV::units::pi * volume);
    e_recip *= escale;
    stress[0] = s_xx * escale;
    stress[1] = s_yy * escale;
    stress[2] = s_zz * escale;
    stress[3] = s_xy * escale;
    stress[4] = s_xz * escale;
    stress[5] = s_yz * escale;
}

/// @brief Interpolate reciprocal forces onto owned atoms from the inverse-FFT
///        convolution grid.
///
/// F_{i,a} = −(units::coulomb/πV) q_i Σ_stencil C(grid) Σ_d dθ_d Π_{d'≠d} θ_{d'}
///           · K_d · inv_cell(a,d)
///
/// @param q       Charges [e], length @p n_atoms.
/// @param n_atoms Atom count.
/// @param order   B-spline order.
/// @param K       Grid dimensions.
/// @param inv_cell Inverse cell [1/Å].
/// @param volume  Cell volume [Å³].
/// @param base,theta,dtheta From #pme_fill_atom_splines (same layout).
/// @param grid    Real-space convolution grid (post inverse FFT); only .real()
///                is used.
/// @param[in,out] fx,fy,fz Force accumulators [eV/Å], length @p n_atoms.
inline void pme_interpolate_forces(const StorageReal* q, int n_atoms, int order,
                                   const int K[3], const Mat3& inv_cell,
                                   double volume, const int* base,
                                   const double* theta, const double* dtheta,
                                   const std::complex<double>* grid, AccReal* fx,
                                   AccReal* fy, AccReal* fz) {
    const int K0 = K[0], K1 = K[1], K2 = K[2];
    const int n = order;
    const int nm1 = n - 1;
    const double fscale = ATV::units::coulomb / (ATV::units::pi * volume);
    for (int i = 0; i < n_atoms; ++i) {
        const double qi = q[i];
        if (qi == 0.0)
            continue;
        const double* thx = &theta[(3 * i + 0) * n];
        const double* thy = &theta[(3 * i + 1) * n];
        const double* thz = &theta[(3 * i + 2) * n];
        const double* dthx = &dtheta[(3 * i + 0) * n];
        const double* dthy = &dtheta[(3 * i + 1) * n];
        const double* dthz = &dtheta[(3 * i + 2) * n];
        double g_d[3] = {0, 0, 0}; // gradient in (scaled fractional) dims
        for (int a = 0; a < n; ++a) {
            const int g0 = ((base[3 * i + 0] - nm1 + a) % K0 + K0) % K0;
            for (int b = 0; b < n; ++b) {
                const int g1 = ((base[3 * i + 1] - nm1 + b) % K1 + K1) % K1;
                for (int c = 0; c < n; ++c) {
                    const int g2 = ((base[3 * i + 2] - nm1 + c) % K2 + K2) % K2;
                    const int idx = g2 + K2 * (g1 + K1 * g0);
                    const double cval = grid[idx].real();
                    g_d[0] += cval * dthx[a] * thy[b] * thz[c];
                    g_d[1] += cval * thx[a] * dthy[b] * thz[c];
                    g_d[2] += cval * thx[a] * thy[b] * dthz[c];
                }
            }
        }
        // Chain rule: ∂E/∂r_a = fscale·q_i·Σ_d g_d·K_d·ic(a,d). Force = −∂E/∂r.
        double fxi = 0, fyi = 0, fzi = 0;
        for (int d = 0; d < 3; ++d) {
            const double gk = g_d[d] * K[d];
            fxi += gk * inv_cell(0, d);
            fyi += gk * inv_cell(1, d);
            fzi += gk * inv_cell(2, d);
        }
        fx[i] += -fscale * qi * fxi;
        fy[i] += -fscale * qi * fyi;
        fz[i] += -fscale * qi * fzi;
    }
}

} // namespace longrange
} // namespace ATV
