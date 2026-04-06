/// @file cpu/scf/basis_eval.h
/// @brief Inline basis function evaluation on grid points (s+p type Gaussians).
///
/// Shared by scf_solve.cpp (V_xc construction) and gradient.cpp (XC gradient).
/// All coordinates in Bohr.
///
/// For s-type (l=0): phi = sum c_k exp(-alpha_k r^2)
/// For p-type (l=1): phi = d_c * sum c_k exp(-alpha_k r^2)
///   where d_c is the Cartesian displacement along sh.cart (0=x,1=y,2=z).
#pragma once

#include "atomiverse/cpu/scf/native_integrals.h"
#include <cmath>

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Evaluate contracted basis function at a grid point.
/// @param sh Shell descriptor (exponents, normalized coefficients, center, l, cart).
/// @param rx Grid point x [Bohr].
/// @param ry Grid point y [Bohr].
/// @param rz Grid point z [Bohr].
/// @return Basis function value phi(r).
inline double eval_basis(const native::Shell& sh, double rx, double ry, double rz) {
    double dx = rx - sh.center[0];
    double dy = ry - sh.center[1];
    double dz = rz - sh.center[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    double min_alpha = 1e30;
    for (int k = 0; k < sh.n_prim; ++k)
        if (sh.alpha[k] < min_alpha) min_alpha = sh.alpha[k];
    if (min_alpha * r2 > 36.0) return 0.0;
    double radial = 0;
    for (int k = 0; k < sh.n_prim; ++k)
        radial += sh.coeff[k] * std::exp(-sh.alpha[k] * r2);
    if (sh.l == 0) return radial;
    // p-type: phi = d_c * radial
    const double d[3] = {dx, dy, dz};
    return d[sh.cart] * radial;
}

/// @brief Evaluate spatial gradient of contracted basis function.
///
/// s-type: dphi/dr_i = sum_k c_k * (-2 alpha_k) * d_i * exp(...)
/// p-type: dphi/dr_i = delta(i,c)*f + d_c * df/dr_i
///   where f = radial part, df/dr_i = s-type gradient
///
/// @param sh Shell descriptor.
/// @param rx Grid point x [Bohr].
/// @param ry Grid point y [Bohr].
/// @param rz Grid point z [Bohr].
/// @param[out] gx dphi/dx.
/// @param[out] gy dphi/dy.
/// @param[out] gz dphi/dz.
inline void eval_basis_grad(const native::Shell& sh, double rx, double ry, double rz,
                            double& gx, double& gy, double& gz) {
    double dx = rx - sh.center[0];
    double dy = ry - sh.center[1];
    double dz = rz - sh.center[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    gx = gy = gz = 0;
    double min_alpha = 1e30;
    for (int k = 0; k < sh.n_prim; ++k)
        if (sh.alpha[k] < min_alpha) min_alpha = sh.alpha[k];
    if (min_alpha * r2 > 36.0) return;
    // Compute radial part f and its gradient df/dr_i (= s-type gradient)
    double f = 0, dfx = 0, dfy = 0, dfz = 0;
    for (int k = 0; k < sh.n_prim; ++k) {
        double e = sh.coeff[k] * std::exp(-sh.alpha[k] * r2);
        f += e;
        double fac = -2.0 * sh.alpha[k] * e;
        dfx += fac * dx; dfy += fac * dy; dfz += fac * dz;
    }
    if (sh.l == 0) {
        gx = dfx; gy = dfy; gz = dfz;
    } else {
        // p-type: d/dr_i [d_c * f] = delta(i,c)*f + d_c*df/dr_i
        const double d[3] = {dx, dy, dz};
        const double df[3] = {dfx, dfy, dfz};
        gx = d[sh.cart] * dfx;
        gy = d[sh.cart] * dfy;
        gz = d[sh.cart] * dfz;
        // Add delta(i,c)*f contribution
        double* g_arr[3] = {&gx, &gy, &gz};
        *g_arr[sh.cart] += f;
    }
}

/// @brief Evaluate spatial Hessian of contracted basis function.
///
/// s-type Hessian: d^2f/(dr_i dr_j) = sum_k c_k exp(...)
///       * [-2 alpha_k delta_{ij} + 4 alpha_k^2 d_i d_j]
/// p-type Hessian: d^2phi/(dr_i dr_j) = delta(i,c)*df/dr_j + delta(j,c)*df/dr_i
///       + d_c * d^2f/(dr_i dr_j)
///
/// @param sh  Shell descriptor.
/// @param rx  Grid point x [Bohr].
/// @param ry  Grid point y [Bohr].
/// @param rz  Grid point z [Bohr].
/// @param[out] hxx d^2 phi / dx^2.
/// @param[out] hxy d^2 phi / (dx dy).
/// @param[out] hxz d^2 phi / (dx dz).
/// @param[out] hyy d^2 phi / dy^2.
/// @param[out] hyz d^2 phi / (dy dz).
/// @param[out] hzz d^2 phi / dz^2.
inline void eval_basis_hessian(const native::Shell& sh, double rx, double ry, double rz,
                               double& hxx, double& hxy, double& hxz,
                               double& hyy, double& hyz, double& hzz) {
    double dx = rx - sh.center[0];
    double dy = ry - sh.center[1];
    double dz = rz - sh.center[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    hxx = hxy = hxz = hyy = hyz = hzz = 0;
    double min_alpha = 1e30;
    for (int k = 0; k < sh.n_prim; ++k)
        if (sh.alpha[k] < min_alpha) min_alpha = sh.alpha[k];
    if (min_alpha * r2 > 36.0) return;
    // Compute s-type radial quantities: f, df[3], d2f[6]
    double f = 0;
    double df[3] = {0, 0, 0};
    double d2f_xx = 0, d2f_xy = 0, d2f_xz = 0;
    double d2f_yy = 0, d2f_yz = 0, d2f_zz = 0;
    for (int k = 0; k < sh.n_prim; ++k) {
        double a = sh.alpha[k];
        double e = sh.coeff[k] * std::exp(-a * r2);
        double a2 = a * a;
        f += e;
        double fac = -2.0 * a * e;
        df[0] += fac * dx; df[1] += fac * dy; df[2] += fac * dz;
        d2f_xx += e * (-2.0*a + 4.0*a2*dx*dx);
        d2f_yy += e * (-2.0*a + 4.0*a2*dy*dy);
        d2f_zz += e * (-2.0*a + 4.0*a2*dz*dz);
        d2f_xy += e * 4.0*a2*dx*dy;
        d2f_xz += e * 4.0*a2*dx*dz;
        d2f_yz += e * 4.0*a2*dy*dz;
    }
    if (sh.l == 0) {
        hxx = d2f_xx; hxy = d2f_xy; hxz = d2f_xz;
        hyy = d2f_yy; hyz = d2f_yz; hzz = d2f_zz;
    } else {
        // p-type: d^2/(dr_i dr_j)[d_c*f] = delta(i,c)*df_j + delta(j,c)*df_i + d_c*d2f_ij
        const double d[3] = {dx, dy, dz};
        const int c = sh.cart;
        double d2f[3][3] = {
            {d2f_xx, d2f_xy, d2f_xz},
            {d2f_xy, d2f_yy, d2f_yz},
            {d2f_xz, d2f_yz, d2f_zz}
        };
        // h[i][j] = delta(i,c)*df[j] + delta(j,c)*df[i] + d[c]*d2f[i][j]
        double h[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                h[i][j] = d[c] * d2f[i][j]
                         + (i == c ? df[j] : 0.0)
                         + (j == c ? df[i] : 0.0);
        hxx = h[0][0]; hxy = h[0][1]; hxz = h[0][2];
        hyy = h[1][1]; hyz = h[1][2]; hzz = h[2][2];
    }
}

/// @brief Fused evaluation of basis function value, gradient, and Hessian.
///
/// Computes all three in a single pass over primitives, sharing the exp() call.
/// Handles both s-type (l=0) and p-type (l=1) basis functions.
///
/// @param sh  Shell descriptor.
/// @param rx  Grid point x [Bohr].
/// @param ry  Grid point y [Bohr].
/// @param rz  Grid point z [Bohr].
/// @param[out] val Basis function value phi(r).
/// @param[out] gx dphi/dx.
/// @param[out] gy dphi/dy.
/// @param[out] gz dphi/dz.
/// @param[out] hxx d^2 phi / dx^2.
/// @param[out] hxy d^2 phi / (dx dy).
/// @param[out] hxz d^2 phi / (dx dz).
/// @param[out] hyy d^2 phi / dy^2.
/// @param[out] hyz d^2 phi / (dy dz).
/// @param[out] hzz d^2 phi / dz^2.
inline void eval_basis_all(const native::Shell& sh, double rx, double ry, double rz,
                           double& val,
                           double& gx, double& gy, double& gz,
                           double& hxx, double& hxy, double& hxz,
                           double& hyy, double& hyz, double& hzz) {
    double dx = rx - sh.center[0];
    double dy = ry - sh.center[1];
    double dz = rz - sh.center[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    val = gx = gy = gz = 0;
    hxx = hxy = hxz = hyy = hyz = hzz = 0;
    double min_alpha = 1e30;
    for (int k = 0; k < sh.n_prim; ++k)
        if (sh.alpha[k] < min_alpha) min_alpha = sh.alpha[k];
    if (min_alpha * r2 > 36.0) return;
    // Compute s-type radial quantities
    double f = 0;
    double df_arr[3] = {0, 0, 0};
    double d2f_arr[6] = {0, 0, 0, 0, 0, 0}; // xx,xy,xz,yy,yz,zz
    for (int k = 0; k < sh.n_prim; ++k) {
        double a = sh.alpha[k];
        double e = sh.coeff[k] * std::exp(-a * r2);
        double a2 = a * a;
        f += e;
        double fac = -2.0 * a * e;
        df_arr[0] += fac * dx; df_arr[1] += fac * dy; df_arr[2] += fac * dz;
        d2f_arr[0] += e * (-2.0*a + 4.0*a2*dx*dx); // xx
        d2f_arr[1] += e * 4.0*a2*dx*dy;             // xy
        d2f_arr[2] += e * 4.0*a2*dx*dz;             // xz
        d2f_arr[3] += e * (-2.0*a + 4.0*a2*dy*dy); // yy
        d2f_arr[4] += e * 4.0*a2*dy*dz;             // yz
        d2f_arr[5] += e * (-2.0*a + 4.0*a2*dz*dz); // zz
    }
    if (sh.l == 0) {
        val = f;
        gx = df_arr[0]; gy = df_arr[1]; gz = df_arr[2];
        hxx = d2f_arr[0]; hxy = d2f_arr[1]; hxz = d2f_arr[2];
        hyy = d2f_arr[3]; hyz = d2f_arr[4]; hzz = d2f_arr[5];
    } else {
        const double d[3] = {dx, dy, dz};
        const int c = sh.cart;
        val = d[c] * f;
        // Gradient: delta(i,c)*f + d_c*df_i
        double g[3];
        for (int i = 0; i < 3; ++i)
            g[i] = d[c] * df_arr[i] + (i == c ? f : 0.0);
        gx = g[0]; gy = g[1]; gz = g[2];
        // Hessian: delta(i,c)*df_j + delta(j,c)*df_i + d_c*d2f_ij
        // d2f indexed as: [0]=xx [1]=xy [2]=xz [3]=yy [4]=yz [5]=zz
        double d2f[3][3] = {
            {d2f_arr[0], d2f_arr[1], d2f_arr[2]},
            {d2f_arr[1], d2f_arr[3], d2f_arr[4]},
            {d2f_arr[2], d2f_arr[4], d2f_arr[5]}
        };
        double h[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                h[i][j] = d[c] * d2f[i][j]
                         + (i == c ? df_arr[j] : 0.0)
                         + (j == c ? df_arr[i] : 0.0);
        hxx = h[0][0]; hxy = h[0][1]; hxz = h[0][2];
        hyy = h[1][1]; hyz = h[1][2]; hzz = h[2][2];
    }
}

/// @brief Fused evaluation of basis function value and gradient.
///
/// Single pass over primitives, sharing the exp() call.
/// Handles both s-type (l=0) and p-type (l=1) basis functions.
///
/// @param sh  Shell descriptor.
/// @param rx  Grid point x [Bohr].
/// @param ry  Grid point y [Bohr].
/// @param rz  Grid point z [Bohr].
/// @param[out] val Basis function value phi(r).
/// @param[out] gx dphi/dx.
/// @param[out] gy dphi/dy.
/// @param[out] gz dphi/dz.
inline void eval_basis_val_grad(const native::Shell& sh, double rx, double ry, double rz,
                                double& val, double& gx, double& gy, double& gz) {
    double dx = rx - sh.center[0];
    double dy = ry - sh.center[1];
    double dz = rz - sh.center[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    val = gx = gy = gz = 0;
    double min_alpha = 1e30;
    for (int k = 0; k < sh.n_prim; ++k)
        if (sh.alpha[k] < min_alpha) min_alpha = sh.alpha[k];
    if (min_alpha * r2 > 36.0) return;
    double f = 0, dfx = 0, dfy = 0, dfz = 0;
    for (int k = 0; k < sh.n_prim; ++k) {
        double a = sh.alpha[k];
        double e = sh.coeff[k] * std::exp(-a * r2);
        f += e;
        double fac = -2.0 * a * e;
        dfx += fac * dx; dfy += fac * dy; dfz += fac * dz;
    }
    if (sh.l == 0) {
        val = f;
        gx = dfx; gy = dfy; gz = dfz;
    } else {
        const double d[3] = {dx, dy, dz};
        const int c = sh.cart;
        val = d[c] * f;
        gx = d[c] * dfx + (c == 0 ? f : 0.0);
        gy = d[c] * dfy + (c == 1 ? f : 0.0);
        gz = d[c] * dfz + (c == 2 ? f : 0.0);
    }
}

} // namespace scf
} // namespace cpu
} // namespace ATV
