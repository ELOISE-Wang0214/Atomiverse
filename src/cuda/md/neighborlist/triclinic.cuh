#pragma once

/// @file triclinic.cuh
/// @brief Fractional coordinate transforms, minimum image, stencil computation.
///
/// Cell convention: row-major 3x3 matrix, r_cart = cell * r_frac.
/// Rows of cell are lattice vectors a, b, c.

#include <cmath>
#include <cstdint>

namespace ATV::md::nlist {

/// @brief 3-component vector for device-side math.
template <typename Real>
struct Vec3 {
    Real x, y, z;
};

/// @brief Compute the inverse of a 3x3 row-major matrix.
/// @param m Input 3x3 matrix (row-major)
/// @param inv Output inverse matrix (row-major)
/// @return determinant
template <typename Real>
__host__ __device__ inline Real invert_3x3(const Real* m, Real* inv) {
    Real det = m[0] * (m[4] * m[8] - m[5] * m[7])
             - m[1] * (m[3] * m[8] - m[5] * m[6])
             + m[2] * (m[3] * m[7] - m[4] * m[6]);

    Real inv_det = Real(1) / det;
    inv[0] = (m[4] * m[8] - m[5] * m[7]) * inv_det;
    inv[1] = (m[2] * m[7] - m[1] * m[8]) * inv_det;
    inv[2] = (m[1] * m[5] - m[2] * m[4]) * inv_det;
    inv[3] = (m[5] * m[6] - m[3] * m[8]) * inv_det;
    inv[4] = (m[0] * m[8] - m[2] * m[6]) * inv_det;
    inv[5] = (m[2] * m[3] - m[0] * m[5]) * inv_det;
    inv[6] = (m[3] * m[7] - m[4] * m[6]) * inv_det;
    inv[7] = (m[1] * m[6] - m[0] * m[7]) * inv_det;
    inv[8] = (m[0] * m[4] - m[1] * m[3]) * inv_det;
    return det;
}

/// @brief Cartesian to fractional: r_frac = inv_cell * r_cart
/// where cell is row-major, inv_cell is inverse of cell.
template <typename Real>
__host__ __device__ inline Vec3<Real> cart_to_frac(Real rx, Real ry, Real rz,
                                                    const Real* inv_cell) {
    return {inv_cell[0] * rx + inv_cell[1] * ry + inv_cell[2] * rz,
            inv_cell[3] * rx + inv_cell[4] * ry + inv_cell[5] * rz,
            inv_cell[6] * rx + inv_cell[7] * ry + inv_cell[8] * rz};
}

/// @brief Fractional to Cartesian: r_cart = cell * r_frac
template <typename Real>
__host__ __device__ inline Vec3<Real> frac_to_cart(Real fx, Real fy, Real fz,
                                                    const Real* cell) {
    return {cell[0] * fx + cell[1] * fy + cell[2] * fz,
            cell[3] * fx + cell[4] * fy + cell[5] * fz,
            cell[6] * fx + cell[7] * fy + cell[8] * fz};
}

/// @brief Clamp fractional coordinate to [0, 1).
/// @brief Boundary epsilon for clamping fractional coordinates to [0,1).
template <typename Real>
inline constexpr Real kFracBoundaryEps = Real(1e-6);

template <typename Real>
__device__ inline Real clamp01(Real v) {
    if (v < Real(0)) return Real(0);
    if (v >= Real(1)) return Real(1) - kFracBoundaryEps<Real>;
    return v;
}

/// @brief Apply minimum image convention in fractional space.
/// @param[in,out] dx, dy, dz Cartesian displacement, modified in place
template <typename Real>
__host__ __device__ inline void minimum_image(Real& dx, Real& dy, Real& dz,
                                               const Real* cell,
                                               const Real* inv_cell) {
    // To fractional
    Real fx = inv_cell[0] * dx + inv_cell[1] * dy + inv_cell[2] * dz;
    Real fy = inv_cell[3] * dx + inv_cell[4] * dy + inv_cell[5] * dz;
    Real fz = inv_cell[6] * dx + inv_cell[7] * dy + inv_cell[8] * dz;

    // Nearest image
    fx -= nearbyint(fx);
    fy -= nearbyint(fy);
    fz -= nearbyint(fz);

    // Back to Cartesian
    dx = cell[0] * fx + cell[1] * fy + cell[2] * fz;
    dy = cell[3] * fx + cell[4] * fy + cell[5] * fz;
    dz = cell[6] * fx + cell[7] * fy + cell[8] * fz;
}

/// @brief Wrap fractional coordinate to [0, 1).
template <typename Real>
__host__ __device__ inline Real wrap_frac(Real f) {
    f -= floor(f);
    // Guard against numerical edge case where floor gives exactly 1.0
    if (f >= Real(1)) f -= Real(1);
    if (f < Real(0)) f += Real(1);
    return f;
}

} // namespace ATV::md::nlist
