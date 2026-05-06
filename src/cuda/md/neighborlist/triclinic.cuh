#pragma once

/// @file triclinic.cuh
/// @brief Fractional coordinate transforms, minimum image, stencil computation.
///
/// Cell convention: row-major 3x3 matrix, r_cart = cell * r_frac.
/// Rows of cell are lattice vectors a, b, c.
///
/// **Migration note.** `ATV::cuda::Box` (in
/// `include/atomiverse/cuda/box.cuh`) is the canonical Box at the CUDA
/// backend root. New callers should use `ATV::cuda::min_image(dx, dy,
/// dz, box)` directly instead of the raw-pointer overloads here.
/// The pointer overloads stay so existing kernel signatures
/// (`(cell, inv_cell)` pairs) compile unchanged during the staged
/// migration.

#include <cmath>
#include <cstdint>

#include "atomiverse/cuda/box.cuh"

namespace ATV::md::nlist {

/// @brief Forward to `ATV::cuda::min_image(dx, dy, dz, box)`.
///
/// Box-aware overload of the older pointer-based `minimum_image`. New
/// code should call this; the pointer-based templated overload below
/// remains for compatibility with kernel signatures that still pass
/// `(cell, inv_cell)` directly.
template <typename work_t>
__host__ __device__ inline void
minimum_image(work_t& dx, work_t& dy, work_t& dz,
              const ATV::cuda::Box& box) {
    ATV::cuda::min_image(dx, dy, dz, box);
}

/// @brief 3-component vector for device-side math.
template <typename work_t>
struct Vec3 {
    work_t x, y, z;
};

/// @brief Compute the inverse of a 3x3 row-major matrix.
/// @param m Input 3x3 matrix (row-major)
/// @param inv Output inverse matrix (row-major)
/// @return determinant
template <typename work_t>
__host__ __device__ inline work_t invert_3x3(const real_t* m, work_t* inv) {
    work_t det = m[0] * (m[4] * m[8] - m[5] * m[7])
             - m[1] * (m[3] * m[8] - m[5] * m[6])
             + m[2] * (m[3] * m[7] - m[4] * m[6]);

    work_t inv_det = work_t(1) / det;
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
template <typename work_t>
__host__ __device__ inline Vec3<work_t> cart_to_frac(work_t rx, work_t ry, work_t rz,
                                                    const real_t* inv_cell) {
    return {inv_cell[0] * rx + inv_cell[1] * ry + inv_cell[2] * rz,
            inv_cell[3] * rx + inv_cell[4] * ry + inv_cell[5] * rz,
            inv_cell[6] * rx + inv_cell[7] * ry + inv_cell[8] * rz};
}

/// @brief Fractional to Cartesian: r_cart = cell * r_frac
template <typename work_t>
__host__ __device__ inline Vec3<work_t> frac_to_cart(work_t fx, work_t fy, work_t fz,
                                                    const real_t* cell) {
    return {cell[0] * fx + cell[1] * fy + cell[2] * fz,
            cell[3] * fx + cell[4] * fy + cell[5] * fz,
            cell[6] * fx + cell[7] * fy + cell[8] * fz};
}

/// @brief Clamp fractional coordinate to [0, 1).
/// @brief Boundary epsilon for clamping fractional coordinates to [0,1).
template <typename work_t>
inline constexpr work_t kFracBoundaryEps = work_t(1e-6);

template <typename work_t>
__device__ inline work_t clamp01(work_t v) {
    if (v < work_t(0)) return work_t(0);
    if (v >= work_t(1)) return work_t(1) - kFracBoundaryEps<work_t>;
    return v;
}

/// @brief Apply minimum image convention in fractional space.
/// @param[in,out] dx, dy, dz Cartesian displacement, modified in place
template <typename work_t>
__host__ __device__ inline void minimum_image(work_t& dx, work_t& dy, work_t& dz,
                                               const real_t* cell,
                                               const real_t* inv_cell) {
    // To fractional
    work_t fx = inv_cell[0] * dx + inv_cell[1] * dy + inv_cell[2] * dz;
    work_t fy = inv_cell[3] * dx + inv_cell[4] * dy + inv_cell[5] * dz;
    work_t fz = inv_cell[6] * dx + inv_cell[7] * dy + inv_cell[8] * dz;

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
template <typename work_t>
__host__ __device__ inline work_t wrap_frac(work_t f) {
    f -= floor(f);
    // Guard against numerical edge case where floor gives exactly 1.0
    if (f >= work_t(1)) f -= work_t(1);
    if (f < work_t(0)) f += work_t(1);
    return f;
}

} // namespace ATV::md::nlist
