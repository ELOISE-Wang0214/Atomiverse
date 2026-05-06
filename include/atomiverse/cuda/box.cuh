/// @file cuda/box.cuh
/// @brief cuda::Box — canonical simulation cell + min-image + wrap on the GPU.
///
/// Plain trivially-copyable struct so the type can be uploaded to
/// `__constant__` memory or passed by value through a kernel argument.
/// Free `__host__ __device__` helpers form the parity API with `cpu::Box`.
///
/// **Migration note.** `cuda::Frame::box` is currently a device pointer
/// (`prec_t*`) holding a 9-element row-major cell. Switching that field
/// to a host-side `cuda::Box` value requires updating every kernel that
/// today reads `frame.box[i]`; that migration is tracked as a follow-up
/// to this spec. The struct + helpers below are introduced now so new
/// CUDA components can adopt the unified API and existing helpers in
/// `src/cuda/md/neighborlist/triclinic.cuh` can be folded into it.
///
/// @see include/atomiverse/cpu/box.h for the CPU mirror.
#pragma once

#include <type_traits>

#include "atomiverse/cuda/types.cuh"

namespace ATV {
namespace cuda {

/// @brief Plain device-uploadable cell + cached inverse + ortho fast-path scalars.
///
/// Row-major 3×3 storage matching `triclinic.cuh`. Trivially copyable
/// so the struct can travel through `cudaMemcpyToSymbol` or kernel-arg
/// pass-by-value without surgery.
struct Box {
    real_t cell[9]    = {};
    real_t inv_cell[9] = {};
    real_t half_l[3]  = {};
    real_t inv_l[3]   = {};
    bool   periodic[3] = {true, true, true};
    bool   ortho      = true;
};

static_assert(std::is_trivially_copyable<Box>::value,
              "cuda::Box must be trivially copyable for kernel-arg pass-by-value");

/// @brief Edge length along axis `d`.
__host__ __device__ inline real_t length(const Box& b, int d) noexcept {
    return b.cell[d * 3 + d];
}

/// @brief Cell volume.
///
/// Orthorhombic fast path: product of diagonal entries. Triclinic:
/// |det(cell)|.
__host__ __device__ inline real_t volume(const Box& b) noexcept {
    if (b.ortho) {
        return b.cell[0] * b.cell[4] * b.cell[8];
    }
    const real_t a = b.cell[0], bb = b.cell[1], c = b.cell[2];
    const real_t d = b.cell[3], e  = b.cell[4], f = b.cell[5];
    const real_t g = b.cell[6], h  = b.cell[7], i = b.cell[8];
    const real_t det = a * (e * i - f * h) - bb * (d * i - f * g) + c * (d * h - e * g);
    return det < real_t(0) ? -det : det;
}

/// @brief Apply minimum-image convention, in place.
///
/// Orthorhombic path uses the `nearbyint` form (branch-free, identical
/// across host/device); triclinic uses the inv_cell roundtrip. Matches
/// the pattern in `src/cuda/md/neighborlist/triclinic.cuh`.
template <typename T>
__host__ __device__ inline void
min_image(T& dx, T& dy, T& dz, const Box& b) noexcept {
    if (b.ortho) {
        dx -= T(b.cell[0]) * nearbyint(dx * T(b.inv_l[0]));
        dy -= T(b.cell[4]) * nearbyint(dy * T(b.inv_l[1]));
        dz -= T(b.cell[8]) * nearbyint(dz * T(b.inv_l[2]));
    } else {
        T fx = T(b.inv_cell[0]) * dx + T(b.inv_cell[1]) * dy + T(b.inv_cell[2]) * dz;
        T fy = T(b.inv_cell[3]) * dx + T(b.inv_cell[4]) * dy + T(b.inv_cell[5]) * dz;
        T fz = T(b.inv_cell[6]) * dx + T(b.inv_cell[7]) * dy + T(b.inv_cell[8]) * dz;
        fx -= nearbyint(fx);
        fy -= nearbyint(fy);
        fz -= nearbyint(fz);
        dx = T(b.cell[0]) * fx + T(b.cell[1]) * fy + T(b.cell[2]) * fz;
        dy = T(b.cell[3]) * fx + T(b.cell[4]) * fy + T(b.cell[5]) * fz;
        dz = T(b.cell[6]) * fx + T(b.cell[7]) * fy + T(b.cell[8]) * fz;
    }
}

/// @brief Wrap one axis into `[0, L_d)`. Orthorhombic only.
template <typename T>
__host__ __device__ inline T
wrap_axis(T q, int d, const Box& b) noexcept {
    const T L = T(b.cell[d * 3 + d]);
    if (L <= T(0)) return q;
    T w = q - L * floor(q * T(b.inv_l[d]));
    if (w < T(0))  w += L;
    if (w >= L)    w -= L;
    return w;
}

} // namespace cuda
} // namespace ATV
