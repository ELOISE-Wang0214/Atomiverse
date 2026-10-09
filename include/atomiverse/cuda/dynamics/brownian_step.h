/// @file cuda/dynamics/brownian_step.h
/// @brief Overdamped Brownian displacement shared by the device kernel
///        and the host unit test.
///
///   r ← r + (D / kT) F Δt + √(2 D Δt) ξ,   ξ ~ N(0, 1)
///
/// Free particle (F = 0): per-coordinate variance is 2 D Δt, so
/// ⟨r²⟩ = 6 D t in 3D. Same update as `cpu::md::Brownian`.
#pragma once

#include <cmath>
#include <cstdint>

#include "atomiverse/core/units.h"

#if !defined(__CUDACC__)
#if !defined(__host__)
#define __host__
#endif
#if !defined(__device__)
#define __device__
#endif
#if !defined(__forceinline__)
#define __forceinline__ inline
#endif
#endif

namespace ATV {
namespace cuda {
namespace md {
namespace detail {

template <typename T>
__host__ __device__ __forceinline__ T bd_sqrt(T x) {
#if defined(__CUDA_ARCH__)
    return sqrt(x);
#else
    return std::sqrt(x);
#endif
}

template <typename T>
__host__ __device__ __forceinline__ T bd_log(T x) {
#if defined(__CUDA_ARCH__)
    return log(x);
#else
    return std::log(x);
#endif
}

template <typename T>
__host__ __device__ __forceinline__ T bd_cos(T x) {
#if defined(__CUDA_ARCH__)
    return cos(x);
#else
    return std::cos(x);
#endif
}

template <typename T>
__host__ __device__ __forceinline__ T bd_sin(T x) {
#if defined(__CUDA_ARCH__)
    return sin(x);
#else
    return std::sin(x);
#endif
}

template <typename T>
__host__ __device__ __forceinline__ T bd_floor(T x) {
#if defined(__CUDA_ARCH__)
    return floor(x);
#else
    return std::floor(x);
#endif
}

template <typename T>
__host__ __device__ __forceinline__ bool bd_finite(T v) {
#if defined(__CUDA_ARCH__)
    return isfinite(v);
#else
    return std::isfinite(static_cast<double>(v));
#endif
}

/// @brief splitmix64. Same mixer as `cuda::md::FixLangevin`.
__host__ __device__ __forceinline__ std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

/// @brief Map a hash into (0, 1]. Never zero, so `log` in Box–Muller is safe.
template <typename T>
__host__ __device__ __forceinline__ T u01(std::uint64_t h) {
    const std::uint64_t mant = (h >> 11) | std::uint64_t(1);
    return T(mant) * (T(1) / T(9007199254740992.0)); // 2^53
}

template <typename T>
__host__ __device__ __forceinline__ void box_muller(T u1, T u2, T& z1, T& z2) {
    const T r = bd_sqrt(T(-2) * bd_log(u1));
    const T theta = T(6.28318530717958647692) * u2;
    z1 = r * bd_cos(theta);
    z2 = r * bd_sin(theta);
}

} // namespace detail

/// @brief Mobility μ = D / kT [Å²/fs / eV].
/// @param D Diffusion constant [Å²/fs].
/// @param temperature Temperature [K].
/// @return μ such that μ F Δt is a length when F is eV/Å.
template <typename T>
__host__ __device__ __forceinline__ T brownian_mobility(T D, T temperature) {
    return D / (T(ATV::units::kB_ev_per_k) * temperature);
}

/// @brief Noise amplitude √(2 D Δt) [Å].
/// @param D Diffusion constant [Å²/fs].
/// @param dt Timestep [fs].
/// @return √(2 D Δt) [Å].
template <typename T>
__host__ __device__ __forceinline__ T brownian_noise_amplitude(T D, T dt) {
    return detail::bd_sqrt(T(2) * D * dt);
}

/// @brief True when every Cartesian force component is finite.
/// @param fx Force x [eV/Å].
/// @param fy Force y [eV/Å].
/// @param fz Force z [eV/Å].
/// @return False when any component is NaN or infinite.
template <typename T>
__host__ __device__ __forceinline__ bool brownian_forces_finite(T fx, T fy,
                                                                T fz) {
    return detail::bd_finite(fx) && detail::bd_finite(fy)
           && detail::bd_finite(fz);
}

/// @brief One coordinate of the overdamped step, before PBC.
/// @param x Position [Å].
/// @param f Force [eV/Å]. Ignored when @p force_finite is false.
/// @param force_finite Whether the force may enter the drift term.
/// @param mobility D/kT [Å²/fs / eV].
/// @param dt Timestep [fs].
/// @param noise √(2 D Δt) [Å].
/// @param xi Standard normal sample.
/// @return Updated position [Å], not wrapped.
template <typename T>
__host__ __device__ __forceinline__ T
brownian_advance(T x, T f, bool force_finite, T mobility, T dt, T noise, T xi) {
    T y = x;
    if (force_finite)
        y += mobility * f * dt;
    y += noise * xi;
    return y;
}

/// @brief Wrap one axis into [0, L). No-op for a free sandbox (L ≥ 1e5 Å)
///        or a non-positive edge, matching `cpu::md::Brownian`.
/// @param q Coordinate [Å].
/// @param L Edge length [Å].
/// @return Wrapped coordinate [Å].
template <typename T>
__host__ __device__ __forceinline__ T brownian_wrap(T q, T L) {
    if (!(L > T(0)) || !(L < T(1.0e5)))
        return q;
    return q - L * detail::bd_floor(q / L);
}

/// @brief True when the orthorhombic cell is small enough to wrap.
///        The CPU integrator keys this decision on `length(0)` only.
/// @param ortho Cell is orthorhombic.
/// @param length_x Edge length along x [Å].
/// @return False for a triclinic cell or a free sandbox (L ≥ 1e5 Å).
__host__ __device__ __forceinline__ bool brownian_should_wrap(bool ortho,
                                                              double length_x) {
    return ortho && length_x > 0.0 && length_x < 1.0e5;
}

/// @brief Three standard normals keyed by `(seed, step, tag)`.
///        Counter-based, so the device kernel has no RNG state.
/// @param seed User seed.
/// @param step MD step index.
/// @param tag Global atom tag (or a stand-in index).
/// @param gx Standard normal, x.
/// @param gy Standard normal, y.
/// @param gz Standard normal, z.
template <typename T>
__host__ __device__ __forceinline__ void
brownian_gaussians(std::uint64_t seed, std::uint64_t step, std::uint64_t tag,
                   T& gx, T& gy, T& gz) {
    const std::uint64_t base =
        seed * 0x9E3779B97F4A7C15ULL + step * 0xBF58476D1CE4E5B9ULL;
    const std::uint64_t key = base + tag * 0xD6E8FEB86659FD93ULL;
    const std::uint64_t h1 = detail::splitmix64(key);
    const std::uint64_t h2 = detail::splitmix64(h1);
    const std::uint64_t h3 = detail::splitmix64(h2);
    const std::uint64_t h4 = detail::splitmix64(h3);
    T z1, z2, z3, z_unused;
    detail::box_muller(detail::u01<T>(h1), detail::u01<T>(h2), z1, z2);
    detail::box_muller(detail::u01<T>(h3), detail::u01<T>(h4), z3, z_unused);
    gx = z1;
    gy = z2;
    gz = z3;
}

} // namespace md
} // namespace cuda
} // namespace ATV
