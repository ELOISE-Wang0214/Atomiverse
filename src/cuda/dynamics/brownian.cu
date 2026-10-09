/// @file src/cuda/dynamics/brownian.cu
/// @brief Device overdamped Brownian step (`brownian_integrate`).

#include "atomiverse/core/stage.h"
#include "atomiverse/cuda/device/launch.cuh"
#include "atomiverse/cuda/dynamics/brownian.cuh"
#include "atomiverse/cuda/dynamics/brownian_step.h"
#include "atomiverse/cuda/frame/types.cuh"

#include <cstdint>
#include <cuda_runtime.h>

namespace ATV {

using cuda::real_t;
using cuda::work_t;

namespace cuda {
namespace md {

namespace {

__global__ void brownian_integrate_kernel(
    real_t* __restrict__ x, real_t* __restrict__ y, real_t* __restrict__ z,
    const real_t* __restrict__ fx, const real_t* __restrict__ fy,
    const real_t* __restrict__ fz, const std::int64_t* __restrict__ tag,
    work_t mobility, work_t dt, work_t noise, work_t Lx, work_t Ly, work_t Lz,
    int wrap, std::uint64_t seed, int step, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;

    const std::uint64_t tag_i = tag ? static_cast<std::uint64_t>(tag[i])
                                    : static_cast<std::uint64_t>(i) + 1ull;
    work_t gx, gy, gz;
    brownian_gaussians(seed, static_cast<std::uint64_t>(step), tag_i, gx, gy,
                       gz);

    const work_t fxi = static_cast<work_t>(fx[i]);
    const work_t fyi = static_cast<work_t>(fy[i]);
    const work_t fzi = static_cast<work_t>(fz[i]);
    const bool finite = brownian_forces_finite(fxi, fyi, fzi);

    work_t nx = brownian_advance(static_cast<work_t>(x[i]), fxi, finite,
                                 mobility, dt, noise, gx);
    work_t ny = brownian_advance(static_cast<work_t>(y[i]), fyi, finite,
                                 mobility, dt, noise, gy);
    work_t nz = brownian_advance(static_cast<work_t>(z[i]), fzi, finite,
                                 mobility, dt, noise, gz);
    if (wrap) {
        nx = brownian_wrap(nx, Lx);
        ny = brownian_wrap(ny, Ly);
        nz = brownian_wrap(nz, Lz);
    }
    x[i] = static_cast<real_t>(nx);
    y[i] = static_cast<real_t>(ny);
    z[i] = static_cast<real_t>(nz);
}

} // namespace

void Brownian::launch(const LaunchContext& ctx) {
    if (!ctx.at(ATV::md::Phase::initial_integrate))
        return;

    auto& frame = ctx.frame;
    const int n = frame.n_atoms();
    const work_t dt = static_cast<work_t>(dt_);
    const work_t D = static_cast<work_t>(D_);
    const work_t temperature = static_cast<work_t>(temperature_);
    const work_t mobility = brownian_mobility(D, temperature);
    const work_t noise = brownian_noise_amplitude(D, dt);
    const cuda::Box& box = ctx.runtime.host_box;
    const int wrap =
        brownian_should_wrap(box.ortho, static_cast<double>(box.length(0))) ? 1
                                                                            : 0;
    const int block = cuda::kDefaultBlock;
    const int grid = cuda::grid_size(n);
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    brownian_integrate_kernel<<<grid, block, 0, stream>>>(
        frame.atoms.x, frame.atoms.y, frame.atoms.z, frame.atoms.fx,
        frame.atoms.fy, frame.atoms.fz, frame.atoms.tag, mobility, dt, noise,
        static_cast<work_t>(box.length(0)), static_cast<work_t>(box.length(1)),
        static_cast<work_t>(box.length(2)), wrap, seed_, ctx.step, n);
}

} // namespace md
} // namespace cuda
} // namespace ATV
