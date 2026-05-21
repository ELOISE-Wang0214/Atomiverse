/// @file minimizer/sd.cu
/// @brief SteepestDescentStep kernel: x += alpha * f / m * kFtm2v.

#include "atomiverse/core/units.h"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/minimizer/sd.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

namespace minimizer {

namespace {

/// @brief Steepest descent position update.
///
/// Moves each atom along the force direction:
///   x[i] += alpha * fx[i] / mass[i] * kFtm2v
///
/// Same unit conversion as VelocityVerlet. At step 0, forces are zero
/// (from initialization), so positions remain unchanged.
__global__ void sd_step(real_t* __restrict__ x, real_t* __restrict__ y,
                        real_t* __restrict__ z, const real_t* __restrict__ fx,
                        const real_t* __restrict__ fy,
                        const real_t* __restrict__ fz,
                        const real_t* __restrict__ mass,
                        const prec_t* __restrict__ box, work_t alpha, int n) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;

    work_t inv_m = work_t(1) / static_cast<work_t>(mass[i]);
    x[i] = static_cast<real_t>(static_cast<work_t>(x[i])
                               + alpha * static_cast<work_t>(fx[i]) * inv_m
                                     * kFtm2v);
    y[i] = static_cast<real_t>(static_cast<work_t>(y[i])
                               + alpha * static_cast<work_t>(fy[i]) * inv_m
                                     * kFtm2v);
    z[i] = static_cast<real_t>(static_cast<work_t>(z[i])
                               + alpha * static_cast<work_t>(fz[i]) * inv_m
                                     * kFtm2v);
    // PBC wrapping
    work_t Lx = static_cast<work_t>(box[0]), Ly = static_cast<work_t>(box[4]),
           Lz = static_cast<work_t>(box[8]);
    work_t inv_Lx = work_t(1) / Lx, inv_Ly = work_t(1) / Ly,
           inv_Lz = work_t(1) / Lz;
    work_t x_w = static_cast<work_t>(x[i]);
    work_t y_w = static_cast<work_t>(y[i]);
    work_t z_w = static_cast<work_t>(z[i]);
    x[i] = static_cast<real_t>(x_w - Lx * floor(x_w * inv_Lx));
    y[i] = static_cast<real_t>(y_w - Ly * floor(y_w * inv_Ly));
    z[i] = static_cast<real_t>(z_w - Lz * floor(z_w * inv_Lz));
}

} // anonymous namespace

void SteepestDescentStep::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    int n = ctx.frame.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));
    sd_step<<<grid, block, 0, stream>>>(
        frame.atoms.x, frame.atoms.y, frame.atoms.z, frame.atoms.fx,
        frame.atoms.fy, frame.atoms.fz, frame.atoms.mass, frame.box, alpha_, n);
}

} // namespace minimizer
} // namespace ATV
