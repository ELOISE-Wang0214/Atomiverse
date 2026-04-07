/// @file minimizer/sd.cu
/// @brief SteepestDescentStep kernel: x += alpha * f / m * kFtm2v.

#include "atomiverse/cuda/minimizer/sd.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/core/units.h"
#include <cuda_runtime.h>
#include "atomiverse/cuda/types.cuh"

namespace ATV {

using cuda::real_t;
using cuda::work_t;
using cuda::prec_t;

namespace minimizer {

namespace {

/// @brief Steepest descent position update.
///
/// Moves each atom along the force direction:
///   x[i] += alpha * fx[i] / mass[i] * kFtm2v
///
/// Same unit conversion as VelocityVerlet. At step 0, forces are zero
/// (from initialization), so positions remain unchanged.
__global__ void sd_step(
    real_t* __restrict__ x, real_t* __restrict__ y, real_t* __restrict__ z,
    const real_t* __restrict__ fx, const real_t* __restrict__ fy,
    const real_t* __restrict__ fz, const real_t* __restrict__ mass,
    const prec_t* __restrict__ box,
    work_t alpha, int n) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    work_t inv_m = work_t(1) / mass[i];
    x[i] += alpha * fx[i] * inv_m * kFtm2v;
    y[i] += alpha * fy[i] * inv_m * kFtm2v;
    z[i] += alpha * fz[i] * inv_m * kFtm2v;
    // PBC wrapping
    work_t Lx = box[0], Ly = box[4], Lz = box[8];
    work_t inv_Lx = work_t(1) / Lx, inv_Ly = work_t(1) / Ly, inv_Lz = work_t(1) / Lz;
    x[i] -= Lx * floor(x[i] * inv_Lx);
    y[i] -= Ly * floor(y[i] * inv_Ly);
    z[i] -= Lz * floor(z[i] * inv_Lz);
}

} // anonymous namespace

void SteepestDescentStep::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    int n = ctx.frame.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));
    sd_step<<<grid, block, 0, stream>>>(
        frame.atoms.x, frame.atoms.y, frame.atoms.z,
        frame.atoms.fx, frame.atoms.fy, frame.atoms.fz, frame.atoms.mass,
        frame.box, alpha_, n);
}

} // namespace minimizer
} // namespace ATV
