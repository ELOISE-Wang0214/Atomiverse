/// @file minimizer/sd.cu
/// @brief SteepestDescentStep kernel: x += alpha * f / m * kFtm2v.

#include "atomiverse/cuda/minimizer/sd.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/core/units.h"
#include <cuda_runtime.h>

namespace ATV {
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
    Real* __restrict__ x, Real* __restrict__ y, Real* __restrict__ z,
    const Real* __restrict__ fx, const Real* __restrict__ fy,
    const Real* __restrict__ fz, const Real* __restrict__ mass,
    const Real* __restrict__ box,
    Real alpha, Int n) {

    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    Real inv_m = Real(1) / mass[i];
    x[i] += alpha * fx[i] * inv_m * kFtm2v;
    y[i] += alpha * fy[i] * inv_m * kFtm2v;
    z[i] += alpha * fz[i] * inv_m * kFtm2v;
    // PBC wrapping
    Real Lx = box[0], Ly = box[4], Lz = box[8];
    Real inv_Lx = Real(1) / Lx, inv_Ly = Real(1) / Ly, inv_Lz = Real(1) / Lz;
    x[i] -= Lx * floor(x[i] * inv_Lx);
    y[i] -= Ly * floor(y[i] * inv_Ly);
    z[i] -= Lz * floor(z[i] * inv_Lz);
}

} // anonymous namespace

void SteepestDescentStep::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    Int n = ctx.frame.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));
    sd_step<<<grid, block, 0, stream>>>(
        frame.x, frame.y, frame.z,
        frame.fx, frame.fy, frame.fz, frame.mass,
        frame.box, alpha_, n);
}

} // namespace minimizer
} // namespace ATV
