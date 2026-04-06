/// @file md/kinetic_energy.cu
/// @brief KineticEnergy kernel: per-atom KE with shared-memory block reduction.

#include "atomiverse/cuda/md/kinetic_energy.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/core/units.h"
#include <cuda_runtime.h>

namespace ATV {
namespace md {

namespace {

/// @brief Per-atom kinetic energy with block-level reduction.
///
/// Each thread computes ke_i = 0.5 * kMvv2e * mass[i] * (vx^2 + vy^2 + vz^2).
/// Block-level shared-memory reduction sums contributions.
/// Thread 0 of each block atomicAdds the block sum to the global accumulator.
__global__ void kinetic_energy(
    const Real* __restrict__ vx, const Real* __restrict__ vy,
    const Real* __restrict__ vz, const Real* __restrict__ mass,
    Real* __restrict__ kinetic, Int n) {

    Int i = blockIdx.x * blockDim.x + threadIdx.x;

    Real my_ke = Real(0);
    if (i < n) {
        Real v2 = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        my_ke = Real(0.5) * kMvv2e * mass[i] * v2;
    }

    extern __shared__ Real s_ke[];
    s_ke[threadIdx.x] = my_ke;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) s_ke[threadIdx.x] += s_ke[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0 && s_ke[0] != Real(0))
        atomicAdd(kinetic, s_ke[0]);
}

} // anonymous namespace

void KineticEnergy::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    Int n = ctx.frame.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));
    kinetic_energy<<<grid, block, block * sizeof(Real), stream>>>(
        frame.atoms.vx, frame.atoms.vy, frame.atoms.vz, frame.atoms.mass,
        frame.kinetic, n);
}

} // namespace md
} // namespace ATV
