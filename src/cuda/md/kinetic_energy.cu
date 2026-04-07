/// @file md/kinetic_energy.cu
/// @brief KineticEnergy kernel: per-atom KE with shared-memory block reduction.

#include "atomiverse/cuda/md/kinetic_energy.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/core/units.h"
#include <cuda_runtime.h>
#include "atomiverse/cuda/types.cuh"

namespace ATV {

using cuda::real_t;
using cuda::work_t;
using cuda::prec_t;

namespace md {

namespace {

/// @brief Per-atom kinetic energy with block-level reduction.
///
/// Each thread computes ke_i = 0.5 * kMvv2e * mass[i] * (vx^2 + vy^2 + vz^2).
/// Block-level shared-memory reduction sums contributions.
/// Thread 0 of each block atomicAdds the block sum to the global accumulator.
__global__ void kinetic_energy(
    const real_t* __restrict__ vx, const real_t* __restrict__ vy,
    const real_t* __restrict__ vz, const real_t* __restrict__ mass,
    real_t* __restrict__ kinetic, int n) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    work_t my_ke = work_t(0);
    if (i < n) {
        work_t v2 = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        my_ke = work_t(0.5) * kMvv2e * mass[i] * v2;
    }

    extern __shared__ work_t s_ke[];
    s_ke[threadIdx.x] = my_ke;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) s_ke[threadIdx.x] += s_ke[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0 && s_ke[0] != work_t(0))
        atomicAdd(kinetic, s_ke[0]);
}

} // anonymous namespace

void KineticEnergy::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    int n = ctx.frame.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));
    kinetic_energy<<<grid, block, block * sizeof(work_t), stream>>>(
        frame.atoms.vx, frame.atoms.vy, frame.atoms.vz, frame.atoms.mass,
        frame.kinetic, n);
}

} // namespace md
} // namespace ATV
