/// @file md/lj.cu

#include "atomiverse/cuda/potentials/lj.cuh"
#include "atomiverse/cuda/potentials/lj_kernel.cuh"
#include "atomiverse/cuda/launch.cuh"
#include <cuda_runtime.h>

namespace ATV {

namespace cuda { namespace kernels {

__global__ void lj_force(
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ pi, const Int* __restrict__ pj,
    const Real* __restrict__ pdx, const Real* __restrict__ pdy,
    const Real* __restrict__ pdz, const Real* __restrict__ pdist,
    Real epsilon, Real sigma, Real cutoff, Real v_shift,
    const Int* __restrict__ n_pairs_ptr) {

    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    Int n_pairs = *n_pairs_ptr;

    Real my_energy = Real(0);
    if (idx < n_pairs) {
        Real r = pdist[idx];
        if (r < cutoff && r >= kMinPairDist) {
            Real inv_r = Real(1) / r;
            Real sr = sigma * inv_r;
            Real sr2 = sr * sr;
            Real sr6 = sr2 * sr2 * sr2;
            Real sr12 = sr6 * sr6;

            Real f_over_r = Real(24) * epsilon * (Real(2) * sr12 - sr6) * inv_r * inv_r;
            Real dx = pdx[idx], dy = pdy[idx], dz = pdz[idx];

            Int i = pi[idx], j = pj[idx];
            atomicAdd(&fx[i], -f_over_r * dx);
            atomicAdd(&fy[i], -f_over_r * dy);
            atomicAdd(&fz[i], -f_over_r * dz);
            atomicAdd(&fx[j],  f_over_r * dx);
            atomicAdd(&fy[j],  f_over_r * dy);
            atomicAdd(&fz[j],  f_over_r * dz);

            my_energy = Real(4) * epsilon * (sr12 - sr6) - v_shift;
        }
    }

    extern __shared__ Real s_energy[];
    s_energy[threadIdx.x] = my_energy;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) s_energy[threadIdx.x] += s_energy[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0 && s_energy[0] != Real(0))
        atomicAdd(energy, s_energy[0]);
}

}} // namespace cuda::kernels

namespace {
/// @brief Estimated max neighbors per atom for grid-sizing heuristic.
constexpr int kEstMaxNeighbors = 100;
} // anonymous namespace

LJ::LJ(Real epsilon, Real sigma, Real cutoff)
    : epsilon_(epsilon), sigma_(sigma), cutoff_(cutoff) {
    v_shift_ = cuda::kernels::lj_shift(epsilon, sigma, cutoff);
}

void LJ::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    Int n = ctx.frame.n_atoms();
    int max_pairs = static_cast<int>(n) * kEstMaxNeighbors / 2;
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(max_pairs);
    cuda::kernels::lj_force<<<grid, block, block * sizeof(Real), stream>>>(
        frame.fx, frame.fy, frame.fz, frame.energy,
        frame.pairs.i, frame.pairs.j, frame.pairs.dx, frame.pairs.dy,
        frame.pairs.dz, frame.pairs.dist,
        epsilon_, sigma_, cutoff_, v_shift_, frame.pairs.n_pairs);
}

} // namespace ATV
