/// @file cuda/potentials/pair_lj_cut.cu
/// @brief CUDA PairLJCut — Lennard-Jones 12-6 shifted (kernel + class).

#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/potentials/lj_kernel.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

namespace cuda {
namespace kernels {

__global__ void lj_force(real_t* __restrict__ fx, real_t* __restrict__ fy,
                         real_t* __restrict__ fz, prec_t* __restrict__ energy,
                         const int* __restrict__ pi, const int* __restrict__ pj,
                         const real_t* __restrict__ pdx,
                         const real_t* __restrict__ pdy,
                         const real_t* __restrict__ pdz,
                         const real_t* __restrict__ pdist, work_t epsilon,
                         work_t sigma, work_t cutoff, work_t v_shift,
                         const int* __restrict__ n_pairs_ptr) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int n_pairs = *n_pairs_ptr;

    work_t my_energy = work_t(0);
    if (idx < n_pairs) {
        work_t r = pdist[idx];
        if (r < cutoff && r >= kMinPairDist) {
            work_t inv_r = work_t(1) / r;
            work_t sr = sigma * inv_r;
            work_t sr2 = sr * sr;
            work_t sr6 = sr2 * sr2 * sr2;
            work_t sr12 = sr6 * sr6;

            work_t f_over_r =
                work_t(24) * epsilon * (work_t(2) * sr12 - sr6) * inv_r * inv_r;
            work_t dx = pdx[idx], dy = pdy[idx], dz = pdz[idx];

            int i = pi[idx], j = pj[idx];
            atomicAdd(&fx[i], static_cast<real_t>(-f_over_r * dx));
            atomicAdd(&fy[i], static_cast<real_t>(-f_over_r * dy));
            atomicAdd(&fz[i], static_cast<real_t>(-f_over_r * dz));
            atomicAdd(&fx[j], static_cast<real_t>(f_over_r * dx));
            atomicAdd(&fy[j], static_cast<real_t>(f_over_r * dy));
            atomicAdd(&fz[j], static_cast<real_t>(f_over_r * dz));

            my_energy = work_t(4) * epsilon * (sr12 - sr6) - v_shift;
        }
    }

    extern __shared__ work_t s_energy[];
    s_energy[threadIdx.x] = my_energy;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s)
            s_energy[threadIdx.x] += s_energy[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0 && s_energy[0] != work_t(0))
        atomicAdd(energy, static_cast<prec_t>(s_energy[0]));
}

} // namespace kernels
} // namespace cuda

namespace {
/// @brief Estimated max neighbors per atom for grid-sizing heuristic.
constexpr int kEstMaxNeighbors = 100;
} // anonymous namespace

namespace cuda {

PairLJCut::PairLJCut(double epsilon, double sigma, double cutoff)
    : Pair(cutoff), epsilon_(epsilon), sigma_(sigma) {
    v_shift_ = kernels::lj_shift(epsilon, sigma, cutoff);
}

void PairLJCut::compute(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    int n = ctx.frame.n_atoms();
    int max_pairs = static_cast<int>(n) * kEstMaxNeighbors / 2;
    int block = kDefaultBlock;
    int grid = grid_size(max_pairs);
    kernels::lj_force<<<grid, block, block * sizeof(work_t), stream>>>(
        frame.atoms.fx, frame.atoms.fy, frame.atoms.fz, frame.energy,
        frame.pairs.i, frame.pairs.j, frame.pairs.dx, frame.pairs.dy,
        frame.pairs.dz, frame.pairs.dist, epsilon_, sigma_, cutoff_, v_shift_,
        frame.pairs.n_pairs);
}

} // namespace cuda
} // namespace ATV
