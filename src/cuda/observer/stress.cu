/// @file observer/stress.cu
/// @brief VirialStress -- kinetic + pair virial stress tensor.

#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/observer/stress.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

struct VirialStress::State {
    prec_t* d_stress_tmp = nullptr; // [6] accumulation buffer
    int n_atoms = 0;
};

VirialStress::VirialStress() : state_(std::make_unique<State>()) {}
VirialStress::~VirialStress() = default;
VirialStress::VirialStress(VirialStress&&) noexcept = default;
VirialStress& VirialStress::operator=(VirialStress&&) noexcept = default;

namespace {

// Kinetic contribution with inline volume normalization.
// Each block's reduced sum is pre-divided by V = box[0]*box[4]*box[8]
// before the atomicAdd, eliminating the separate <<<1,1>>> volume kernel.
__global__ void kinetic_stress_k(const real_t* __restrict__ vx,
                                 const real_t* __restrict__ vy,
                                 const real_t* __restrict__ vz,
                                 const real_t* __restrict__ mass,
                                 prec_t* __restrict__ stress,
                                 const prec_t* __restrict__ box, int n) {
    extern __shared__ work_t sdata[]; // 6 * blockDim.x
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    work_t s[6] = {};
    if (idx < n) {
        work_t m = static_cast<work_t>(mass[idx]);
        work_t vxi = static_cast<work_t>(vx[idx]),
               vyi = static_cast<work_t>(vy[idx]),
               vzi = static_cast<work_t>(vz[idx]);
        s[0] = m * vxi * vxi; // xx
        s[1] = m * vyi * vyi; // yy
        s[2] = m * vzi * vzi; // zz
        s[3] = m * vxi * vyi; // xy
        s[4] = m * vxi * vzi; // xz
        s[5] = m * vyi * vzi; // yz
    }

    for (int c = 0; c < 6; ++c)
        sdata[c * blockDim.x + threadIdx.x] = s[c];
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride)
            for (int c = 0; c < 6; ++c)
                sdata[c * blockDim.x + threadIdx.x] +=
                    sdata[c * blockDim.x + threadIdx.x + stride];
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        // V = box[0][0] * box[1][1] * box[2][2] for orthorhombic
        work_t V_inv =
            work_t(1)
            / (static_cast<work_t>(box[0]) * static_cast<work_t>(box[4])
               * static_cast<work_t>(box[8]));
        for (int c = 0; c < 6; ++c)
            atomicAdd(&stress[c],
                      static_cast<prec_t>(sdata[c * blockDim.x] * V_inv));
    }
}

} // anonymous namespace

void VirialStress::build(const BuildContext& ctx) {
    state_->n_atoms = ctx.frame.n_atoms();
    cudaMalloc(&state_->d_stress_tmp, 6 * sizeof(prec_t));
}

void VirialStress::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    // Zero accumulation buffer
    cudaMemsetAsync(f.stress, 0, 6 * sizeof(prec_t), stream);

    // Kinetic contribution (volume factor folded into reduction)
    kinetic_stress_k<<<grid, block, 6 * block * sizeof(work_t), stream>>>(
        f.atoms.vx, f.atoms.vy, f.atoms.vz, f.atoms.mass, f.stress, f.box, n);
}

void VirialStress::destroy(const DestroyContext& ctx) {
    cudaFree(state_->d_stress_tmp);
    *state_ = {};
}

} // namespace ATV
