/// @file observer/stress.cu
/// @brief VirialStress -- kinetic + pair virial stress tensor.

#include "atomiverse/cuda/observer/stress.cuh"
#include "atomiverse/cuda/launch.cuh"
#include <cuda_runtime.h>

namespace ATV {

struct VirialStress::State {
    Real* d_stress_tmp = nullptr;  // [6] accumulation buffer
    Int n_atoms = 0;
};

VirialStress::VirialStress() : state_(std::make_unique<State>()) {}
VirialStress::~VirialStress() = default;
VirialStress::VirialStress(VirialStress&&) noexcept = default;
VirialStress& VirialStress::operator=(VirialStress&&) noexcept = default;

namespace {

// Kinetic contribution with inline volume normalization.
// Each block's reduced sum is pre-divided by V = box[0]*box[4]*box[8]
// before the atomicAdd, eliminating the separate <<<1,1>>> volume kernel.
__global__ void kinetic_stress_k(
    const Real* __restrict__ vx, const Real* __restrict__ vy, const Real* __restrict__ vz,
    const Real* __restrict__ mass, Real* __restrict__ stress,
    const Real* __restrict__ box, Int n)
{
    extern __shared__ Real sdata[];  // 6 * blockDim.x
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;

    Real s[6] = {};
    if (idx < n) {
        Real m = mass[idx];
        Real vxi = vx[idx], vyi = vy[idx], vzi = vz[idx];
        s[0] = m * vxi * vxi;  // xx
        s[1] = m * vyi * vyi;  // yy
        s[2] = m * vzi * vzi;  // zz
        s[3] = m * vxi * vyi;  // xy
        s[4] = m * vxi * vzi;  // xz
        s[5] = m * vyi * vzi;  // yz
    }

    for (int c = 0; c < 6; ++c) sdata[c * blockDim.x + threadIdx.x] = s[c];
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride)
            for (int c = 0; c < 6; ++c)
                sdata[c * blockDim.x + threadIdx.x] += sdata[c * blockDim.x + threadIdx.x + stride];
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        // V = box[0][0] * box[1][1] * box[2][2] for orthorhombic
        Real V_inv = Real(1) / (box[0] * box[4] * box[8]);
        for (int c = 0; c < 6; ++c)
            atomicAdd(&stress[c], sdata[c * blockDim.x] * V_inv);
    }
}

} // anonymous namespace

void VirialStress::build(const BuildContext& ctx) {
    state_->n_atoms = ctx.frame.n_atoms();
    cudaMalloc(&state_->d_stress_tmp, 6 * sizeof(Real));
}

void VirialStress::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    Int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    // Zero accumulation buffer
    cudaMemsetAsync(f.stress, 0, 6 * sizeof(Real), stream);

    // Kinetic contribution (volume factor folded into reduction)
    kinetic_stress_k<<<grid, block, 6 * block * sizeof(Real), stream>>>(
        f.vx, f.vy, f.vz, f.mass, f.stress, f.box, n);
}

void VirialStress::destroy(const DestroyContext& ctx) {
    cudaFree(state_->d_stress_tmp);
    *state_ = {};
}

} // namespace ATV
