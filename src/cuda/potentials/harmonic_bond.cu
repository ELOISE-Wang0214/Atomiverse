/// @file bonded/harmonic_bond.cu
/// @brief HarmonicBond implementation -- V(r) = K(r - r0)^2.

#include "atomiverse/cuda/potentials/harmonic_bond.cuh"
#include "atomiverse/cuda/potentials/harmonic_bond_kernel.cuh"
#include "atomiverse/cuda/launch.cuh"
#include <cuda_runtime.h>
#include "atomiverse/cuda/types.cuh"
#include <cmath>

namespace ATV {

using cuda::real_t;
using cuda::work_t;
using cuda::prec_t;


namespace cuda { namespace kernels {

__global__ void harmonic_bond_k(
    const real_t* __restrict__ x, const real_t* __restrict__ y, const real_t* __restrict__ z,
    real_t* __restrict__ fx, real_t* __restrict__ fy, real_t* __restrict__ fz,
    real_t* __restrict__ energy,
    const int* __restrict__ bi, const int* __restrict__ bj, const int* __restrict__ bt,
    const real_t* __restrict__ K, const real_t* __restrict__ r0,
    const real_t* __restrict__ box, int n_bonds)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_bonds) return;

    int i = bi[idx], j = bj[idx], t = bt[idx];

    work_t dx = x[j] - x[i];
    work_t dy = y[j] - y[i];
    work_t dz = z[j] - z[i];

    // Minimum image convention (orthorhombic)
    work_t bx = box[0], by = box[4], bz = box[8];
    dx -= bx * round(dx / bx);
    dy -= by * round(dy / by);
    dz -= bz * round(dz / bz);

    work_t r = sqrt(dx * dx + dy * dy + dz * dz);
    work_t dr = r - r0[t];

    // V = K(r - r0)^2
    work_t v = K[t] * dr * dr;

    // F = -2K(r - r0) * r_hat  (on atom j; opposite on atom i)
    work_t f_over_r = work_t(-2) * K[t] * dr / r;
    work_t f_x = f_over_r * dx;
    work_t f_y = f_over_r * dy;
    work_t f_z = f_over_r * dz;

    atomicAdd(&fx[j], f_x);
    atomicAdd(&fy[j], f_y);
    atomicAdd(&fz[j], f_z);
    atomicAdd(&fx[i], -f_x);
    atomicAdd(&fy[i], -f_y);
    atomicAdd(&fz[i], -f_z);

    // Block-level energy reduction
    extern __shared__ work_t sdata[];
    sdata[threadIdx.x] = v;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sdata[threadIdx.x] += sdata[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) atomicAdd(energy, sdata[0]);
}

}} // namespace cuda::kernels

struct HarmonicBond::State {
    Real* d_K = nullptr;
    Real* d_r0 = nullptr;
    int n_types = 0;
};

HarmonicBond::HarmonicBond(std::vector<work_t> K, std::vector<work_t> r0)
    : K_(std::move(K)), r0_(std::move(r0)), state_(std::make_unique<State>()) {}

HarmonicBond::~HarmonicBond() = default;
HarmonicBond::HarmonicBond(HarmonicBond&&) noexcept = default;
HarmonicBond& HarmonicBond::operator=(HarmonicBond&&) noexcept = default;

void HarmonicBond::build(const BuildContext& ctx) {
    int nt = static_cast<int>(K_.size());
    state_->n_types = nt;
    cudaMalloc(&state_->d_K, nt * sizeof(work_t));
    cudaMalloc(&state_->d_r0, nt * sizeof(work_t));
    cudaMemcpy(state_->d_K, K_.data(), nt * sizeof(work_t), cudaMemcpyHostToDevice);
    cudaMemcpy(state_->d_r0, r0_.data(), nt * sizeof(work_t), cudaMemcpyHostToDevice);
}

void HarmonicBond::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int nb = f.n_bonds();
    if (nb == 0) return;

    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(nb));
    cuda::kernels::harmonic_bond_k<<<grid, block, block * sizeof(work_t), stream>>>(
        f.atoms.x, f.atoms.y, f.atoms.z, f.atoms.fx, f.atoms.fy, f.atoms.fz, f.energy,
        f.bonds.atomi, f.bonds.atomj, f.bonds.type,
        state_->d_K, state_->d_r0, f.box, nb);
}

void HarmonicBond::destroy(const DestroyContext& ctx) {
    cudaFree(state_->d_K);
    cudaFree(state_->d_r0);
    *state_ = {};
}

} // namespace ATV
