/// @file bonded/harmonic_bond.cu
/// @brief HarmonicBond implementation -- V(r) = K(r - r0)^2.

#include "atomiverse/cuda/potentials/harmonic_bond.cuh"
#include "atomiverse/cuda/potentials/harmonic_bond_kernel.cuh"
#include "atomiverse/cuda/launch.cuh"
#include <cuda_runtime.h>
#include <cmath>

namespace ATV {

namespace cuda { namespace kernels {

__global__ void harmonic_bond_k(
    const Real* __restrict__ x, const Real* __restrict__ y, const Real* __restrict__ z,
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ bi, const Int* __restrict__ bj, const Int* __restrict__ bt,
    const Real* __restrict__ K, const Real* __restrict__ r0,
    const Real* __restrict__ box, Int n_bonds)
{
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_bonds) return;

    Int i = bi[idx], j = bj[idx], t = bt[idx];

    Real dx = x[j] - x[i];
    Real dy = y[j] - y[i];
    Real dz = z[j] - z[i];

    // Minimum image convention (orthorhombic)
    Real bx = box[0], by = box[4], bz = box[8];
    dx -= bx * round(dx / bx);
    dy -= by * round(dy / by);
    dz -= bz * round(dz / bz);

    Real r = sqrt(dx * dx + dy * dy + dz * dz);
    Real dr = r - r0[t];

    // V = K(r - r0)^2
    Real v = K[t] * dr * dr;

    // F = -2K(r - r0) * r_hat  (on atom j; opposite on atom i)
    Real f_over_r = Real(-2) * K[t] * dr / r;
    Real f_x = f_over_r * dx;
    Real f_y = f_over_r * dy;
    Real f_z = f_over_r * dz;

    atomicAdd(&fx[j], f_x);
    atomicAdd(&fy[j], f_y);
    atomicAdd(&fz[j], f_z);
    atomicAdd(&fx[i], -f_x);
    atomicAdd(&fy[i], -f_y);
    atomicAdd(&fz[i], -f_z);

    // Block-level energy reduction
    extern __shared__ Real sdata[];
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

HarmonicBond::HarmonicBond(std::vector<Real> K, std::vector<Real> r0)
    : K_(std::move(K)), r0_(std::move(r0)), state_(std::make_unique<State>()) {}

HarmonicBond::~HarmonicBond() = default;
HarmonicBond::HarmonicBond(HarmonicBond&&) noexcept = default;
HarmonicBond& HarmonicBond::operator=(HarmonicBond&&) noexcept = default;

void HarmonicBond::build(const BuildContext& ctx) {
    int nt = static_cast<int>(K_.size());
    state_->n_types = nt;
    cudaMalloc(&state_->d_K, nt * sizeof(Real));
    cudaMalloc(&state_->d_r0, nt * sizeof(Real));
    cudaMemcpy(state_->d_K, K_.data(), nt * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(state_->d_r0, r0_.data(), nt * sizeof(Real), cudaMemcpyHostToDevice);
}

void HarmonicBond::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    Int nb = f.n_bonds();
    if (nb == 0) return;

    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(nb));
    cuda::kernels::harmonic_bond_k<<<grid, block, block * sizeof(Real), stream>>>(
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
