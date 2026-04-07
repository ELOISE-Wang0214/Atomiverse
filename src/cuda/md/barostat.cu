/// @file md/barostat.cu
/// @brief Barostat -- Berendsen isotropic pressure coupling.

#include "atomiverse/cuda/md/barostat.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/facade/stage.h"
#include <cuda_runtime.h>
#include "atomiverse/cuda/types.cuh"
#include <cmath>

namespace ATV {

using cuda::real_t;
using cuda::work_t;
using cuda::prec_t;

namespace md {

struct Barostat::State {
    Real* d_mu = nullptr;  // scale factor [device, 1 element]
    int n_atoms = 0;
};

Barostat::Barostat(work_t P_target, work_t tau_p, work_t dt)
    : P_target_(P_target), tau_p_(tau_p), dt_(dt),
      state_(std::make_unique<State>()) {}

Barostat::~Barostat() = default;
Barostat::Barostat(Barostat&&) noexcept = default;
Barostat& Barostat::operator=(Barostat&&) noexcept = default;

std::vector<DefaultBinding> Barostat::default_bindings() const {
    return {{static_cast<int>(Phase::final_integrate), 200}};
}

namespace {

// Compute μ from pressure and target, then scale box + positions
__global__ void berendsen_scale_k(
    real_t* __restrict__ x, real_t* __restrict__ y, real_t* __restrict__ z,
    real_t* __restrict__ box, const real_t* __restrict__ stress,
    work_t P_target, work_t beta_dt_tau, int n)
{
    // Isotropic pressure = -(σ_xx + σ_yy + σ_zz) / 3
    // Compute P and mu once per block via shared memory
    __shared__ work_t s_mu;
    if (threadIdx.x == 0) {
        work_t P = -(stress[0] + stress[1] + stress[2]) / work_t(3);
        s_mu = cbrt(work_t(1) - beta_dt_tau * (P_target - P));
    }
    __syncthreads();
    work_t mu = s_mu;

    // Thread 0 of block 0 rescales box
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        box[0] *= mu; box[4] *= mu; box[8] *= mu;
    }

    // All threads rescale positions
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    x[idx] *= mu;
    y[idx] *= mu;
    z[idx] *= mu;
}

} // anonymous namespace

void Barostat::build(const BuildContext& ctx) {
    state_->n_atoms = ctx.frame.n_atoms();
}

void Barostat::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    // β is the isothermal compressibility; for Berendsen we use β*dt/τ_p as a single parameter.
    // Typical β for solids ~ 1e-4 GPa^-1. We absorb it into the coupling.
    work_t beta_dt_tau = dt_ / tau_p_;

    berendsen_scale_k<<<grid, block, 0, stream>>>(
        f.atoms.x, f.atoms.y, f.atoms.z, f.box, f.stress,
        P_target_, beta_dt_tau, n);
}

void Barostat::destroy(const DestroyContext& ctx) {
    *state_ = {};
}

} // namespace md
} // namespace ATV
