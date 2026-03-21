/// @file md/nose_hoover.cu
/// @brief NoseHoover -- single-chain Nosé-Hoover thermostat.

#include "atomiverse/cuda/md/nose_hoover.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/facade/stage.h"
#include "atomiverse/core/units.h"
#include <cuda_runtime.h>

namespace ATV {
namespace md {

struct NoseHoover::State {
    Real* d_xi = nullptr;        // thermostat variable [1/fs]
    Real* d_ke_temp = nullptr;   // temporary KE [eV]
    Real* d_scale = nullptr;     // velocity scale factor [1 element]
    unsigned int* d_block_count = nullptr; // last-block detection counter
    Real Q = 0;                  // thermostat mass
    Int ndof = 0;
};

NoseHoover::NoseHoover(Real T_target, Real tau, Real dt)
    : T_target_(T_target), tau_(tau), dt_(dt),
      state_(std::make_unique<State>()) {}

NoseHoover::~NoseHoover() = default;
NoseHoover::NoseHoover(NoseHoover&&) noexcept = default;
NoseHoover& NoseHoover::operator=(NoseHoover&&) noexcept = default;

std::vector<DefaultBinding> NoseHoover::default_bindings() const {
    return {{static_cast<int>(Phase::initial_integrate), 100},
            {static_cast<int>(Phase::final_integrate), 100}};
}

namespace {

/// @brief Fused KE reduction + thermostat update.
///
/// Reduces per-atom kinetic energy via shared-memory block reduction and
/// atomicAdd. The last block to finish (detected via atomicInc counter)
/// updates the thermostat variable xi and computes the velocity scale
/// factor, eliminating the separate <<<1,1>>> kernel launch.
__global__ void reduce_ke_and_update(
    const Real* __restrict__ vx, const Real* __restrict__ vy, const Real* __restrict__ vz,
    const Real* __restrict__ mass, Real* __restrict__ ke_out,
    Int n, Real kMvv2e_val,
    Real* __restrict__ xi, Real* __restrict__ scale,
    Real Q, Real ndof_kbT, Real half_dt,
    unsigned int* __restrict__ block_count)
{
    extern __shared__ Real sdata[];
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    Real val = 0;
    if (idx < n) {
        Real m = mass[idx];
        val = Real(0.5) * kMvv2e_val * m *
              (vx[idx] * vx[idx] + vy[idx] * vy[idx] + vz[idx] * vz[idx]);
    }
    sdata[threadIdx.x] = val;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sdata[threadIdx.x] += sdata[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) atomicAdd(ke_out, sdata[0]);

    // Last-block pattern: ensure atomicAdd is globally visible, then detect
    // the last block to finish. Thread 0 of the last block updates xi/scale.
    __threadfence();
    __shared__ bool is_last;
    if (threadIdx.x == 0) {
        unsigned int old = atomicInc(block_count, gridDim.x - 1);
        is_last = (old == gridDim.x - 1);
    }
    __syncthreads();

    if (is_last && threadIdx.x == 0) {
        Real ke_val = *ke_out;
        *xi += half_dt * (Real(2) * ke_val - ndof_kbT) / Q;
        *scale = exp(-half_dt * (*xi));
    }
}

/// @brief Scale all velocities by the precomputed factor.
__global__ void scale_vel(
    Real* __restrict__ vx, Real* __restrict__ vy, Real* __restrict__ vz,
    const Real* __restrict__ scale, Int n)
{
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Real s = *scale;
    vx[idx] *= s;
    vy[idx] *= s;
    vz[idx] *= s;
}

} // anonymous namespace

void NoseHoover::build(const BuildContext& ctx) {
    Int n = ctx.frame.n_atoms();
    state_->ndof = 3 * n - 3;
    state_->Q = static_cast<Real>(state_->ndof) * kBoltzmann * T_target_ * tau_ * tau_;

    cudaMalloc(&state_->d_xi, sizeof(Real));
    cudaMalloc(&state_->d_ke_temp, sizeof(Real));
    cudaMalloc(&state_->d_scale, sizeof(Real));
    cudaMalloc(&state_->d_block_count, sizeof(unsigned int));
    cudaMemset(state_->d_xi, 0, sizeof(Real));
    cudaMemset(state_->d_block_count, 0, sizeof(unsigned int));

    // Initialize scale to 1.0
    Real one = Real(1);
    cudaMemcpy(state_->d_scale, &one, sizeof(Real), cudaMemcpyHostToDevice);
}

void NoseHoover::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    Int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    Real ndof_kbT = static_cast<Real>(state_->ndof) * kBoltzmann * T_target_;
    Real half_dt = Real(0.5) * dt_;

    // Fused: reduce KE + update xi/scale (last block computes scale)
    cudaMemsetAsync(state_->d_ke_temp, 0, sizeof(Real), stream);
    reduce_ke_and_update<<<grid, block, block * sizeof(Real), stream>>>(
        f.vx, f.vy, f.vz, f.mass, state_->d_ke_temp,
        n, kMvv2e,
        state_->d_xi, state_->d_scale,
        state_->Q, ndof_kbT, half_dt,
        state_->d_block_count);

    // Scale velocities (reads scale written by the last block above)
    scale_vel<<<grid, block, 0, stream>>>(
        f.vx, f.vy, f.vz, state_->d_scale, n);
}

void NoseHoover::destroy(const DestroyContext& ctx) {
    cudaFree(state_->d_xi);
    cudaFree(state_->d_ke_temp);
    cudaFree(state_->d_scale);
    cudaFree(state_->d_block_count);
    *state_ = {};
}

} // namespace md
} // namespace ATV
