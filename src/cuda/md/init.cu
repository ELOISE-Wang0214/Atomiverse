/// @file md/init.cu
/// @brief Initializer implementations: ZeroForces, ZeroEnergy, ZeroKinetic, InitVelocity.

#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/core/units.h"

#include <cuda_runtime.h>
#include "atomiverse/cuda/types.cuh"
#include <curand_kernel.h>

namespace ATV {

using cuda::real_t;
using cuda::work_t;
using cuda::prec_t;

namespace md {

// ── Zero initializers ────────────────────────────────────────────────

void ZeroForces::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    auto n = static_cast<size_t>(ctx.frame.n_atoms());
    cudaMemsetAsync(frame.atoms.fx, 0, n * sizeof(work_t), stream);
    cudaMemsetAsync(frame.atoms.fy, 0, n * sizeof(work_t), stream);
    cudaMemsetAsync(frame.atoms.fz, 0, n * sizeof(work_t), stream);
}

void ZeroEnergy::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    cudaMemsetAsync(frame.energy, 0, sizeof(work_t), stream);
}

void ZeroKinetic::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    cudaMemsetAsync(frame.kinetic, 0, sizeof(work_t), stream);
}

// ── InitVelocity ─────────────────────────────────────────────────────

namespace {

/// @brief Generate Gaussian-distributed velocities at target temperature.
///
/// v_i ~ N(0, sigma) where sigma = sqrt(kB*T / (kMvv2e * m_i)).
/// Uses cuRAND device API (per-thread state, Philox generator).
__global__ void generate_velocities(
    real_t* __restrict__ vx, real_t* __restrict__ vy, real_t* __restrict__ vz,
    const real_t* __restrict__ mass, work_t kbt_over_mvv2e,
    unsigned long long seed, const int* __restrict__ n_ptr) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= *n_ptr) return;

    curandState state;
    curand_init(seed, i, 0, &state);

    work_t sigma = sqrt(kbt_over_mvv2e / mass[i]);
    vx[i] = static_cast<work_t>(curand_normal_double(&state)) * sigma;
    vy[i] = static_cast<work_t>(curand_normal_double(&state)) * sigma;
    vz[i] = static_cast<work_t>(curand_normal_double(&state)) * sigma;
}

/// @brief Accumulate total momentum and total mass via atomicAdd.
///
/// Build-time only (atomicAdd precision acceptable for initialization).
__global__ void accumulate_momentum(
    const real_t* __restrict__ vx, const real_t* __restrict__ vy,
    const real_t* __restrict__ vz, const real_t* __restrict__ mass,
    real_t* __restrict__ px, real_t* __restrict__ py, real_t* __restrict__ pz,
    real_t* __restrict__ total_mass, const int* __restrict__ n_ptr) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= *n_ptr) return;

    atomicAdd(px, mass[i] * vx[i]);
    atomicAdd(py, mass[i] * vy[i]);
    atomicAdd(pz, mass[i] * vz[i]);
    atomicAdd(total_mass, mass[i]);
}

/// @brief Subtract center-of-mass velocity from all atoms.
__global__ void remove_com(
    real_t* __restrict__ vx, real_t* __restrict__ vy, real_t* __restrict__ vz,
    work_t vcm_x, work_t vcm_y, work_t vcm_z, const int* __restrict__ n_ptr) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= *n_ptr) return;

    vx[i] -= vcm_x;
    vy[i] -= vcm_y;
    vz[i] -= vcm_z;
}

} // anonymous namespace

void InitVelocity::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    int n = ctx.frame.n_atoms();

    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    // sigma^2 = kB * T / (kMvv2e * m)  =>  pass kB*T/kMvv2e as common factor
    work_t kbt_over_mvv2e = kBoltzmann * temperature_ / kMvv2e;

    // 1. Generate Gaussian velocities
    generate_velocities<<<grid, block, 0, stream>>>(
        frame.atoms.vx, frame.atoms.vy, frame.atoms.vz, frame.atoms.mass,
        kbt_over_mvv2e, seed_, frame.atoms.d_size);

    // 2. Compute total momentum + total mass
    Real* d_px;  Real* d_py;  Real* d_pz;  Real* d_M;
    cudaMalloc(&d_px, sizeof(work_t));
    cudaMalloc(&d_py, sizeof(work_t));
    cudaMalloc(&d_pz, sizeof(work_t));
    cudaMalloc(&d_M,  sizeof(work_t));
    cudaMemsetAsync(d_px, 0, sizeof(work_t), stream);
    cudaMemsetAsync(d_py, 0, sizeof(work_t), stream);
    cudaMemsetAsync(d_pz, 0, sizeof(work_t), stream);
    cudaMemsetAsync(d_M,  0, sizeof(work_t), stream);

    accumulate_momentum<<<grid, block, 0, stream>>>(
        frame.atoms.vx, frame.atoms.vy, frame.atoms.vz, frame.atoms.mass,
        d_px, d_py, d_pz, d_M, frame.atoms.d_size);

    // 3. D2H momentum + mass (build-time D2H allowed)
    work_t px = 0, py = 0, pz = 0, M = 0;
    cudaStreamSynchronize(stream);
    cudaMemcpy(&px, d_px, sizeof(work_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&py, d_py, sizeof(work_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pz, d_pz, sizeof(work_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&M,  d_M,  sizeof(work_t), cudaMemcpyDeviceToHost);

    cudaFree(d_px); cudaFree(d_py); cudaFree(d_pz); cudaFree(d_M);

    // 4. Remove center-of-mass velocity
    if (M > work_t(0)) {
        work_t vcm_x = px / M;
        work_t vcm_y = py / M;
        work_t vcm_z = pz / M;
        remove_com<<<grid, block, 0, stream>>>(
            frame.atoms.vx, frame.atoms.vy, frame.atoms.vz,
            vcm_x, vcm_y, vcm_z, frame.atoms.d_size);
        cudaStreamSynchronize(stream);
    }
}

void InitVelocity::launch(const LaunchContext& /*ctx*/) {
    // No-op: velocities initialized once in build().
}

} // namespace md
} // namespace ATV
