/// @file md/init.cu
/// @brief Initializer implementations: ZeroForces, ZeroEnergy, ZeroKinetic, InitVelocity.

#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/core/units.h"

#include <cuda_runtime.h>
#include <curand_kernel.h>

namespace ATV {
namespace md {

// ── Zero initializers ────────────────────────────────────────────────

void ZeroForces::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    auto n = static_cast<size_t>(ctx.frame.n_atoms());
    cudaMemsetAsync(frame.atoms.fx, 0, n * sizeof(Real), stream);
    cudaMemsetAsync(frame.atoms.fy, 0, n * sizeof(Real), stream);
    cudaMemsetAsync(frame.atoms.fz, 0, n * sizeof(Real), stream);
}

void ZeroEnergy::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    cudaMemsetAsync(frame.energy, 0, sizeof(Real), stream);
}

void ZeroKinetic::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    cudaMemsetAsync(frame.kinetic, 0, sizeof(Real), stream);
}

// ── InitVelocity ─────────────────────────────────────────────────────

namespace {

/// @brief Generate Gaussian-distributed velocities at target temperature.
///
/// v_i ~ N(0, sigma) where sigma = sqrt(kB*T / (kMvv2e * m_i)).
/// Uses cuRAND device API (per-thread state, Philox generator).
__global__ void generate_velocities(
    Real* __restrict__ vx, Real* __restrict__ vy, Real* __restrict__ vz,
    const Real* __restrict__ mass, Real kbt_over_mvv2e,
    unsigned long long seed, const Int* __restrict__ n_ptr) {

    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= *n_ptr) return;

    curandState state;
    curand_init(seed, i, 0, &state);

    Real sigma = sqrt(kbt_over_mvv2e / mass[i]);
    vx[i] = static_cast<Real>(curand_normal_double(&state)) * sigma;
    vy[i] = static_cast<Real>(curand_normal_double(&state)) * sigma;
    vz[i] = static_cast<Real>(curand_normal_double(&state)) * sigma;
}

/// @brief Accumulate total momentum and total mass via atomicAdd.
///
/// Build-time only (atomicAdd precision acceptable for initialization).
__global__ void accumulate_momentum(
    const Real* __restrict__ vx, const Real* __restrict__ vy,
    const Real* __restrict__ vz, const Real* __restrict__ mass,
    Real* __restrict__ px, Real* __restrict__ py, Real* __restrict__ pz,
    Real* __restrict__ total_mass, const Int* __restrict__ n_ptr) {

    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= *n_ptr) return;

    atomicAdd(px, mass[i] * vx[i]);
    atomicAdd(py, mass[i] * vy[i]);
    atomicAdd(pz, mass[i] * vz[i]);
    atomicAdd(total_mass, mass[i]);
}

/// @brief Subtract center-of-mass velocity from all atoms.
__global__ void remove_com(
    Real* __restrict__ vx, Real* __restrict__ vy, Real* __restrict__ vz,
    Real vcm_x, Real vcm_y, Real vcm_z, const Int* __restrict__ n_ptr) {

    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= *n_ptr) return;

    vx[i] -= vcm_x;
    vy[i] -= vcm_y;
    vz[i] -= vcm_z;
}

} // anonymous namespace

void InitVelocity::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    Int n = ctx.frame.n_atoms();

    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    // sigma^2 = kB * T / (kMvv2e * m)  =>  pass kB*T/kMvv2e as common factor
    Real kbt_over_mvv2e = kBoltzmann * temperature_ / kMvv2e;

    // 1. Generate Gaussian velocities
    generate_velocities<<<grid, block, 0, stream>>>(
        frame.atoms.vx, frame.atoms.vy, frame.atoms.vz, frame.atoms.mass,
        kbt_over_mvv2e, seed_, frame.atoms.d_size);

    // 2. Compute total momentum + total mass
    Real* d_px;  Real* d_py;  Real* d_pz;  Real* d_M;
    cudaMalloc(&d_px, sizeof(Real));
    cudaMalloc(&d_py, sizeof(Real));
    cudaMalloc(&d_pz, sizeof(Real));
    cudaMalloc(&d_M,  sizeof(Real));
    cudaMemsetAsync(d_px, 0, sizeof(Real), stream);
    cudaMemsetAsync(d_py, 0, sizeof(Real), stream);
    cudaMemsetAsync(d_pz, 0, sizeof(Real), stream);
    cudaMemsetAsync(d_M,  0, sizeof(Real), stream);

    accumulate_momentum<<<grid, block, 0, stream>>>(
        frame.atoms.vx, frame.atoms.vy, frame.atoms.vz, frame.atoms.mass,
        d_px, d_py, d_pz, d_M, frame.atoms.d_size);

    // 3. D2H momentum + mass (build-time D2H allowed)
    Real px = 0, py = 0, pz = 0, M = 0;
    cudaStreamSynchronize(stream);
    cudaMemcpy(&px, d_px, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&py, d_py, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pz, d_pz, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&M,  d_M,  sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(d_px); cudaFree(d_py); cudaFree(d_pz); cudaFree(d_M);

    // 4. Remove center-of-mass velocity
    if (M > Real(0)) {
        Real vcm_x = px / M;
        Real vcm_y = py / M;
        Real vcm_z = pz / M;
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
