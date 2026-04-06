/// @file md/velocity_verlet.cu

#include "atomiverse/cuda/md/velocity_verlet.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/facade/stage.h"
#include "atomiverse/core/units.h"
#include <cuda_runtime.h>

namespace ATV {
namespace md {

namespace {

__global__ void initial_integrate(
    Real* __restrict__ x, Real* __restrict__ y, Real* __restrict__ z,
    Real* __restrict__ vx, Real* __restrict__ vy, Real* __restrict__ vz,
    const Real* __restrict__ fx, const Real* __restrict__ fy,
    const Real* __restrict__ fz, const Real* __restrict__ mass,
    const Real* __restrict__ box, Real dt, Int n) {
    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    Real inv_m = Real(1) / mass[i];
    Real hdt = Real(0.5) * dt;
    Real ax = fx[i] * inv_m * kFtm2v;
    Real ay = fy[i] * inv_m * kFtm2v;
    Real az = fz[i] * inv_m * kFtm2v;
    vx[i] += hdt * ax;
    vy[i] += hdt * ay;
    vz[i] += hdt * az;
    x[i] += dt * vx[i];
    y[i] += dt * vy[i];
    z[i] += dt * vz[i];
    // PBC wrapping (orthorhombic: diagonal elements)
    Real Lx = box[0], Ly = box[4], Lz = box[8];
    Real inv_Lx = Real(1) / Lx, inv_Ly = Real(1) / Ly, inv_Lz = Real(1) / Lz;
    x[i] -= Lx * floor(x[i] * inv_Lx);
    y[i] -= Ly * floor(y[i] * inv_Ly);
    z[i] -= Lz * floor(z[i] * inv_Lz);
}

__global__ void final_integrate(
    Real* __restrict__ vx, Real* __restrict__ vy, Real* __restrict__ vz,
    const Real* __restrict__ fx, const Real* __restrict__ fy,
    const Real* __restrict__ fz, const Real* __restrict__ mass,
    Real dt, Int n) {
    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    Real inv_m = Real(1) / mass[i];
    Real hdt = Real(0.5) * dt;
    vx[i] += hdt * fx[i] * inv_m * kFtm2v;
    vy[i] += hdt * fy[i] * inv_m * kFtm2v;
    vz[i] += hdt * fz[i] * inv_m * kFtm2v;
}

} // anonymous namespace

void VelocityVerlet::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    Int n = ctx.frame.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));
    if (ctx.stage == static_cast<int>(Phase::initial_integrate)) {
        initial_integrate<<<grid, block, 0, stream>>>(
            frame.atoms.x, frame.atoms.y, frame.atoms.z, frame.atoms.vx, frame.atoms.vy, frame.atoms.vz,
            frame.atoms.fx, frame.atoms.fy, frame.atoms.fz, frame.atoms.mass,
            frame.box, dt_, n);
    } else {
        final_integrate<<<grid, block, 0, stream>>>(
            frame.atoms.vx, frame.atoms.vy, frame.atoms.vz, frame.atoms.fx, frame.atoms.fy, frame.atoms.fz,
            frame.atoms.mass, dt_, n);
    }
}

} // namespace md
} // namespace ATV
