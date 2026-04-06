/// @file bonded/harmonic_angle.cu
/// @brief HarmonicAngle implementation -- V(θ) = K(θ - θ0)^2.

#include "atomiverse/cuda/potentials/harmonic_angle.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/core/constants.h"
#include <cuda_runtime.h>
#include <cmath>

namespace ATV {

struct HarmonicAngle::State {
    Real* d_K = nullptr;
    Real* d_theta0 = nullptr;  // in radians
    int n_types = 0;
};

HarmonicAngle::HarmonicAngle(std::vector<Real> K, std::vector<Real> theta0)
    : K_(std::move(K)), theta0_(std::move(theta0)),
      state_(std::make_unique<State>()) {}

HarmonicAngle::~HarmonicAngle() = default;
HarmonicAngle::HarmonicAngle(HarmonicAngle&&) noexcept = default;
HarmonicAngle& HarmonicAngle::operator=(HarmonicAngle&&) noexcept = default;

namespace {

__global__ void harmonic_angle_k(
    const Real* __restrict__ x, const Real* __restrict__ y, const Real* __restrict__ z,
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ ai, const Int* __restrict__ aj, const Int* __restrict__ ak,
    const Int* __restrict__ at,
    const Real* __restrict__ K, const Real* __restrict__ theta0,
    const Real* __restrict__ box, Int n_angles)
{
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_angles) return;

    Int i = ai[idx], j = aj[idx], k = ak[idx];
    Int type = at[idx];

    // Vectors from center atom j
    Real rij_x = x[i] - x[j], rij_y = y[i] - y[j], rij_z = z[i] - z[j];
    Real rkj_x = x[k] - x[j], rkj_y = y[k] - y[j], rkj_z = z[k] - z[j];

    // Minimum image
    Real bx = box[0], by = box[4], bz = box[8];
    rij_x -= bx * round(rij_x / bx);
    rij_y -= by * round(rij_y / by);
    rij_z -= bz * round(rij_z / bz);
    rkj_x -= bx * round(rkj_x / bx);
    rkj_y -= by * round(rkj_y / by);
    rkj_z -= bz * round(rkj_z / bz);

    Real rij2 = rij_x * rij_x + rij_y * rij_y + rij_z * rij_z;
    Real rkj2 = rkj_x * rkj_x + rkj_y * rkj_y + rkj_z * rkj_z;
    Real rij_inv = rsqrt(rij2);
    Real rkj_inv = rsqrt(rkj2);

    Real cos_theta = (rij_x * rkj_x + rij_y * rkj_y + rij_z * rkj_z) * rij_inv * rkj_inv;
    cos_theta = fmin(fmax(cos_theta, Real(-1)), Real(1));
    Real theta = acos(cos_theta);
    Real dtheta = theta - theta0[type];

    // V = K(θ - θ0)^2
    Real v = K[type] * dtheta * dtheta;

    // dV/dθ = 2K(θ - θ0)
    Real sin_theta = sqrt(Real(1) - cos_theta * cos_theta);
    constexpr Real kSinEpsilon = Real(1e-6);
    if (sin_theta < kSinEpsilon) sin_theta = kSinEpsilon;
    Real prefactor = -Real(2) * K[type] * dtheta / sin_theta;

    // Force on atom i: prefactor * (rkj_hat / |rij| - cos_theta * rij_hat / |rij|)
    Real rij_inv2 = rij_inv * rij_inv;
    Real rkj_inv2 = rkj_inv * rkj_inv;
    Real cross_inv = rij_inv * rkj_inv;

    Real fi_x = prefactor * (rkj_x * cross_inv - cos_theta * rij_x * rij_inv2);
    Real fi_y = prefactor * (rkj_y * cross_inv - cos_theta * rij_y * rij_inv2);
    Real fi_z = prefactor * (rkj_z * cross_inv - cos_theta * rij_z * rij_inv2);

    Real fk_x = prefactor * (rij_x * cross_inv - cos_theta * rkj_x * rkj_inv2);
    Real fk_y = prefactor * (rij_y * cross_inv - cos_theta * rkj_y * rkj_inv2);
    Real fk_z = prefactor * (rij_z * cross_inv - cos_theta * rkj_z * rkj_inv2);

    atomicAdd(&fx[i], fi_x); atomicAdd(&fy[i], fi_y); atomicAdd(&fz[i], fi_z);
    atomicAdd(&fx[k], fk_x); atomicAdd(&fy[k], fk_y); atomicAdd(&fz[k], fk_z);
    atomicAdd(&fx[j], -(fi_x + fk_x));
    atomicAdd(&fy[j], -(fi_y + fk_y));
    atomicAdd(&fz[j], -(fi_z + fk_z));

    // Energy reduction
    extern __shared__ Real sdata[];
    sdata[threadIdx.x] = v;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sdata[threadIdx.x] += sdata[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) atomicAdd(energy, sdata[0]);
}

} // anonymous namespace

void HarmonicAngle::build(const BuildContext& ctx) {
    int nt = static_cast<int>(K_.size());
    state_->n_types = nt;

    // Convert degrees to radians
    std::vector<Real> theta0_rad(nt);
    for (int i = 0; i < nt; ++i)
        theta0_rad[i] = theta0_[i] * Real(kDegToRad);

    cudaMalloc(&state_->d_K, nt * sizeof(Real));
    cudaMalloc(&state_->d_theta0, nt * sizeof(Real));
    cudaMemcpy(state_->d_K, K_.data(), nt * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(state_->d_theta0, theta0_rad.data(), nt * sizeof(Real), cudaMemcpyHostToDevice);
}

void HarmonicAngle::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    Int na = f.n_angles();
    if (na == 0) return;

    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(na));
    harmonic_angle_k<<<grid, block, block * sizeof(Real), stream>>>(
        f.atoms.x, f.atoms.y, f.atoms.z, f.atoms.fx, f.atoms.fy, f.atoms.fz, f.energy,
        f.angles.atomi, f.angles.atomj, f.angles.atomk, f.angles.type,
        state_->d_K, state_->d_theta0, f.box, na);
}

void HarmonicAngle::destroy(const DestroyContext& ctx) {
    cudaFree(state_->d_K);
    cudaFree(state_->d_theta0);
    *state_ = {};
}

} // namespace ATV
