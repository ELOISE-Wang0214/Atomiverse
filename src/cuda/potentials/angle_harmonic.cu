/// @file cuda/potentials/angle_harmonic.cu
/// @brief CUDA AngleHarmonic — V = K(θ − θ₀)² (kernel + class).

#include "atomiverse/core/constants.h"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/potentials/angle_harmonic.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cmath>
#include <cuda_runtime.h>

namespace ATV {
namespace cuda {

struct AngleHarmonic::State {
    work_t* d_K = nullptr;
    work_t* d_theta0 = nullptr; // in radians
    int n_types = 0;
};

AngleHarmonic::AngleHarmonic(std::vector<double> K, std::vector<double> theta0)
    : K_(std::move(K)), theta0_(std::move(theta0)),
      state_(std::make_unique<State>()) {}

AngleHarmonic::~AngleHarmonic() = default;
AngleHarmonic::AngleHarmonic(AngleHarmonic&&) noexcept = default;
AngleHarmonic& AngleHarmonic::operator=(AngleHarmonic&&) noexcept = default;

namespace {

__global__ void
harmonic_angle_k(const real_t* __restrict__ x, const real_t* __restrict__ y,
                 const real_t* __restrict__ z, real_t* __restrict__ fx,
                 real_t* __restrict__ fy, real_t* __restrict__ fz,
                 prec_t* __restrict__ energy, const int* __restrict__ ai,
                 const int* __restrict__ aj, const int* __restrict__ ak,
                 const int* __restrict__ at, const work_t* __restrict__ K,
                 const work_t* __restrict__ theta0,
                 const prec_t* __restrict__ box, int n_angles) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_angles)
        return;

    int i = ai[idx], j = aj[idx], k = ak[idx];
    int type = at[idx];

    // Vectors from center atom j
    work_t rij_x = x[i] - x[j], rij_y = y[i] - y[j], rij_z = z[i] - z[j];
    work_t rkj_x = x[k] - x[j], rkj_y = y[k] - y[j], rkj_z = z[k] - z[j];

    // Minimum image
    work_t bx = box[0], by = box[4], bz = box[8];
    rij_x -= bx * round(rij_x / bx);
    rij_y -= by * round(rij_y / by);
    rij_z -= bz * round(rij_z / bz);
    rkj_x -= bx * round(rkj_x / bx);
    rkj_y -= by * round(rkj_y / by);
    rkj_z -= bz * round(rkj_z / bz);

    work_t rij2 = rij_x * rij_x + rij_y * rij_y + rij_z * rij_z;
    work_t rkj2 = rkj_x * rkj_x + rkj_y * rkj_y + rkj_z * rkj_z;
    work_t rij_inv = rsqrt(rij2);
    work_t rkj_inv = rsqrt(rkj2);

    work_t cos_theta =
        (rij_x * rkj_x + rij_y * rkj_y + rij_z * rkj_z) * rij_inv * rkj_inv;
    cos_theta = fmin(fmax(cos_theta, work_t(-1)), work_t(1));
    work_t theta = acos(cos_theta);
    work_t dtheta = theta - theta0[type];

    // V = K(θ - θ0)^2
    work_t v = K[type] * dtheta * dtheta;

    // dV/dθ = 2K(θ - θ0)
    work_t sin_theta = sqrt(work_t(1) - cos_theta * cos_theta);
    constexpr work_t kSinEpsilon = work_t(1e-6);
    if (sin_theta < kSinEpsilon)
        sin_theta = kSinEpsilon;
    work_t prefactor = -work_t(2) * K[type] * dtheta / sin_theta;

    // Force on atom i: prefactor * (rkj_hat / |rij| - cos_theta * rij_hat / |rij|)
    work_t rij_inv2 = rij_inv * rij_inv;
    work_t rkj_inv2 = rkj_inv * rkj_inv;
    work_t cross_inv = rij_inv * rkj_inv;

    work_t fi_x =
        prefactor * (rkj_x * cross_inv - cos_theta * rij_x * rij_inv2);
    work_t fi_y =
        prefactor * (rkj_y * cross_inv - cos_theta * rij_y * rij_inv2);
    work_t fi_z =
        prefactor * (rkj_z * cross_inv - cos_theta * rij_z * rij_inv2);

    work_t fk_x =
        prefactor * (rij_x * cross_inv - cos_theta * rkj_x * rkj_inv2);
    work_t fk_y =
        prefactor * (rij_y * cross_inv - cos_theta * rkj_y * rkj_inv2);
    work_t fk_z =
        prefactor * (rij_z * cross_inv - cos_theta * rkj_z * rkj_inv2);

    atomicAdd(&fx[i], static_cast<real_t>(fi_x));
    atomicAdd(&fy[i], static_cast<real_t>(fi_y));
    atomicAdd(&fz[i], static_cast<real_t>(fi_z));
    atomicAdd(&fx[k], static_cast<real_t>(fk_x));
    atomicAdd(&fy[k], static_cast<real_t>(fk_y));
    atomicAdd(&fz[k], static_cast<real_t>(fk_z));
    atomicAdd(&fx[j], static_cast<real_t>(-(fi_x + fk_x)));
    atomicAdd(&fy[j], static_cast<real_t>(-(fi_y + fk_y)));
    atomicAdd(&fz[j], static_cast<real_t>(-(fi_z + fk_z)));

    // Energy reduction
    extern __shared__ work_t sdata[];
    sdata[threadIdx.x] = v;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s)
            sdata[threadIdx.x] += sdata[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0)
        atomicAdd(energy, static_cast<prec_t>(sdata[0]));
}

} // anonymous namespace

void AngleHarmonic::build(const BuildContext& /*ctx*/) {
    int nt = static_cast<int>(K_.size());
    state_->n_types = nt;

    // Convert degrees to radians
    std::vector<work_t> theta0_rad(nt);
    for (int i = 0; i < nt; ++i)
        theta0_rad[i] = theta0_[i] * work_t(kDegToRad);

    cudaMalloc(&state_->d_K, nt * sizeof(work_t));
    cudaMalloc(&state_->d_theta0, nt * sizeof(work_t));
    cudaMemcpy(state_->d_K, K_.data(), nt * sizeof(work_t),
               cudaMemcpyHostToDevice);
    cudaMemcpy(state_->d_theta0, theta0_rad.data(), nt * sizeof(work_t),
               cudaMemcpyHostToDevice);
}

void AngleHarmonic::compute(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int na = f.n_angles();
    if (na == 0)
        return;

    int block = kDefaultBlock;
    int grid = grid_size(static_cast<int>(na));
    harmonic_angle_k<<<grid, block, block * sizeof(work_t), stream>>>(
        f.atoms.x, f.atoms.y, f.atoms.z, f.atoms.fx, f.atoms.fy, f.atoms.fz,
        f.energy, f.angles.atomi, f.angles.atomj, f.angles.atomk, f.angles.type,
        state_->d_K, state_->d_theta0, f.box, na);
}

void AngleHarmonic::destroy(const DestroyContext& /*ctx*/) {
    cudaFree(state_->d_K);
    cudaFree(state_->d_theta0);
    *state_ = {};
}

} // namespace cuda
} // namespace ATV
