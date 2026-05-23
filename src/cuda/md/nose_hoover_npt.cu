/// @file md/nose_hoover_npt.cu
/// @brief CUDA MTK NPT integrator implementation. Mirrors the CPU
///        sibling cpu::md::NoseHooverNPT line for line — see that file
///        for the canonical algorithm documentation block.

#include "atomiverse/core/units.h"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/md/nose_hoover_npt.cuh"
#include "atomiverse/cuda/types.cuh"
#include "atomiverse/facade/stage.h"
#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

namespace cuda {
namespace md {

namespace {

/// @brief Smooth `(e^x - 1)/x` with a Taylor branch around 0 to avoid
///        catastrophic cancellation. Identical to the CPU `expm1c`.
__host__ inline double expm1c(double x) {
    if (std::fabs(x) < 1e-4) {
        // 1 + x/2 + x²/6 — accurate to O(x³)
        return 1.0 + x * (0.5 + x * (1.0 / 6.0));
    }
    return std::expm1(x) / x;
}

constexpr double kEpsDtMaxFraction = 0.005;

// ── Per-atom kernels ────────────────────────────────────────────────

__global__ void scale_vel_uniform_k(real_t* __restrict__ vx,
                                    real_t* __restrict__ vy,
                                    real_t* __restrict__ vz, work_t scale,
                                    int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    vx[i] = static_cast<real_t>(static_cast<work_t>(vx[i]) * scale);
    vy[i] = static_cast<real_t>(static_cast<work_t>(vy[i]) * scale);
    vz[i] = static_cast<real_t>(static_cast<work_t>(vz[i]) * scale);
}

__global__ void scale_vel_per_axis_k(real_t* __restrict__ vx,
                                     real_t* __restrict__ vy,
                                     real_t* __restrict__ vz, work_t sx,
                                     work_t sy, work_t sz, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    vx[i] = static_cast<real_t>(static_cast<work_t>(vx[i]) * sx);
    vy[i] = static_cast<real_t>(static_cast<work_t>(vy[i]) * sy);
    vz[i] = static_cast<real_t>(static_cast<work_t>(vz[i]) * sz);
}

__global__ void kick_v_k(real_t* __restrict__ vx, real_t* __restrict__ vy,
                         real_t* __restrict__ vz, const real_t* __restrict__ fx,
                         const real_t* __restrict__ fy,
                         const real_t* __restrict__ fz,
                         const real_t* __restrict__ mass, work_t half_dt,
                         int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    work_t dtfm = half_dt * kFtm2v / static_cast<work_t>(mass[i]);
    vx[i] = static_cast<real_t>(static_cast<work_t>(vx[i])
                                + dtfm * static_cast<work_t>(fx[i]));
    vy[i] = static_cast<real_t>(static_cast<work_t>(vy[i])
                                + dtfm * static_cast<work_t>(fy[i]));
    vz[i] = static_cast<real_t>(static_cast<work_t>(vz[i])
                                + dtfm * static_cast<work_t>(fz[i]));
}

__global__ void propagate_x_k(real_t* __restrict__ x, real_t* __restrict__ y,
                              real_t* __restrict__ z,
                              const real_t* __restrict__ vx,
                              const real_t* __restrict__ vy,
                              const real_t* __restrict__ vz, work_t ax,
                              work_t ay, work_t az, work_t dx, work_t dy,
                              work_t dz, work_t dt, prec_t* __restrict__ box,
                              work_t lx_new, work_t ly_new, work_t lz_new,
                              int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        // Closed-form per-axis MTK position propagator (matches CPU
        // propagate_x_and_box).
        x[i] = static_cast<real_t>(static_cast<work_t>(x[i]) * ax
                                   + static_cast<work_t>(vx[i]) * dt * dx);
        y[i] = static_cast<real_t>(static_cast<work_t>(y[i]) * ay
                                   + static_cast<work_t>(vy[i]) * dt * dy);
        z[i] = static_cast<real_t>(static_cast<work_t>(z[i]) * az
                                   + static_cast<work_t>(vz[i]) * dt * dz);
        // Orthorhombic PBC wrap into the freshly-scaled box.
        work_t inv_lx = work_t(1) / lx_new;
        work_t inv_ly = work_t(1) / ly_new;
        work_t inv_lz = work_t(1) / lz_new;
        x[i] -= static_cast<real_t>(
            lx_new * floor(static_cast<work_t>(x[i]) * inv_lx));
        y[i] -= static_cast<real_t>(
            ly_new * floor(static_cast<work_t>(y[i]) * inv_ly));
        z[i] -= static_cast<real_t>(
            lz_new * floor(static_cast<work_t>(z[i]) * inv_lz));
    }
    // Thread 0 of block 0 writes the new box diagonal.
    if (i == 0) {
        box[0] = static_cast<prec_t>(lx_new);
        box[4] = static_cast<prec_t>(ly_new);
        box[8] = static_cast<prec_t>(lz_new);
    }
}

// ── KE reductions ───────────────────────────────────────────────────
//
// Kernels write to a 4-double device buffer:
//   d_ke_buf[0] = Σ m·(vx² + vy² + vz²)   (global ke pre-factor)
//   d_ke_buf[1] = Σ m·vx²                  (per-axis x)
//   d_ke_buf[2] = Σ m·vy²
//   d_ke_buf[3] = Σ m·vz²
// Caller multiplies by 0.5·kMvv2e to get eV.

__global__ void reduce_ke_axes_k(const real_t* __restrict__ vx,
                                 const real_t* __restrict__ vy,
                                 const real_t* __restrict__ vz,
                                 const real_t* __restrict__ mass,
                                 double* __restrict__ buf, int n) {
    extern __shared__ double s_buf[]; // 4 * blockDim.x
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + tid;
    double mx = 0.0, my = 0.0, mz = 0.0;
    if (i < n) {
        double m = static_cast<double>(mass[i]);
        double vxi = static_cast<double>(vx[i]);
        double vyi = static_cast<double>(vy[i]);
        double vzi = static_cast<double>(vz[i]);
        mx = m * vxi * vxi;
        my = m * vyi * vyi;
        mz = m * vzi * vzi;
    }
    int blk = blockDim.x;
    s_buf[tid] = mx;
    s_buf[blk + tid] = my;
    s_buf[2 * blk + tid] = mz;
    __syncthreads();
    for (int off = blk / 2; off > 0; off >>= 1) {
        if (tid < off) {
            s_buf[tid] += s_buf[tid + off];
            s_buf[blk + tid] += s_buf[blk + tid + off];
            s_buf[2 * blk + tid] += s_buf[2 * blk + tid + off];
        }
        __syncthreads();
    }
    if (tid == 0) {
        double bx = s_buf[0], by = s_buf[blk], bz = s_buf[2 * blk];
        atomicAdd(&buf[1], bx);
        atomicAdd(&buf[2], by);
        atomicAdd(&buf[3], bz);
        atomicAdd(&buf[0], bx + by + bz);
    }
}

} // anonymous namespace

// ── ctor + lifecycle ────────────────────────────────────────────────

NoseHooverNPT::NoseHooverNPT(Config cfg)
    : T_target_(cfg.T_target),
      T_stop_(cfg.T_stop < 0 ? cfg.T_target : cfg.T_stop), tau_T_(cfg.tau_T),
      P_target_(cfg.P_target), tau_P_(cfg.tau_P), dt_(cfg.dt),
      drag_(cfg.drag), mtchain_(cfg.mtchain), mpchain_(cfg.mpchain),
      couple_axis_(cfg.couple_axis), T_live_(cfg.T_target) {
    pdim_ = 0;
    for (int a = 0; a < 3; ++a)
        if (couple_axis_[a])
            ++pdim_;
}

void NoseHooverNPT::build(const BuildContext& ctx) {
    n_atoms_ = ctx.frame.n_atoms();
    n_dof_ = 3 * n_atoms_ - 3;
    kT_ = kBoltzmann * T_target_;

    Q_t_.assign(static_cast<std::size_t>(mtchain_), 0.0);
    Q_b_.assign(static_cast<std::size_t>(mpchain_), 0.0);
    eta_t_.assign(static_cast<std::size_t>(mtchain_), 0.0);
    eta_t_dot_.assign(static_cast<std::size_t>(mtchain_), 0.0);
    eta_b_.assign(static_cast<std::size_t>(mpchain_), 0.0);
    eta_b_dot_.assign(static_cast<std::size_t>(mpchain_), 0.0);

    const double tau_T2 = tau_T_ * tau_T_;
    const double tau_P2 = tau_P_ * tau_P_;
    Q_t_[0] = static_cast<double>(n_dof_) * kT_ * tau_T2;
    for (int k = 1; k < mtchain_; ++k)
        Q_t_[k] = kT_ * tau_T2;
    Q_b_[0] = kT_ * tau_P2;
    for (int k = 1; k < mpchain_; ++k)
        Q_b_[k] = kT_ * tau_P2;
    W_ = (static_cast<double>(n_dof_) + 3.0) * kT_ * tau_P2;

    tdrag_factor_ = 1.0 - drag_ * dt_ / tau_T_;
    pdrag_factor_ = 1.0 - drag_ * dt_ / tau_P_;

    // Initial KE at equipartition.
    KE_global_ = 0.5 * static_cast<double>(n_dof_) * kT_;

    // Cache box diagonal (orthorhombic only — same restriction as
    // cuda::md::Barostat / VelocityVerlet wrap).
    prec_t host_box[9];
    cudaMemcpy(host_box, ctx.frame.box, 9 * sizeof(prec_t),
               cudaMemcpyDeviceToHost);
    box_lx_ = static_cast<double>(host_box[0]);
    box_ly_ = static_cast<double>(host_box[4]);
    box_lz_ = static_cast<double>(host_box[8]);

    // Device scratch: 4 doubles for KE reductions.
    double* d_ke = nullptr;
    cudaMalloc(&d_ke, 4 * sizeof(double));
    d_ke_buf_ = d_ke;
}

void NoseHooverNPT::destroy(const DestroyContext& /*ctx*/) {
    if (d_ke_buf_) {
        cudaFree(d_ke_buf_);
        d_ke_buf_ = nullptr;
    }
}

// ── KE helpers ──────────────────────────────────────────────────────

double NoseHooverNPT::compute_global_KE_(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();
    auto* d_ke = static_cast<double*>(d_ke_buf_);
    cudaMemsetAsync(d_ke, 0, 4 * sizeof(double), stream);
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(n);
    size_t smem = 3 * block * sizeof(double);
    reduce_ke_axes_k<<<grid, block, smem, stream>>>(
        f.atoms.vx, f.atoms.vy, f.atoms.vz, f.atoms.mass, d_ke, n);
    cudaStreamSynchronize(stream);
    double h[4] = {0, 0, 0, 0};
    cudaMemcpy(h, d_ke, 4 * sizeof(double), cudaMemcpyDeviceToHost);
    return 0.5 * kMvv2e * h[0];
}

std::array<double, 3>
NoseHooverNPT::compute_axis_KE_(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();
    auto* d_ke = static_cast<double*>(d_ke_buf_);
    cudaMemsetAsync(d_ke, 0, 4 * sizeof(double), stream);
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(n);
    size_t smem = 3 * block * sizeof(double);
    reduce_ke_axes_k<<<grid, block, smem, stream>>>(
        f.atoms.vx, f.atoms.vy, f.atoms.vz, f.atoms.mass, d_ke, n);
    cudaStreamSynchronize(stream);
    double h[4] = {0, 0, 0, 0};
    cudaMemcpy(h, d_ke, 4 * sizeof(double), cudaMemcpyDeviceToHost);
    return {0.5 * kMvv2e * h[1], 0.5 * kMvv2e * h[2], 0.5 * kMvv2e * h[3]};
}

// ── Temperature ramp ────────────────────────────────────────────────

void NoseHooverNPT::compute_temp_target_(int step) {
    double delta = 0.0;
    if (nlast_ > nfirst_) {
        delta = static_cast<double>(step - nfirst_)
                / static_cast<double>(nlast_ - nfirst_);
        if (delta < 0.0)
            delta = 0.0;
        if (delta > 1.0)
            delta = 1.0;
    }
    T_live_ = T_target_ + delta * (T_stop_ - T_target_);
    kT_ = kBoltzmann * T_live_;
}

// ── NHC chain helpers ───────────────────────────────────────────────

void NoseHooverNPT::nhc_atoms_half_(const LaunchContext& ctx, double half_dt) {
    const double quarter_dt = 0.5 * half_dt;
    const double eighth_dt = 0.5 * quarter_dt;
    const int M = mtchain_;

    auto G_t = [&](int k) -> double {
        if (k == 0)
            return (2.0 * KE_global_ - n_dof_ * kT_) / Q_t_[0];
        return (Q_t_[k - 1] * eta_t_dot_[k - 1] * eta_t_dot_[k - 1] - kT_)
               / Q_t_[k];
    };

    // Backward pass with per-link drag.
    for (int k = M - 1; k >= 0; --k) {
        double expfac =
            (k < M - 1) ? std::exp(-eta_t_dot_[k + 1] * eighth_dt) : 1.0;
        eta_t_dot_[k] *= expfac;
        eta_t_dot_[k] += G_t(k) * quarter_dt;
        eta_t_dot_[k] *= tdrag_factor_;
        eta_t_dot_[k] *= expfac;
    }

    // Scale atom velocities by the innermost thermostat momentum (GPU kernel).
    double scale = std::exp(-eta_t_dot_[0] * half_dt);
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(n);
    scale_vel_uniform_k<<<grid, block, 0, stream>>>(
        f.atoms.vx, f.atoms.vy, f.atoms.vz, static_cast<work_t>(scale), n);
    KE_global_ *= scale * scale;

    for (int k = 0; k < M; ++k)
        eta_t_[k] += eta_t_dot_[k] * half_dt;

    // Forward pass, no drag.
    for (int k = 0; k < M; ++k) {
        double expfac =
            (k < M - 1) ? std::exp(-eta_t_dot_[k + 1] * eighth_dt) : 1.0;
        eta_t_dot_[k] *= expfac;
        eta_t_dot_[k] += G_t(k) * quarter_dt;
        eta_t_dot_[k] *= expfac;
    }
}

void NoseHooverNPT::nhc_baro_half_(double half_dt) {
    const double quarter_dt = 0.5 * half_dt;
    const double eighth_dt = 0.5 * quarter_dt;
    const int M = mpchain_;
    const double W_axis = W_ / 3.0;

    auto baro_2ke = [&]() {
        double s = 0.0;
        for (int a = 0; a < 3; ++a)
            if (couple_axis_[a])
                s += W_axis * eps_dot_[a] * eps_dot_[a];
        return s;
    };
    auto G_b = [&](int k) -> double {
        if (k == 0)
            return (baro_2ke() - kT_) / Q_b_[0];
        return (Q_b_[k - 1] * eta_b_dot_[k - 1] * eta_b_dot_[k - 1] - kT_)
               / Q_b_[k];
    };

    // Backward pass with per-link drag.
    for (int k = M - 1; k >= 0; --k) {
        double expfac =
            (k < M - 1) ? std::exp(-eta_b_dot_[k + 1] * eighth_dt) : 1.0;
        eta_b_dot_[k] *= expfac;
        eta_b_dot_[k] += G_b(k) * quarter_dt;
        eta_b_dot_[k] *= pdrag_factor_;
        eta_b_dot_[k] *= expfac;
    }

    double baro_scale = std::exp(-eta_b_dot_[0] * half_dt);
    for (int a = 0; a < 3; ++a)
        if (couple_axis_[a])
            eps_dot_[a] *= baro_scale;

    for (int k = 0; k < M; ++k)
        eta_b_[k] += eta_b_dot_[k] * half_dt;

    // Forward pass, no drag.
    for (int k = 0; k < M; ++k) {
        double expfac =
            (k < M - 1) ? std::exp(-eta_b_dot_[k + 1] * eighth_dt) : 1.0;
        eta_b_dot_[k] *= expfac;
        eta_b_dot_[k] += G_b(k) * quarter_dt;
        eta_b_dot_[k] *= expfac;
    }
}

// ── Pressure update ─────────────────────────────────────────────────

void NoseHooverNPT::update_eps_dot_(const LaunchContext& ctx, double half_dt) {
    const auto& f = ctx.frame;
    double V = box_lx_ * box_ly_ * box_lz_;
    if (V <= 0.0)
        return;

    // Per-axis kinetic energy (one reduction kernel + D2H).
    auto ke_axis = compute_axis_KE_(ctx);

    // D2H the pair-virial diagonal stress[0..2].
    double stress_diag[3] = {0.0, 0.0, 0.0};
    {
        double tmp[6];
        cudaMemcpy(tmp, &f.scalars->stress[0], 6 * sizeof(double),
                   cudaMemcpyDeviceToHost);
        stress_diag[0] = tmp[0];
        stress_diag[1] = tmp[1];
        stress_diag[2] = tmp[2];
    }

    double diag[3] = {stress_diag[0] + 2.0 * ke_axis[0],
                      stress_diag[1] + 2.0 * ke_axis[1],
                      stress_diag[2] + 2.0 * ke_axis[2]};
    if (pdim_ == 3) {
        double ave = (1.0 / 3.0) * (diag[0] + diag[1] + diag[2]);
        diag[0] = diag[1] = diag[2] = ave;
    }

    double mtk_active_ke = 0.0;
    for (int a = 0; a < 3; ++a)
        if (couple_axis_[a])
            mtk_active_ke += ke_axis[a];
    const int N_atoms = (n_dof_ + 3) / 3;
    const double mtk_term =
        (pdim_ > 0)
            ? 6.0 * mtk_active_ke / (static_cast<double>(pdim_) * N_atoms)
            : 0.0;

    for (int a = 0; a < 3; ++a) {
        if (!couple_axis_[a])
            continue;
        double P_axis = diag[a] / V;
        double G_eps_a = (3.0 * V * (P_axis - P_target_) + mtk_term) / W_;
        eps_dot_[a] += G_eps_a * half_dt;
    }
}

// ── Velocity scale by baro ──────────────────────────────────────────

void NoseHooverNPT::scale_v_by_baro_(const LaunchContext& ctx, double half_dt) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();

    double sum_eps_dot = 0.0;
    for (int a = 0; a < 3; ++a)
        if (couple_axis_[a])
            sum_eps_dot += eps_dot_[a];
    const int N_atoms = (n_dof_ + 3) / 3;
    double mtk_term2 =
        (pdim_ > 0) ? sum_eps_dot / (static_cast<double>(pdim_) * N_atoms)
                    : 0.0;

    double sx = std::exp(-(eps_dot_[0] + mtk_term2) * half_dt);
    double sy = std::exp(-(eps_dot_[1] + mtk_term2) * half_dt);
    double sz = std::exp(-(eps_dot_[2] + mtk_term2) * half_dt);

    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(n);
    if (sx == sy && sy == sz) {
        scale_vel_uniform_k<<<grid, block, 0, stream>>>(
            f.atoms.vx, f.atoms.vy, f.atoms.vz, static_cast<work_t>(sx), n);
        KE_global_ *= sx * sx;
    } else {
        scale_vel_per_axis_k<<<grid, block, 0, stream>>>(
            f.atoms.vx, f.atoms.vy, f.atoms.vz, static_cast<work_t>(sx),
            static_cast<work_t>(sy), static_cast<work_t>(sz), n);
        KE_global_ = compute_global_KE_(ctx);
    }
}

// ── Force kick ──────────────────────────────────────────────────────

void NoseHooverNPT::kick_v_by_force_(const LaunchContext& ctx, double half_dt) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(n);
    kick_v_k<<<grid, block, 0, stream>>>(
        f.atoms.vx, f.atoms.vy, f.atoms.vz, f.atoms.fx, f.atoms.fy, f.atoms.fz,
        f.atoms.mass, static_cast<work_t>(half_dt), n);
}

// ── Drift + box scale ───────────────────────────────────────────────

void NoseHooverNPT::propagate_x_and_box_(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();

    double eps_dt[3], a_arr[3], drift[3];
    for (int d = 0; d < 3; ++d) {
        eps_dt[d] = std::max(-kEpsDtMaxFraction,
                             std::min(eps_dot_[d] * dt_, kEpsDtMaxFraction));
        a_arr[d] = std::exp(eps_dt[d]);
        drift[d] = expm1c(eps_dt[d] * 0.5);
    }

    // Refresh cached box from device — a cooperative component bound earlier
    // in `initial_integrate` (e.g. `cuda::md::Deform`, which runs at order
    // -100 vs this integrator's order 0) may have rescaled f.box[0,4,8]
    // since the last NPT step. Without this D2H the cached values would
    // overwrite Deform's modification: stage-5 tension showed box_x frozen
    // at the initial 64.665 Å and strain stuck at 0 instead of ramping to
    // 1.7 over 172 000 steps. The sync is one D2H of three doubles per
    // step (negligible vs the force eval cost) and only matters when a
    // peer integrator-phase component mutates the box; otherwise it's a
    // round-trip of the value we cached on the previous step.
    cudaStreamSynchronize(stream);
    prec_t box_diag_dev[3] = {prec_t(0), prec_t(0), prec_t(0)};
    cudaMemcpy(&box_diag_dev[0], &f.box[0], sizeof(prec_t),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&box_diag_dev[1], &f.box[4], sizeof(prec_t),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&box_diag_dev[2], &f.box[8], sizeof(prec_t),
               cudaMemcpyDeviceToHost);
    box_lx_ = static_cast<double>(box_diag_dev[0]);
    box_ly_ = static_cast<double>(box_diag_dev[1]);
    box_lz_ = static_cast<double>(box_diag_dev[2]);

    double lx_new = box_lx_ * a_arr[0];
    double ly_new = box_ly_ * a_arr[1];
    double lz_new = box_lz_ * a_arr[2];

    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(n);
    propagate_x_k<<<grid, block, 0, stream>>>(
        f.atoms.x, f.atoms.y, f.atoms.z, f.atoms.vx, f.atoms.vy, f.atoms.vz,
        static_cast<work_t>(a_arr[0]), static_cast<work_t>(a_arr[1]),
        static_cast<work_t>(a_arr[2]), static_cast<work_t>(drift[0]),
        static_cast<work_t>(drift[1]), static_cast<work_t>(drift[2]),
        static_cast<work_t>(dt_), f.box, static_cast<work_t>(lx_new),
        static_cast<work_t>(ly_new), static_cast<work_t>(lz_new), n);

    box_lx_ = lx_new;
    box_ly_ = ly_new;
    box_lz_ = lz_new;
}

// ── Conserved MTK Hamiltonian ───────────────────────────────────────

double NoseHooverNPT::extended_hamiltonian_(const LaunchContext& ctx) const {
    const auto& f = ctx.frame;
    // D2H read of energy + kinetic.
    double E_pot = 0.0, KE = 0.0;
    cudaMemcpy(&E_pot, &f.scalars->energy, sizeof(double),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&KE, &f.scalars->kinetic, sizeof(double),
               cudaMemcpyDeviceToHost);
    double H = KE + E_pot;
    const double W_axis = W_ / 3.0;
    for (int a = 0; a < 3; ++a)
        if (couple_axis_[a])
            H += 0.5 * W_axis * eps_dot_[a] * eps_dot_[a];
    for (int k = 0; k < mtchain_; ++k)
        H += 0.5 * Q_t_[k] * eta_t_dot_[k] * eta_t_dot_[k];
    for (int k = 0; k < mpchain_; ++k)
        H += 0.5 * Q_b_[k] * eta_b_dot_[k] * eta_b_dot_[k];
    H += static_cast<double>(n_dof_) * kT_ * eta_t_[0];
    for (int k = 1; k < mtchain_; ++k)
        H += kT_ * eta_t_[k];
    for (int k = 0; k < mpchain_; ++k)
        H += kT_ * eta_b_[k];
    H += P_target_ * (box_lx_ * box_ly_ * box_lz_);
    return H;
}

// ── launch() dispatch ───────────────────────────────────────────────

void NoseHooverNPT::launch(const LaunchContext& ctx) {
    const double half_dt = 0.5 * dt_;

    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate)) {
        compute_temp_target_(ctx.step);
        const double tau_T2 = tau_T_ * tau_T_;
        Q_t_[0] = static_cast<double>(n_dof_) * kT_ * tau_T2;
        for (int k = 1; k < mtchain_; ++k)
            Q_t_[k] = kT_ * tau_T2;

        nhc_baro_half_(half_dt);          // P1
        nhc_atoms_half_(ctx, half_dt);    // P2
        update_eps_dot_(ctx, half_dt);    // P3 (uses prev step's stress)
        scale_v_by_baro_(ctx, half_dt);   // P4
        kick_v_by_force_(ctx, half_dt);   // P5
        propagate_x_and_box_(ctx);        // P6

    } else if (ctx.stage
               == static_cast<int>(ATV::md::Phase::final_integrate)) {
        kick_v_by_force_(ctx, half_dt);   // P7

    } else if (ctx.stage == static_cast<int>(ATV::md::Phase::end_of_step)) {
        // ReduceEnergy has run (single-GPU: no-op, but the order matches CPU).
        KE_global_ = compute_global_KE_(ctx);

        scale_v_by_baro_(ctx, half_dt);   // P8
        update_eps_dot_(ctx, half_dt);    // P9
        nhc_atoms_half_(ctx, half_dt);    // P10
        nhc_baro_half_(half_dt);          // P11

        // Rewrite frame.scalars.kinetic so downstream observers see
        // the post-NPT KE. D2H write.
        cudaMemcpy(&ctx.frame.scalars->kinetic, &KE_global_, sizeof(double),
                   cudaMemcpyHostToDevice);
        double H = extended_hamiltonian_(ctx);
        cudaMemcpy(&ctx.frame.scalars->conserved, &H, sizeof(double),
                   cudaMemcpyHostToDevice);
    }
}

} // namespace md
} // namespace cuda
} // namespace ATV
