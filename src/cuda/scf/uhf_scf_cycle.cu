/// @file cuda/scf/uhf_scf_cycle.cu
/// @brief One-stage UHF helpers. The loop is cuda::scf::UHF.

#include "atomiverse/cuda/scf/internal/uhf_scf_cycle.cuh"

#include "atomiverse/cuda/device/launch.cuh"
#include "atomiverse/cuda/frame/types.cuh"
#include "atomiverse/cuda/scf/internal/diis_mixer.cuh"
#include "atomiverse/cuda/scf/internal/eigen_solver.cuh"
#include "atomiverse/cuda/scf/internal/uhf_density.cuh"
#include "atomiverse/cuda/scf/internal/uhf_fock.cuh"
#include "cublas_dispatch.cuh"

#include <cmath>
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <utility>

namespace ATV {
namespace cuda {
namespace scf {
namespace internal {
namespace {

__global__ void uhf_deshift_virtuals_kernel(prec_t* __restrict__ eig,
                                            prec_t lambda, int n_occ, int n) {
    const int i =
        static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x) + n_occ;
    if (i < n)
        eig[i] -= lambda;
}

__global__ void uhf_residual_kernel(prec_t* __restrict__ residual,
                                    int* __restrict__ flag,
                                    const prec_t* __restrict__ norms,
                                    prec_t tol, int step) {
    const prec_t r = sqrt((*norms) * (*norms) + norms[1] * norms[1]);
    *residual = r;
    if (step > 0 && r < tol)
        *flag = 1;
}

void swap_beta(SCFWorkspace& ws, UhfBeta& b) {
    std::swap(ws.H, b.H);
    std::swap(ws.P, b.P);
    std::swap(ws.P_old, b.P_old);
    std::swap(ws.C, b.C);
    std::swap(ws.C_prime, b.C_prime);
    std::swap(ws.eig, b.eig);
    std::swap(ws.work, b.work);
    std::swap(ws.F_prime, b.F_prime);
    std::swap(ws.n_occ, b.n_occ);
}

struct BetaGuard {
    SCFWorkspace& ws;
    UhfBeta& b;
    explicit BetaGuard(SCFWorkspace& ws_in, UhfBeta& b_in)
        : ws(ws_in), b(b_in) {
        swap_beta(ws, b);
    }
    ~BetaGuard() { swap_beta(ws, b); }
};

void free_beta(UhfBeta& b) {
    cudaFree(b.H);
    cudaFree(b.P);
    cudaFree(b.P_old);
    cudaFree(b.C);
    cudaFree(b.C_prime);
    cudaFree(b.eig);
    cudaFree(b.work);
    cudaFree(b.F_prime);
    b = {};
}

void free_workspace(SCFWorkspace& ws) {
    cudaFree(ws.H0);
    cudaFree(ws.H);
    cudaFree(ws.P);
    cudaFree(ws.P_old);
    cudaFree(ws.C);
    cudaFree(ws.eig);
    cudaFree(ws.work);
    cudaFree(ws.Hcore);
    cudaFree(ws.S);
    cudaFree(ws.X);
    cudaFree(ws.F_prime);
    cudaFree(ws.ERI);
    cudaFree(ws.d_shell_alpha);
    cudaFree(ws.d_shell_coeff);
    cudaFree(ws.d_shell_center);
    cudaFree(ws.d_shell_atom);
    cudaFree(ws.d_shell_L);
    cudaFree(ws.d_ao_shell);
    cudaFree(ws.d_ao_cart);
    cudaFree(ws.d_atom_Z);
    cudaFree(ws.d_atom_pos);
    cudaFree(ws.d_V_nuc);
    cudaFree(ws.d_residual);
    if (ws.h_converged_mapped) {
        cudaFreeHost(ws.h_converged_mapped);
        ws.h_converged_mapped = nullptr;
        ws.d_converged = nullptr;
    }
    cudaFree(ws.d_scalar_scratch);
    cudaFree(ws.d_W);
    cudaFree(ws.C_prime);
    ws = {};
}

} // namespace

UhfScfCycle::UhfScfCycle(int n_shells, int n_ao, int n_prim,
                         const double* alpha_mol, const double* coeff_mol,
                         const int* L_mol, const int* atom_of_mol, int n_alpha,
                         int n_beta, double tol, int max_iter,
                         double level_shift)
    : tol_(tol), level_shift_(level_shift) {
    (void)max_iter;
    if (n_alpha < 0 || n_beta < 0) {
        throw std::invalid_argument(
            "cuda::scf::UHF: n_alpha and n_beta must be >= 0 (got "
            + std::to_string(n_alpha) + ", " + std::to_string(n_beta) + ").");
    }
    if (n_shells < 1 || n_ao < 1 || n_prim < 1) {
        throw std::invalid_argument(
            "cuda::scf::UHF: n_shells, n_ao, n_prim must be >= 1");
    }
    eigen_ = std::make_unique<EigenSolver>();
    diis_alpha_ = std::make_unique<DIISMixer>(8, 1, level_shift_);
    diis_beta_ = std::make_unique<DIISMixer>(8, 1, level_shift_);
    fock_ = std::make_unique<UHFFockDevice>(n_shells, n_ao, n_prim, alpha_mol,
                                            coeff_mol, L_mol, atom_of_mol,
                                            n_alpha, n_beta);
    density_ = std::make_unique<UHFDensity>();
}

UhfScfCycle::~UhfScfCycle() = default;

void UhfScfCycle::build(void* stream, cuda::Frame& frame) {
    SolverContext ctx{stream, frame, 0};
    fock_->build(ctx, workspace_, beta_);
    diis_alpha_->build(ctx, workspace_);
    {
        BetaGuard guard(workspace_, beta_);
        diis_beta_->build(ctx, workspace_);
    }
    eigen_->build(ctx, workspace_);
    density_->build(ctx, workspace_);
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasSetStream(handle, static_cast<cudaStream_t>(stream));
    check_blas_ = handle;
    built_ = true;
}

void UhfScfCycle::begin(void* stream) {
    if (!built_) {
        throw std::runtime_error(
            "cuda::scf::UHF: UhfScfCycle::begin before build()");
    }
    auto cu = static_cast<cudaStream_t>(stream);
    workspace_.converged = false;
    workspace_.iteration = 0;
    if (workspace_.h_converged_mapped)
        *workspace_.h_converged_mapped = 0;
    const size_t bytes =
        static_cast<size_t>(workspace_.n * workspace_.n) * sizeof(prec_t);
    if (workspace_.P)
        cudaMemsetAsync(workspace_.P, 0, bytes, cu);
    if (beta_.P)
        cudaMemsetAsync(beta_.P, 0, bytes, cu);
    if (diis_alpha_)
        diis_alpha_->reset();
    if (diis_beta_)
        diis_beta_->reset();
}

void UhfScfCycle::step_fock(void* stream, cuda::Frame& frame, int step) {
    SolverContext ctx{stream, frame, step};
    fock_->step(ctx, workspace_, beta_);
}

void UhfScfCycle::step_mix(void* stream, cuda::Frame& frame, int step) {
    SolverContext ctx{stream, frame, step};
    diis_alpha_->step(ctx, workspace_);
    BetaGuard guard(workspace_, beta_);
    diis_beta_->step(ctx, workspace_);
}

void UhfScfCycle::step_diagonalize(void* stream, cuda::Frame& frame, int step) {
    SolverContext ctx{stream, frame, step};
    eigen_->step(ctx, workspace_);
    BetaGuard guard(workspace_, beta_);
    eigen_->step(ctx, workspace_);
}

void UhfScfCycle::step_density(void* stream, cuda::Frame& frame, int step) {
    SolverContext ctx{stream, frame, step};
    density_->step(ctx, workspace_, beta_);
}

void UhfScfCycle::step_check(void* stream, cuda::Frame& frame, int step) {
    (void)frame;
    auto& ws = workspace_;
    if (ws.converged || ws.n <= 0)
        return;
    auto cu = static_cast<cudaStream_t>(stream);
    auto handle = static_cast<cublasHandle_t>(check_blas_);
    cublasSetStream(handle, cu);
    const int nn = ws.n * ws.n;
    const size_t bytes = static_cast<size_t>(nn) * sizeof(prec_t);
    cudaMemcpyAsync(ws.work, ws.P, bytes, cudaMemcpyDeviceToDevice, cu);
    cudaMemcpyAsync(beta_.work, beta_.P, bytes, cudaMemcpyDeviceToDevice, cu);
    cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
    const prec_t neg = prec_t(-1);
    detail::axpy(handle, nn, &neg, ws.P_old, 1, ws.work, 1);
    detail::axpy(handle, nn, &neg, beta_.P_old, 1, beta_.work, 1);
    cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
    detail::nrm2(handle, nn, ws.work, 1, ws.d_scalar_scratch);
    detail::nrm2(handle, nn, beta_.work, 1, ws.d_scalar_scratch + 1);
    cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
    uhf_residual_kernel<<<1, 1, 0, cu>>>(ws.d_residual, ws.d_converged,
                                         ws.d_scalar_scratch,
                                         static_cast<prec_t>(tol_), step);
}

bool UhfScfCycle::poll_converged(void* stream, int step) {
    cudaStreamSynchronize(static_cast<cudaStream_t>(stream));
    if (workspace_.h_converged_mapped)
        workspace_.converged = (*workspace_.h_converged_mapped != 0);
    if (workspace_.converged)
        workspace_.iteration = step;
    return workspace_.converged;
}

void UhfScfCycle::finish(void* stream) {
    auto cu = static_cast<cudaStream_t>(stream);
    auto& ws = workspace_;
    if (ws.h_converged_mapped)
        ws.converged = (*ws.h_converged_mapped != 0);
    if (ws.d_residual) {
        prec_t h_res = 0;
        cudaMemcpyAsync(&h_res, ws.d_residual, sizeof(prec_t),
                        cudaMemcpyDeviceToHost, cu);
        cudaStreamSynchronize(cu);
        ws.residual = static_cast<double>(h_res);
    }
    if (ws.d_V_nuc) {
        cudaMemcpyAsync(&ws.V_nuc, ws.d_V_nuc, sizeof(double),
                        cudaMemcpyDeviceToHost, cu);
        cudaStreamSynchronize(cu);
    }
    if (ws.converged && std::abs(level_shift_) > 1e-4) {
        auto deshift = [&](prec_t* eig, int n_occ) {
            if (!eig || n_occ >= ws.n || n_occ <= 0)
                return;
            const int n_virt = ws.n - n_occ;
            uhf_deshift_virtuals_kernel<<<cuda::grid_size(n_virt,
                                                          cuda::kDefaultBlock),
                                          cuda::kDefaultBlock, 0, cu>>>(
                eig, static_cast<prec_t>(level_shift_), n_occ, ws.n);
        };
        deshift(ws.eig, ws.n_occ);
        deshift(beta_.eig, beta_.n_occ);
        cudaStreamSynchronize(cu);
    }
}

void UhfScfCycle::destroy(void* stream, cuda::Frame& frame) {
    SolverContext ctx{stream, frame, 0};
    if (density_)
        density_->destroy(ctx, workspace_);
    if (eigen_)
        eigen_->destroy(ctx, workspace_);
    if (diis_beta_) {
        BetaGuard guard(workspace_, beta_);
        diis_beta_->destroy(ctx, workspace_);
    }
    if (diis_alpha_)
        diis_alpha_->destroy(ctx, workspace_);
    if (fock_)
        fock_->destroy(ctx, workspace_);
    if (check_blas_) {
        cublasDestroy(static_cast<cublasHandle_t>(check_blas_));
        check_blas_ = nullptr;
    }
    free_beta(beta_);
    free_workspace(workspace_);
    built_ = false;
}

} // namespace internal
} // namespace scf
} // namespace cuda
} // namespace ATV
