/// @file cuda/scf/convergence_check.cu
/// @brief ConvergenceCheck: ||P - P_old||_F, energy = trace(H*P), convergence flag.

#include "atomiverse/cuda/scf/convergence_check.cuh"
#include "cublas_dispatch.cuh"

#include <cuda_runtime.h>

namespace ATV {
namespace scf {

void ConvergenceCheck::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);

    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasSetStream(handle, stream);
    cublas_handle_ = handle;
}

void ConvergenceCheck::launch(const LaunchContext& ctx) {
    auto& scf = ctx.frame.scf;
    if (scf.converged) return;

    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto handle = static_cast<cublasHandle_t>(cublas_handle_);
    Int n = scf.n;
    Int nn = n * n;

    // work = P - P_old
    cudaMemcpyAsync(scf.work, scf.P, nn * sizeof(Real),
                    cudaMemcpyDeviceToDevice, stream);
    Real neg_one = Real(-1);
    detail::axpy(handle, nn, &neg_one, scf.P_old, 1, scf.work, 1);

    // ||P - P_old||_F = nrm2(work)
    Real norm = 0;
    detail::nrm2(handle, nn, scf.work, 1, &norm);

    // Optionally compute energy = dot(H, P) (toy-scf mode)
    if (compute_energy_) {
        Real energy = 0;
        detail::dot(handle, nn, scf.H, 1, scf.P, 1, &energy);
        cudaMemcpyAsync(ctx.frame.energy, &energy, sizeof(Real),
                        cudaMemcpyHostToDevice, stream);
    }
    cudaStreamSynchronize(stream);

    // Update host-cached convergence state
    scf.residual = norm;
    scf.iteration = ctx.step;

    if (ctx.step > 0 && norm < tol_) {
        scf.converged = true;
    }
}

void ConvergenceCheck::destroy(const DestroyContext& /*ctx*/) {
    if (cublas_handle_) {
        cublasDestroy(static_cast<cublasHandle_t>(cublas_handle_));
        cublas_handle_ = nullptr;
    }
}

} // namespace scf
} // namespace ATV
