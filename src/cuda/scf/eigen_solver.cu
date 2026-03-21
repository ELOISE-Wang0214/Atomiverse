/// @file cuda/scf/eigen_solver.cu
/// @brief EigenSolver: dense symmetric eigendecomposition via cuSOLVER syevd.

#include "atomiverse/cuda/scf/eigen_solver.cuh"
#include "cublas_dispatch.cuh"

#include <cuda_runtime.h>

namespace ATV {
namespace scf {

void EigenSolver::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    Int n = ctx.frame.scf.n;

    // Create cuSOLVER handle
    cusolverDnHandle_t handle;
    cusolverDnCreate(&handle);
    cusolverDnSetStream(handle, stream);
    cusolver_handle_ = handle;

    // Query workspace size for syevd
    int lwork = 0;
    detail::syevd_bufferSize(
        handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
        n, static_cast<Real*>(nullptr), n, static_cast<Real*>(nullptr), &lwork);

    work_size_ = lwork;
    cudaMalloc(&d_work_, lwork * sizeof(Real));
    cudaMalloc(&d_info_, sizeof(int));

    cudaStreamSynchronize(stream);
}

void EigenSolver::launch(const LaunchContext& ctx) {
    auto& scf = ctx.frame.scf;
    if (scf.converged) return;

    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto handle = static_cast<cusolverDnHandle_t>(cusolver_handle_);
    Int n = scf.n;

    // Copy H -> C (syevd solves in-place, overwrites input with eigenvectors)
    cudaMemcpyAsync(scf.C, scf.H, n * n * sizeof(Real),
                    cudaMemcpyDeviceToDevice, stream);

    // Solve: C -> eigenvectors (column-major), eig -> eigenvalues (ascending)
    detail::syevd(
        handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
        n, scf.C, n, scf.eig, d_work_, work_size_, d_info_);
}

void EigenSolver::destroy(const DestroyContext& /*ctx*/) {
    if (cusolver_handle_) {
        cusolverDnDestroy(static_cast<cusolverDnHandle_t>(cusolver_handle_));
        cusolver_handle_ = nullptr;
    }
    cudaFree(d_work_);  d_work_ = nullptr;
    cudaFree(d_info_);  d_info_ = nullptr;
}

} // namespace scf
} // namespace ATV
