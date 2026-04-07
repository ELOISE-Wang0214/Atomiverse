/// @file cuda/scf/hamiltonian.cu
/// @brief Hamiltonian: allocate SCF workspace + build H = H0 + alpha * P.

#include "atomiverse/cuda/scf/hamiltonian.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "cublas_dispatch.cuh"

#include <cuda_runtime.h>
#include "atomiverse/cuda/types.cuh"

namespace ATV {

using cuda::real_t;
using cuda::work_t;
using cuda::prec_t;

namespace scf {

// ── Kernel: initialize H0 as a deterministic symmetric matrix ────────

/// Generate a diagonally-dominant symmetric H0.
/// H0[i,j] = 1.0 / (1 + |i-j|)  for off-diagonal
/// H0[i,i] = n + 1 - i           for diagonal (ensures positive-definite)
/// Column-major layout.
__global__ void init_h0_kernel(Real* H0, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n * n;
    for (int k = idx; k < total; k += blockDim.x * gridDim.x) {
        int col = k / n;
        int row = k % n;
        if (row == col) {
            H0[k] = work_t(n + 1 - row);
        } else {
            int diff = row > col ? row - col : col - row;
            H0[k] = work_t(1) / work_t(1 + diff);
        }
    }
}

// ── Lifecycle ────────────────────────────────────────────────────────

void Hamiltonian::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& scf = ctx.frame.scf;
    int n = ctx.frame.n_atoms();

    // Store metadata
    scf.n = n;
    scf.n_occ = n_occ_;
    scf.converged = false;
    scf.iteration = 0;
    scf.residual = work_t(0);

    int nn = n * n;

    // Allocate all SCF workspace matrices
    cudaMalloc(&scf.H0,    nn * sizeof(work_t));
    cudaMalloc(&scf.H,     nn * sizeof(work_t));
    cudaMalloc(&scf.P,     nn * sizeof(work_t));
    cudaMalloc(&scf.P_old, nn * sizeof(work_t));
    cudaMalloc(&scf.C,     nn * sizeof(work_t));
    cudaMalloc(&scf.eig,   n * sizeof(work_t));
    cudaMalloc(&scf.work,  nn * sizeof(work_t));

    // Zero-initialize P
    cudaMemsetAsync(scf.P, 0, nn * sizeof(work_t), stream);
    cudaMemsetAsync(scf.P_old, 0, nn * sizeof(work_t), stream);

    // Initialize H0
    int threads = cuda::kDefaultBlock;
    int blocks = cuda::grid_size(nn, threads);
    init_h0_kernel<<<blocks, threads, 0, stream>>>(scf.H0, n);

    // Create cuBLAS handle
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasSetStream(handle, stream);
    cublas_handle_ = handle;

    cudaStreamSynchronize(stream);
}

void Hamiltonian::launch(const LaunchContext& ctx) {
    auto& scf = ctx.frame.scf;
    if (scf.converged) return;

    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto handle = static_cast<cublasHandle_t>(cublas_handle_);
    int n = scf.n;
    int nn = n * n;

    // H = H0 (copy)
    cudaMemcpyAsync(scf.H, scf.H0, nn * sizeof(work_t), cudaMemcpyDeviceToDevice, stream);

    // H += alpha * P
    work_t alpha = alpha_;
    detail::axpy(handle, nn, &alpha, scf.P, 1, scf.H, 1);
}

void Hamiltonian::destroy(const DestroyContext& /*ctx*/) {
    if (cublas_handle_) {
        cublasDestroy(static_cast<cublasHandle_t>(cublas_handle_));
        cublas_handle_ = nullptr;
    }
}

} // namespace scf
} // namespace ATV
