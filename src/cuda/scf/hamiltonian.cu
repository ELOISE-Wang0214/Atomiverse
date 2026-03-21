/// @file cuda/scf/hamiltonian.cu
/// @brief Hamiltonian: allocate SCF workspace + build H = H0 + alpha * P.

#include "atomiverse/cuda/scf/hamiltonian.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "cublas_dispatch.cuh"

#include <cuda_runtime.h>

namespace ATV {
namespace scf {

// ── Kernel: initialize H0 as a deterministic symmetric matrix ────────

/// Generate a diagonally-dominant symmetric H0.
/// H0[i,j] = 1.0 / (1 + |i-j|)  for off-diagonal
/// H0[i,i] = n + 1 - i           for diagonal (ensures positive-definite)
/// Column-major layout.
__global__ void init_h0_kernel(Real* H0, Int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n * n;
    for (int k = idx; k < total; k += blockDim.x * gridDim.x) {
        int col = k / n;
        int row = k % n;
        if (row == col) {
            H0[k] = Real(n + 1 - row);
        } else {
            int diff = row > col ? row - col : col - row;
            H0[k] = Real(1) / Real(1 + diff);
        }
    }
}

// ── Lifecycle ────────────────────────────────────────────────────────

void Hamiltonian::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& scf = ctx.frame.scf;
    Int n = ctx.frame.n_atoms();

    // Store metadata
    scf.n = n;
    scf.n_occ = n_occ_;
    scf.converged = false;
    scf.iteration = 0;
    scf.residual = Real(0);

    Int nn = n * n;

    // Allocate all SCF workspace matrices
    cudaMalloc(&scf.H0,    nn * sizeof(Real));
    cudaMalloc(&scf.H,     nn * sizeof(Real));
    cudaMalloc(&scf.P,     nn * sizeof(Real));
    cudaMalloc(&scf.P_old, nn * sizeof(Real));
    cudaMalloc(&scf.C,     nn * sizeof(Real));
    cudaMalloc(&scf.eig,   n * sizeof(Real));
    cudaMalloc(&scf.work,  nn * sizeof(Real));

    // Zero-initialize P
    cudaMemsetAsync(scf.P, 0, nn * sizeof(Real), stream);
    cudaMemsetAsync(scf.P_old, 0, nn * sizeof(Real), stream);

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
    Int n = scf.n;
    Int nn = n * n;

    // H = H0 (copy)
    cudaMemcpyAsync(scf.H, scf.H0, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);

    // H += alpha * P
    Real alpha = alpha_;
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
