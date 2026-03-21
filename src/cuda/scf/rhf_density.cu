/// @file cuda/scf/rhf_density.cu
/// @brief RHFDensity: C = X*C', P = 2*C_occ*C_occ^T, E = 0.5*Tr[(Hcore+F)*P]+V_nuc.

#ifdef ATV_HAS_RHF

#include "atomiverse/cuda/scf/rhf_density.cuh"
#include "cublas_dispatch.cuh"

#include <cuda_runtime.h>

namespace ATV {
namespace scf {

void RHFDensity::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& scf = ctx.frame.scf;
    Int n = scf.n;

    // Store n_occ in workspace
    scf.n_occ = n_occ_;

    // Allocate Hcore+F temporary
    cudaMalloc(&d_HplusF_, n * n * sizeof(Real));

    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasSetStream(handle, stream);
    cublas_handle_ = handle;
}

void RHFDensity::launch(const LaunchContext& ctx) {
    auto& scf = ctx.frame.scf;
    if (scf.converged) return;

    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto handle = static_cast<cublasHandle_t>(cublas_handle_);
    Int n = scf.n;
    Int nn = n * n;

    // After EigenSolver: scf.C contains eigenvectors of F' (in orthogonal basis),
    // scf.eig contains eigenvalues. scf.work contains the real Fock matrix F.

    // Step 1: Back-transform C = X * C'
    // C' is in scf.C (from EigenSolver). We need C = X * C'.
    // Use d_HplusF_ as temp to hold C' before overwriting scf.C.
    cudaMemcpyAsync(d_HplusF_, scf.C, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    Real one = Real(1);
    Real zero = Real(0);
    detail::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                 n, n, n, &one, scf.X, n, d_HplusF_, n, &zero, scf.C, n);

    // Step 2: P_old = P; P = 2 * C_occ * C_occ^T
    cudaMemcpyAsync(scf.P_old, scf.P, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    Real two = Real(2);
    detail::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_T,
                 n, n, n_occ_,
                 &two, scf.C, n, scf.C, n, &zero, scf.P, n);

    // Step 3: Energy = 0.5 * Tr[(Hcore + F) * P] + V_nuc
    // scf.work = F (saved by RHFFock). Hcore = scf.Hcore.
    // d_HplusF_ = Hcore + F
    cudaMemcpyAsync(d_HplusF_, scf.Hcore, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    detail::axpy(handle, nn, &one, scf.work, 1, d_HplusF_, 1);

    // E_elec = 0.5 * dot(HplusF, P)
    Real dot_val = 0;
    detail::dot(handle, nn, d_HplusF_, 1, scf.P, 1, &dot_val);
    Real energy = Real(0.5) * dot_val + scf.V_nuc;

    // Write to device energy scalar
    cudaMemcpyAsync(ctx.frame.energy, &energy, sizeof(Real),
                    cudaMemcpyHostToDevice, stream);
    cudaStreamSynchronize(stream);
}

void RHFDensity::destroy(const DestroyContext& /*ctx*/) {
    if (cublas_handle_) {
        cublasDestroy(static_cast<cublasHandle_t>(cublas_handle_));
        cublas_handle_ = nullptr;
    }
    cudaFree(d_HplusF_); d_HplusF_ = nullptr;
}

} // namespace scf
} // namespace ATV

#endif // ATV_HAS_RHF
