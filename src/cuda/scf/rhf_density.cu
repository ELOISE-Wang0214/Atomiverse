/// @file cuda/scf/rhf_density.cu
/// @brief internal::RHFDensity: C = X*C', P = 2*C_occ*C_occ^T,
///        E = 0.5*Tr[(Hcore+F)*P]+V_nuc — fully on device (no D2H).

#ifdef ATV_HAS_HF

#include "atomiverse/cuda/scf/internal/rhf_density.cuh"
#include "cublas_dispatch.cuh"

#include "atomiverse/core/units.h"
#include "atomiverse/cuda/frame/types.cuh"
#include <cublas_v2.h>
#include <cuda_runtime.h>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

namespace cuda {
namespace scf {
namespace internal {

/// energy_ev = (0.5 * dot + V_nuc) * Ha→eV.
/// `V_nuc` is passed by value (host mirror filled at end of device Fock
/// build) so energy never depends on a mid-SCF D2H of the scalar.
__global__ void write_rhf_energy_ev_kernel(prec_t* __restrict__ energy_ev,
                                           const prec_t* __restrict__ dot_HFP,
                                           double V_nuc, prec_t ha2ev) {
    const prec_t e_ha =
        prec_t(0.5) * (*dot_HFP) + static_cast<prec_t>(V_nuc);
    *energy_ev = e_ha * ha2ev;
}

void RHFDensity::build(const SolverContext& ctx, SCFWorkspace& ws) {
    auto stream = static_cast<cudaStream_t>(ctx.stream);
    auto& scf = ws;
    int n = scf.n;

    scf.n_occ = n_occ_;

    cudaMalloc(&d_HplusF_, n * n * sizeof(prec_t));

    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasSetStream(handle, stream);
    // Device pointer mode: dot result stays on GPU (no implicit D2H).
    cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
    cublas_handle_ = handle;
}

void RHFDensity::step(const SolverContext& ctx, SCFWorkspace& ws) {
    auto& scf = ws;
    if (scf.converged)
        return;

    auto stream = static_cast<cudaStream_t>(ctx.stream);
    auto handle = static_cast<cublasHandle_t>(cublas_handle_);
    int n = scf.n;
    int nn = n * n;

    cudaMemcpyAsync(scf.C_prime, scf.C, nn * sizeof(prec_t),
                    cudaMemcpyDeviceToDevice, stream);

    cudaMemcpyAsync(d_HplusF_, scf.C, nn * sizeof(prec_t),
                    cudaMemcpyDeviceToDevice, stream);
    // gemm needs host-side alpha/beta scalars — cublas device pointer mode
    // only applies to result of reductions; alpha/beta stay host for gemm.
    // Temporarily switch to HOST for gemm constants, then back to DEVICE.
    cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
    prec_t one = prec_t(1);
    prec_t zero = prec_t(0);
    detail::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, scf.X, n,
                 d_HplusF_, n, &zero, scf.C, n);

    cudaMemcpyAsync(scf.P_old, scf.P, nn * sizeof(prec_t),
                    cudaMemcpyDeviceToDevice, stream);
    prec_t two = prec_t(2);
    detail::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, n, n, n_occ_, &two, scf.C, n,
                 scf.C, n, &zero, scf.P, n);

    // Hcore + F → d_HplusF_
    cudaMemcpyAsync(d_HplusF_, scf.Hcore, nn * sizeof(prec_t),
                    cudaMemcpyDeviceToDevice, stream);
    detail::axpy(handle, nn, &one, scf.work, 1, d_HplusF_, 1);

    // E = 0.5 Tr[(Hcore+F)P] + V_nuc → frame energy [eV]. Prefer fully
    // device path (no host DOT sync) when d_scalar_scratch is available.
    if (ctx.frame.scalars) {
        if (scf.d_scalar_scratch) {
            cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
            detail::dot(handle, nn, d_HplusF_, 1, scf.P, 1,
                        scf.d_scalar_scratch);
            write_rhf_energy_ev_kernel<<<1, 1, 0, stream>>>(
                &ctx.frame.scalars->energy, scf.d_scalar_scratch, scf.V_nuc,
                static_cast<prec_t>(units::kHartreeInEv));
            cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
        } else {
            cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
            prec_t dot_val = 0;
            detail::dot(handle, nn, d_HplusF_, 1, scf.P, 1, &dot_val);
            prec_t energy_ha =
                prec_t(0.5) * dot_val + static_cast<prec_t>(scf.V_nuc);
            prec_t energy_ev =
                energy_ha * static_cast<prec_t>(units::kHartreeInEv);
            cudaMemcpyAsync(&ctx.frame.scalars->energy, &energy_ev,
                            sizeof(prec_t), cudaMemcpyHostToDevice, stream);
        }
    }
}

void RHFDensity::destroy(const SolverContext& /*ctx*/, SCFWorkspace& /*ws*/) {
    if (cublas_handle_) {
        cublasDestroy(static_cast<cublasHandle_t>(cublas_handle_));
        cublas_handle_ = nullptr;
    }
    cudaFree(d_HplusF_);
    d_HplusF_ = nullptr;
}

} // namespace internal
} // namespace scf
} // namespace cuda
} // namespace ATV

#endif // ATV_HAS_HF
