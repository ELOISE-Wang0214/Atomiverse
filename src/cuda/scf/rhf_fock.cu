/// @file cuda/scf/rhf_fock.cu
/// @brief RHFFock: upload integrals, build G(P), F=Hcore+G, F'=X^T*F*X.

#ifdef ATV_HAS_RHF

#include "atomiverse/cuda/scf/rhf_fock.cuh"
#include "atomiverse/cuda/scf/device_integrals.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/core/gto.h"
#include "atomiverse/core/constants.h"
#include "cublas_dispatch.cuh"

#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <xtensor/xnorm.hpp>
#include <xtensor/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

namespace ATV {
namespace scf {

// ── Kernel: build G(P) from ERI ──────────────────────────────────────
// G_ij = sum_kl P_kl * [(ij|kl) - 0.5*(ik|jl)]
// P already includes factor 2 from RHF closed-shell occupancy.
__global__ void build_G_kernel(Real* __restrict__ G, const Real* __restrict__ P, const Real* __restrict__ ERI, Int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n * n;
    for (int ij = idx; ij < total; ij += blockDim.x * gridDim.x) {
        int j = ij / n;  // column-major: ij = i + j*n
        int i = ij % n;
        Real val = Real(0);
        for (int k = 0; k < n; ++k)
            for (int l = 0; l < n; ++l) {
                int ij_kl = i + j*n + k*n*n + l*n*n*n;
                int ik_jl = i + k*n + j*n*n + l*n*n*n;
                val += P[k + l*n] * (ERI[ij_kl] - Real(0.5) * ERI[ik_jl]);
            }
        G[ij] = val;
    }
}

void RHFFock::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& scf = ctx.frame.scf;
    Int n = integrals_.n_basis;
    Int nn = n * n;
    Int n4 = nn * nn;

    // Set workspace metadata
    scf.n = n;
    scf.n_occ = integrals_.n_occ();
    scf.V_nuc = integrals_.V_nuc;
    scf.converged = false;

    // Allocate device matrices
    cudaMalloc(&scf.H,       nn * sizeof(Real));
    cudaMalloc(&scf.P,       nn * sizeof(Real));
    cudaMalloc(&scf.P_old,   nn * sizeof(Real));
    cudaMalloc(&scf.C,       nn * sizeof(Real));
    cudaMalloc(&scf.eig,     n * sizeof(Real));
    cudaMalloc(&scf.work,    nn * sizeof(Real));
    cudaMalloc(&scf.Hcore,   nn * sizeof(Real));
    cudaMalloc(&scf.X,       nn * sizeof(Real));
    cudaMalloc(&scf.F_prime, nn * sizeof(Real));
    cudaMalloc(&scf.ERI,     n4 * sizeof(Real));
    cudaMalloc(&d_temp_,     nn * sizeof(Real));

    // Zero P
    cudaMemsetAsync(scf.P, 0, nn * sizeof(Real), stream);
    cudaMemsetAsync(scf.P_old, 0, nn * sizeof(Real), stream);

    // Upload Hcore (row-major xtensor -> column-major device)
    // For symmetric matrices, row-major == column-major, so direct copy works.
    std::vector<Real> buf(nn);
    for (Int i = 0; i < n; ++i)
        for (Int j = 0; j < n; ++j)
            buf[i + j * n] = integrals_.Hcore(i, j);  // col-major
    cudaMemcpy(scf.Hcore, buf.data(), nn * sizeof(Real), cudaMemcpyHostToDevice);

    // Upload ERI (row-major 4D xtensor -> column-major 4D device)
    std::vector<Real> eri_buf(n4);
    for (Int i = 0; i < n; ++i)
        for (Int j = 0; j < n; ++j)
            for (Int k = 0; k < n; ++k)
                for (Int l = 0; l < n; ++l)
                    eri_buf[i + j*n + k*nn + l*nn*n] = integrals_.ERI(i, j, k, l);
    cudaMemcpy(scf.ERI, eri_buf.data(), n4 * sizeof(Real), cudaMemcpyHostToDevice);

    // Compute X = S^{-1/2} on host using xtensor-blas
    auto [s_eig, U] = xt::linalg::eigh(integrals_.S);
    xt::xtensor<Real, 2> s_inv_sqrt = xt::zeros<Real>({n, n});
    for (Int i = 0; i < n; ++i)
        s_inv_sqrt(i, i) = Real(1) / std::sqrt(s_eig(i));
    auto X_host = xt::linalg::dot(xt::linalg::dot(U, s_inv_sqrt), xt::transpose(U));

    // Upload X (col-major)
    for (Int i = 0; i < n; ++i)
        for (Int j = 0; j < n; ++j)
            buf[i + j * n] = X_host(i, j);
    cudaMemcpy(scf.X, buf.data(), nn * sizeof(Real), cudaMemcpyHostToDevice);

    // cuBLAS handle
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasSetStream(handle, stream);
    cublas_handle_ = handle;

    cudaStreamSynchronize(stream);
}

void RHFFock::launch(const LaunchContext& ctx) {
    auto& scf = ctx.frame.scf;
    if (scf.converged) return;

    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto handle = static_cast<cublasHandle_t>(cublas_handle_);
    Int n = scf.n;
    Int nn = n * n;

    // Step 1: G(P) from ERI
    // Use work buffer for G
    cudaMemsetAsync(scf.work, 0, nn * sizeof(Real), stream);
    int threads = cuda::kDefaultBlock;
    int blocks = cuda::grid_size(nn, threads);
    build_G_kernel<<<blocks, threads, 0, stream>>>(scf.work, scf.P, scf.ERI, n);

    // Step 2: F = Hcore + G (H serves as F)
    cudaMemcpyAsync(scf.H, scf.Hcore, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    Real one = Real(1);
    detail::axpy(handle, nn, &one, scf.work, 1, scf.H, 1);

    // Step 3: F' = X^T * F * X
    // temp = F * X
    Real zero = Real(0);
    detail::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                 n, n, n, &one, scf.H, n, scf.X, n, &zero, d_temp_, n);
    // F' = X^T * temp
    detail::gemm(handle, CUBLAS_OP_T, CUBLAS_OP_N,
                 n, n, n, &one, scf.X, n, d_temp_, n, &zero, scf.F_prime, n);

    // EigenSolver will diagonalize F_prime (stored where H usually goes for toy-scf).
    // Copy F' -> H so EigenSolver picks it up via scf.H
    // Actually, EigenSolver reads from scf.H and writes to scf.C/eig.
    // We need it to read F_prime. Swap: copy F_prime into the slot EigenSolver reads.
    // EigenSolver copies scf.H -> scf.C then diag in-place. So set scf.H = F_prime.
    // But we need scf.H (Fock matrix F) later for energy. Save F in work temporarily.
    // Simpler: F_prime is separate. Copy F_prime -> H for EigenSolver.
    // After EigenSolver runs, H contains F_prime (overwritten by copy to C).
    // The real Fock F is in d_temp_ ... no, d_temp_ was overwritten.
    // Let me just swap: save F somewhere, put F_prime in H.
    // We can use P_old as temp (it's only used later in build_density).
    // Actually, simpler: EigenSolver reads scf.H. So we just overwrite scf.H with F_prime.
    // We need F for energy later → store it in work.
    cudaMemcpyAsync(scf.work, scf.H, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    cudaMemcpyAsync(scf.H, scf.F_prime, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    // Now: scf.H = F' (for EigenSolver), scf.work = F (for energy in RHFDensity)
}

void RHFFock::destroy(const DestroyContext& /*ctx*/) {
    if (cublas_handle_) {
        cublasDestroy(static_cast<cublasHandle_t>(cublas_handle_));
        cublas_handle_ = nullptr;
    }
    cudaFree(d_temp_); d_temp_ = nullptr;
}

// ══════════════════════════════════════════════════════════════════════
// RHFFockDevice — CUDA-only path: GPU integral kernels, zero CPU integrals
// ══════════════════════════════════════════════════════════════════════

// kAngstromToBohr is in core/constants.h

} // namespace scf (temporarily close for cuda::scf access)

static void call_1e_integrals(ATV::Real* S, ATV::Real* T, ATV::Real* V,
    const ATV::Real* a, const ATV::Real* c, const ATV::Real* ctr,
    const ATV::Int* Z, const ATV::Real* pos, int ns, int np, int na, void* st) {
    ATV::cuda::scf::launch_1e_integrals(S, T, V, a, c, ctr, Z, pos, ns, np, na, st);
}
static void call_eri(ATV::Real* ERI, const ATV::Real* a, const ATV::Real* c,
    const ATV::Real* ctr, int ns, int np, void* st) {
    ATV::cuda::scf::launch_eri(ERI, a, c, ctr, ns, np, st);
}

namespace scf { // reopen

RHFFockDevice::RHFFockDevice(int n_prim, const double* exponents,
                             const double* coefficients, int n_occ)
    : n_prim_(n_prim), n_occ_(n_occ),
      exponents_(exponents, exponents + n_prim),
      coefficients_(coefficients, coefficients + n_prim) {}

void RHFFockDevice::build(const BuildContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& scf = ctx.frame.scf;
    Int natm = ctx.frame.n_atoms();
    Int n = natm;  // 1 basis function per atom (STO-3G 1s for H)
    Int nn = n * n;
    Int n4 = nn * nn;

    scf.n = n;
    scf.n_occ = n_occ_;
    scf.converged = false;
    scf.n_shells = n;
    scf.n_prim = n_prim_;
    scf.natm = natm;

    // ── 1. Upload basis parameters to device (tiny: ~100 bytes) ───────
    // Build normalized coefficients with physical norm
    std::vector<Real> h_alpha(n * n_prim_), h_coeff(n * n_prim_);
    std::vector<Real> h_centers(n * 3), h_pos(natm * 3);
    std::vector<Int> h_Z(natm);

    // Read atom positions from cuda::Frame (they're on device — need D2H)
    // But positions are also in the host-cached Frame. We get them from the
    // host Frame that Backend::materialize() copied. The positions are stored
    // in ctx.frame (device). We'll reconstruct from the known atom count and
    // read from device.
    Real* h_x = new Real[natm]; Real* h_y = new Real[natm]; Real* h_z = new Real[natm];
    Int* h_type = new Int[natm];
    cudaMemcpy(h_x, ctx.frame.x, natm * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_y, ctx.frame.y, natm * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_z, ctx.frame.z, natm * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_type, ctx.frame.type_id, natm * sizeof(Int), cudaMemcpyDeviceToHost);

    for (Int a = 0; a < natm; ++a) {
        double cx = static_cast<double>(h_x[a]) * kAngstromToBohr;
        double cy = static_cast<double>(h_y[a]) * kAngstromToBohr;
        double cz = static_cast<double>(h_z[a]) * kAngstromToBohr;
        h_centers[a*3+0] = Real(cx); h_centers[a*3+1] = Real(cy); h_centers[a*3+2] = Real(cz);
        h_pos[a*3+0] = Real(cx); h_pos[a*3+1] = Real(cy); h_pos[a*3+2] = Real(cz);
        h_Z[a] = h_type[a];
        for (int p = 0; p < n_prim_; ++p) {
            h_alpha[a * n_prim_ + p] = Real(exponents_[p]);
            h_coeff[a * n_prim_ + p] = Real(coefficients_[p]
                * ATV::phys_norm_s(exponents_[p]));
        }
    }

    // Nuclear repulsion (host computation — trivial)
    scf.V_nuc = 0;
    for (Int a = 0; a < natm; ++a)
        for (Int b = a + 1; b < natm; ++b) {
            double dx = static_cast<double>(h_centers[a*3]) - static_cast<double>(h_centers[b*3]);
            double dy = static_cast<double>(h_centers[a*3+1]) - static_cast<double>(h_centers[b*3+1]);
            double dz = static_cast<double>(h_centers[a*3+2]) - static_cast<double>(h_centers[b*3+2]);
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            scf.V_nuc += Real(h_Z[a] * h_Z[b] / r);
        }

    delete[] h_x; delete[] h_y; delete[] h_z; delete[] h_type;

    // Upload basis data
    cudaMalloc(&scf.d_shell_alpha, n * n_prim_ * sizeof(Real));
    cudaMalloc(&scf.d_shell_coeff, n * n_prim_ * sizeof(Real));
    cudaMalloc(&scf.d_shell_center, n * 3 * sizeof(Real));
    cudaMalloc(&scf.d_atom_Z, natm * sizeof(Int));
    cudaMalloc(&scf.d_atom_pos, natm * 3 * sizeof(Real));
    cudaMemcpy(scf.d_shell_alpha, h_alpha.data(), n * n_prim_ * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(scf.d_shell_coeff, h_coeff.data(), n * n_prim_ * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(scf.d_shell_center, h_centers.data(), n * 3 * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(scf.d_atom_Z, h_Z.data(), natm * sizeof(Int), cudaMemcpyHostToDevice);
    cudaMemcpy(scf.d_atom_pos, h_pos.data(), natm * 3 * sizeof(Real), cudaMemcpyHostToDevice);

    // ── 2. Allocate SCF workspace matrices ────────────────────────────
    cudaMalloc(&scf.H, nn * sizeof(Real));
    cudaMalloc(&scf.P, nn * sizeof(Real));
    cudaMalloc(&scf.P_old, nn * sizeof(Real));
    cudaMalloc(&scf.C, nn * sizeof(Real));
    cudaMalloc(&scf.eig, n * sizeof(Real));
    cudaMalloc(&scf.work, nn * sizeof(Real));
    cudaMalloc(&scf.Hcore, nn * sizeof(Real));
    cudaMalloc(&scf.X, nn * sizeof(Real));
    cudaMalloc(&scf.F_prime, nn * sizeof(Real));
    cudaMalloc(&scf.ERI, n4 * sizeof(Real));
    cudaMalloc(&d_temp_, nn * sizeof(Real));
    cudaMalloc(&d_S_, nn * sizeof(Real));
    cudaMalloc(&d_T_, nn * sizeof(Real));
    cudaMalloc(&d_V_, nn * sizeof(Real));

    cudaMemsetAsync(scf.P, 0, nn * sizeof(Real), stream);
    cudaMemsetAsync(scf.P_old, 0, nn * sizeof(Real), stream);

    // ── 3. Launch GPU integral kernels ────────────────────────────────
    // Fused 1e: S, T, V computed simultaneously on device
    call_1e_integrals(
        d_S_, d_T_, d_V_,
        scf.d_shell_alpha, scf.d_shell_coeff, scf.d_shell_center,
        scf.d_atom_Z, scf.d_atom_pos,
        n, n_prim_, natm, stream);

    // ERI: full n^4 tensor on device
    call_eri(
        scf.ERI, scf.d_shell_alpha, scf.d_shell_coeff, scf.d_shell_center,
        n, n_prim_, stream);

    // ── 4. Hcore = T + V on device (cuBLAS) ──────────────────────────
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasSetStream(handle, stream);
    cublas_handle_ = handle;

    cudaMemcpyAsync(scf.Hcore, d_T_, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    Real one = Real(1);
    detail::axpy(handle, nn, &one, d_V_, 1, scf.Hcore, 1);

    // ── 5. X = S^{-1/2} on device (cuSOLVER) ─────────────────────────
    // Eigendecompose S: S = U * diag(lambda) * U^T
    // X = U * diag(1/sqrt(lambda)) * U^T
    cusolverDnHandle_t solver;
    cusolverDnCreate(&solver);
    cusolverDnSetStream(solver, stream);
    cusolver_handle_ = solver;

    // Copy S -> X (cuSOLVER overwrites input with eigenvectors)
    cudaMemcpyAsync(scf.X, d_S_, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);

    // Query workspace and solve
    int lwork = 0;
    Real* d_eigenvalues;
    cudaMalloc(&d_eigenvalues, n * sizeof(Real));
    int* d_info;
    cudaMalloc(&d_info, sizeof(int));

    detail::syevd_bufferSize(solver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                             n, scf.X, n, d_eigenvalues, &lwork);
    Real* d_work_solver;
    cudaMalloc(&d_work_solver, lwork * sizeof(Real));

    detail::syevd(solver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                  n, scf.X, n, d_eigenvalues, d_work_solver, lwork, d_info);
    // Now scf.X = eigenvectors (U), d_eigenvalues = eigenvalues

    // Build X = U * diag(1/sqrt(lambda)) * U^T on device
    // Step: d_temp = U * diag(1/sqrt(lambda))
    // We need to scale columns of U by 1/sqrt(lambda).
    // For small n: read eigenvalues back and use a simple kernel.
    cudaStreamSynchronize(stream);

    // Read eigenvalues to host and compute diag(1/sqrt(lambda))
    std::vector<Real> h_eig(n);
    cudaMemcpy(h_eig.data(), d_eigenvalues, n * sizeof(Real), cudaMemcpyDeviceToHost);
    std::vector<Real> h_diag(nn, 0);
    for (Int i = 0; i < n; ++i)
        h_diag[i + i * n] = Real(1) / std::sqrt(h_eig[i]);

    // Upload diag to device (into d_T_ as temp)
    cudaMemcpy(d_T_, h_diag.data(), nn * sizeof(Real), cudaMemcpyHostToDevice);

    // X_final = U * diag * U^T  →  temp = U * diag, X = temp * U^T
    Real zero = Real(0);
    detail::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                 n, n, n, &one, scf.X, n, d_T_, n, &zero, d_temp_, n);
    detail::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_T,
                 n, n, n, &one, d_temp_, n, scf.X, n, &zero, d_V_, n);

    // Copy result back to scf.X
    cudaMemcpyAsync(scf.X, d_V_, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);

    cudaFree(d_eigenvalues);
    cudaFree(d_work_solver);
    cudaFree(d_info);

    cudaStreamSynchronize(stream);
}

void RHFFockDevice::launch(const LaunchContext& ctx) {
    // Identical to RHFFock::launch — same Fock build logic
    auto& scf = ctx.frame.scf;
    if (scf.converged) return;

    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto handle = static_cast<cublasHandle_t>(cublas_handle_);
    Int n = scf.n;
    Int nn = n * n;

    cudaMemsetAsync(scf.work, 0, nn * sizeof(Real), stream);
    int threads = cuda::kDefaultBlock;
    int blocks = cuda::grid_size(nn, threads);
    build_G_kernel<<<blocks, threads, 0, stream>>>(scf.work, scf.P, scf.ERI, n);

    cudaMemcpyAsync(scf.H, scf.Hcore, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    Real one = Real(1);
    detail::axpy(handle, nn, &one, scf.work, 1, scf.H, 1);

    Real zero = Real(0);
    detail::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                 n, n, n, &one, scf.H, n, scf.X, n, &zero, d_temp_, n);
    detail::gemm(handle, CUBLAS_OP_T, CUBLAS_OP_N,
                 n, n, n, &one, scf.X, n, d_temp_, n, &zero, scf.F_prime, n);

    cudaMemcpyAsync(scf.work, scf.H, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    cudaMemcpyAsync(scf.H, scf.F_prime, nn * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
}

void RHFFockDevice::destroy(const DestroyContext& /*ctx*/) {
    if (cublas_handle_) { cublasDestroy(static_cast<cublasHandle_t>(cublas_handle_)); cublas_handle_ = nullptr; }
    if (cusolver_handle_) { cusolverDnDestroy(static_cast<cusolverDnHandle_t>(cusolver_handle_)); cusolver_handle_ = nullptr; }
    cudaFree(d_temp_); d_temp_ = nullptr;
    cudaFree(d_S_); d_S_ = nullptr;
    cudaFree(d_T_); d_T_ = nullptr;
    cudaFree(d_V_); d_V_ = nullptr;
}

} // namespace scf
} // namespace ATV

#endif // ATV_HAS_RHF
