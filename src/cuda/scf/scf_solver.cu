/// @file cuda/scf/scf_solver.cu
/// @brief cuda::scf::SCFSolver implementation — the inner SCF
///        self-consistency loop.
///
/// `SCFSolver` orchestrates its owned `internal::*` step helpers through
/// the fixed-order SCF cycle. It owns the `SCFWorkspace` (no longer part
/// of `cuda::Frame`) and passes it by reference into every helper call.
/// It absorbs the former `ZeroDensity` infrastructure component: the
/// density matrix is zeroed via `cudaMemsetAsync` as the first step of
/// `run()`. The solver is not a `Driver` and carries no phase enum, and
/// never touches `cuda::Context::outer_status` — it manages its own
/// convergence internally via the workspace's `converged` flag.

#include "atomiverse/cuda/scf/scf_solver.cuh"

#include <cuda_runtime.h>

#include "atomiverse/cuda/frame/types.cuh"

namespace ATV {
namespace cuda {
namespace scf {

SCFSolver::SCFSolver(int n_occ, double alpha, double tol, int max_iter)
    : path_(Path::Synthetic), tol_(tol), max_iter_(max_iter) {
    eigen_ = std::make_unique<internal::EigenSolver>();
    // Synthetic path: ConvergenceCheck owns the energy (E = trace(H*P)).
    convergence_ = std::make_unique<internal::ConvergenceCheck>(tol, true);
    hamiltonian_ = std::make_unique<internal::Hamiltonian>(n_occ, alpha);
    density_ = std::make_unique<internal::DensityBuilder>();
}

#ifdef ATV_HAS_RHF

SCFSolver::SCFSolver(const cpu::scf::IntegralPack& integrals, double tol,
                     int max_iter)
    : path_(Path::RHFHostIntegrals), tol_(tol), max_iter_(max_iter) {
    eigen_ = std::make_unique<internal::EigenSolver>();
    // RHF path: RHFDensity owns the energy; the checker is residual-only.
    convergence_ = std::make_unique<internal::ConvergenceCheck>(tol, false);
    rhf_fock_ = std::make_unique<internal::RHFFock>(integrals);
    rhf_density_ = std::make_unique<internal::RHFDensity>(integrals.n_occ());
}

SCFSolver::SCFSolver(int n_prim, const double* exponents,
                     const double* coefficients, int n_occ, double tol,
                     int max_iter)
    : path_(Path::RHFDeviceIntegrals), tol_(tol), max_iter_(max_iter) {
    eigen_ = std::make_unique<internal::EigenSolver>();
    convergence_ = std::make_unique<internal::ConvergenceCheck>(tol, false);
    rhf_fock_device_ = std::make_unique<internal::RHFFockDevice>(
        n_prim, exponents, coefficients, n_occ);
    rhf_density_ = std::make_unique<internal::RHFDensity>(n_occ);
}

#endif // ATV_HAS_RHF

void SCFSolver::build(cuda::Context& runtime, cuda::Frame& frame) {
    BuildContext ctx{runtime, frame};

    // The Fock/Hamiltonian builder allocates the SCF workspace (H, P, ...)
    // and must build first; the remaining helpers size off workspace_.n.
    switch (path_) {
    case Path::Synthetic:
        hamiltonian_->build(ctx, workspace_);
        eigen_->build(ctx, workspace_);
        density_->build(ctx, workspace_);
        convergence_->build(ctx, workspace_);
        break;
#ifdef ATV_HAS_RHF
    case Path::RHFHostIntegrals:
        rhf_fock_->build(ctx, workspace_);
        eigen_->build(ctx, workspace_);
        rhf_density_->build(ctx, workspace_);
        convergence_->build(ctx, workspace_);
        break;
    case Path::RHFDeviceIntegrals:
        rhf_fock_device_->build(ctx, workspace_);
        eigen_->build(ctx, workspace_);
        rhf_density_->build(ctx, workspace_);
        convergence_->build(ctx, workspace_);
        break;
#endif
    }
    built_ = true;
}

void SCFSolver::run(cuda::Context& runtime, cuda::Frame& frame) {
    if (!built_)
        build(runtime, frame);

    auto stream = static_cast<cudaStream_t>(runtime.stream);
    auto& scf = workspace_;

    // ── Absorbed ZeroDensity: zero P before the first iteration ───────
    scf.converged = false;
    scf.iteration = 0;
    if (scf.n > 0)
        cudaMemsetAsync(scf.P, 0, scf.n * scf.n * sizeof(cuda::prec_t),
                        stream);

    // ── Internal SCF self-consistency loop ────────────────────────────
    // Fixed-order cycle: build Fock/Hamiltonian → diagonalize →
    // build density → convergence check. Helpers no-op once
    // workspace_.converged is set; the loop also breaks early on it.
    iterations_ = 0;
    for (int iter = 0; iter < max_iter_; ++iter) {
        LaunchContext ctx{runtime, frame, /*stage=*/0, /*step=*/iter};

        switch (path_) {
        case Path::Synthetic:
            hamiltonian_->step(ctx, workspace_);
            eigen_->step(ctx, workspace_);
            density_->step(ctx, workspace_);
            convergence_->step(ctx, workspace_);
            break;
#ifdef ATV_HAS_RHF
        case Path::RHFHostIntegrals:
            rhf_fock_->step(ctx, workspace_);
            eigen_->step(ctx, workspace_);
            rhf_density_->step(ctx, workspace_);
            convergence_->step(ctx, workspace_);
            break;
        case Path::RHFDeviceIntegrals:
            rhf_fock_device_->step(ctx, workspace_);
            eigen_->step(ctx, workspace_);
            rhf_density_->step(ctx, workspace_);
            convergence_->step(ctx, workspace_);
            break;
#endif
        }

        ++iterations_;
        if (scf.converged)
            break;
    }
}

void SCFSolver::destroy(cuda::Context& runtime, cuda::Frame& frame) {
    DestroyContext ctx{runtime, frame};

    switch (path_) {
    case Path::Synthetic:
        convergence_->destroy(ctx, workspace_);
        density_->destroy(ctx, workspace_);
        eigen_->destroy(ctx, workspace_);
        hamiltonian_->destroy(ctx, workspace_);
        break;
#ifdef ATV_HAS_RHF
    case Path::RHFHostIntegrals:
        convergence_->destroy(ctx, workspace_);
        rhf_density_->destroy(ctx, workspace_);
        eigen_->destroy(ctx, workspace_);
        rhf_fock_->destroy(ctx, workspace_);
        break;
    case Path::RHFDeviceIntegrals:
        convergence_->destroy(ctx, workspace_);
        rhf_density_->destroy(ctx, workspace_);
        eigen_->destroy(ctx, workspace_);
        rhf_fock_device_->destroy(ctx, workspace_);
        break;
#endif
    }

    // Free the device-resident SCF workspace arrays owned by this solver.
    // (Formerly freed by cuda::destroy(frame) when the workspace lived on
    // the frame; the helpers above release only their private scratch.)
    cudaFree(workspace_.H0);
    cudaFree(workspace_.H);
    cudaFree(workspace_.P);
    cudaFree(workspace_.P_old);
    cudaFree(workspace_.C);
    cudaFree(workspace_.eig);
    cudaFree(workspace_.work);
    cudaFree(workspace_.Hcore);
    cudaFree(workspace_.X);
    cudaFree(workspace_.F_prime);
    cudaFree(workspace_.ERI);
    cudaFree(workspace_.d_shell_alpha);
    cudaFree(workspace_.d_shell_coeff);
    cudaFree(workspace_.d_shell_center);
    cudaFree(workspace_.d_atom_Z);
    cudaFree(workspace_.d_atom_pos);
    workspace_ = {};

    built_ = false;
}

} // namespace scf
} // namespace cuda
} // namespace ATV
