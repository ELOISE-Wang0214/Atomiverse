/// @file cuda/scf/scf_solver.cu
/// @brief cuda::scf::SCFSolver implementation — the inner SCF
///        self-consistency loop.
///
/// `SCFSolver` orchestrates its owned `internal::*` step helpers through
/// the fixed-order SCF cycle. It absorbs the former `ZeroDensity`
/// infrastructure component: the density matrix is zeroed via
/// `cudaMemsetAsync` as the first step of `run()`. The solver is not a
/// `Driver` and carries no phase enum, and never touches
/// `cuda::Context::outer_status` — it manages its own convergence
/// internally via `frame.scf.converged`.

#include "atomiverse/cuda/scf/scf_solver.cuh"

#include <cuda_runtime.h>

#include "atomiverse/cuda/types.cuh"

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
    // and must build first; the remaining helpers size off frame.scf.n.
    switch (path_) {
    case Path::Synthetic:
        hamiltonian_->build(ctx);
        eigen_->build(ctx);
        density_->build(ctx);
        convergence_->build(ctx);
        break;
#ifdef ATV_HAS_RHF
    case Path::RHFHostIntegrals:
        rhf_fock_->build(ctx);
        eigen_->build(ctx);
        rhf_density_->build(ctx);
        convergence_->build(ctx);
        break;
    case Path::RHFDeviceIntegrals:
        rhf_fock_device_->build(ctx);
        eigen_->build(ctx);
        rhf_density_->build(ctx);
        convergence_->build(ctx);
        break;
#endif
    }
    built_ = true;
}

void SCFSolver::run(cuda::Context& runtime, cuda::Frame& frame) {
    if (!built_)
        build(runtime, frame);

    auto stream = static_cast<cudaStream_t>(runtime.stream);
    auto& scf = frame.scf;

    // ── Absorbed ZeroDensity: zero P before the first iteration ───────
    scf.converged = false;
    scf.iteration = 0;
    if (scf.n > 0)
        cudaMemsetAsync(scf.P, 0, scf.n * scf.n * sizeof(cuda::prec_t),
                        stream);

    // ── Internal SCF self-consistency loop ────────────────────────────
    // Fixed-order cycle: build Fock/Hamiltonian → diagonalize →
    // build density → convergence check. Helpers no-op once
    // frame.scf.converged is set; the loop also breaks early on it.
    iterations_ = 0;
    for (int iter = 0; iter < max_iter_; ++iter) {
        LaunchContext ctx{runtime, frame, /*stage=*/0, /*step=*/iter};

        switch (path_) {
        case Path::Synthetic:
            hamiltonian_->step(ctx);
            eigen_->step(ctx);
            density_->step(ctx);
            convergence_->step(ctx);
            break;
#ifdef ATV_HAS_RHF
        case Path::RHFHostIntegrals:
            rhf_fock_->step(ctx);
            eigen_->step(ctx);
            rhf_density_->step(ctx);
            convergence_->step(ctx);
            break;
        case Path::RHFDeviceIntegrals:
            rhf_fock_device_->step(ctx);
            eigen_->step(ctx);
            rhf_density_->step(ctx);
            convergence_->step(ctx);
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
        convergence_->destroy(ctx);
        density_->destroy(ctx);
        eigen_->destroy(ctx);
        hamiltonian_->destroy(ctx);
        break;
#ifdef ATV_HAS_RHF
    case Path::RHFHostIntegrals:
        convergence_->destroy(ctx);
        rhf_density_->destroy(ctx);
        eigen_->destroy(ctx);
        rhf_fock_->destroy(ctx);
        break;
    case Path::RHFDeviceIntegrals:
        convergence_->destroy(ctx);
        rhf_density_->destroy(ctx);
        eigen_->destroy(ctx);
        rhf_fock_device_->destroy(ctx);
        break;
#endif
    }
    built_ = false;
}

} // namespace scf
} // namespace cuda
} // namespace ATV
