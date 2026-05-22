/// @file cuda/scf/scf_solver.cuh
/// @brief cuda::scf::SCFSolver — public low-level inner SCF solver.
///
/// **NOT a `Driver<...>` subclass** — `SCFSolver` is the *inner* solver
/// (the SCF self-consistency loop) that lives behind force-providing
/// evaluator components. The user-facing way to compute single-point
/// CUDA electronic structure is `cuda::Evaluate + cuda::scf::SCFEvaluator`;
/// `SCFSolver` is the internal machinery `SCFEvaluator` calls.
///
/// `SCFSolver` is also a public API for low-level use (debugging, testing,
/// CPU/CUDA parity). It owns its `internal::*` step helpers, each with a
/// `build` / `step` / `destroy` lifecycle and pre-allocated device
/// workspaces. The helpers are **not** registered into any `Driver` Plan
/// and do not reference `scf::Phase` or `cuda::Context::outer_status`.
///
/// Three run paths are surfaced as distinct constructors:
///  - **Synthetic** — `H(P) = H0 + alpha*P` toy self-consistency cycle.
///  - **RHF (host integrals)** — Fock build from a precomputed
///    `cpu::scf::IntegralPack` (`ATV_HAS_RHF` only).
///  - **RHF (device integrals)** — Fock build from GPU-native integral
///    kernels, zero CPU integral computation (`ATV_HAS_RHF` only).
///
/// Usage (low-level direct SCF):
/// @code
///   cuda::scf::SCFSolver solver{n_occ, alpha, tol, max_iter};
///   solver.build(runtime, frame);
///   solver.run(runtime, frame);
///   // frame.scf.converged / frame.scf.residual / frame.energy populated.
/// @endcode
///
/// Usage (single-point through the outer facade):
/// @code
///   cuda::Evaluate eval;
///   eval.add(cuda::scf::SCFEvaluator{...});  // owns an SCFSolver internally
///   eval.build();
///   eval.run(frame);
/// @endcode
#pragma once

#include <memory>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/cuda/scf/internal/convergence_check.cuh"
#include "atomiverse/cuda/scf/internal/density_builder.cuh"
#include "atomiverse/cuda/scf/internal/eigen_solver.cuh"
#include "atomiverse/cuda/scf/internal/hamiltonian.cuh"

#ifdef ATV_HAS_RHF
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cuda/scf/internal/rhf_density.cuh"
#include "atomiverse/cuda/scf/internal/rhf_fock.cuh"
#endif

namespace ATV {
namespace cuda {
namespace scf {

/// @brief Inner CUDA SCF self-consistency loop. Public low-level surface.
///
/// Construction → `build(runtime, frame)` (one-time device-workspace
/// allocation) → `run(runtime, frame)` which loops the SCF cycle
/// internally to convergence or a max-iteration cap. The density matrix
/// is zeroed internally as the first step of `run()`.
class SCFSolver {
  public:
    /// @brief Construct a synthetic-path solver `H(P) = H0 + alpha*P`.
    /// @param n_occ Number of occupied orbitals [dimensionless].
    /// @param alpha Density coupling parameter [dimensionless].
    /// @param tol Convergence tolerance, Frobenius norm of P change
    ///        [dimensionless].
    /// @param max_iter Maximum number of SCF iterations [dimensionless].
    SCFSolver(int n_occ, double alpha, double tol, int max_iter);

#ifdef ATV_HAS_RHF
    /// @brief Construct an RHF solver driven by host-computed integrals.
    /// @param integrals Precomputed molecular integral pack
    ///        (Hcore, ERI, S, V_nuc).
    /// @param tol Convergence tolerance, Frobenius norm of P change
    ///        [dimensionless].
    /// @param max_iter Maximum number of SCF iterations [dimensionless].
    SCFSolver(const cpu::scf::IntegralPack& integrals, double tol,
              int max_iter);

    /// @brief Construct an RHF solver driven by GPU-native integral kernels.
    /// @param n_prim Primitives per shell (3 for STO-3G) [dimensionless].
    /// @param exponents Gaussian exponents [n_prim] [Bohr^-2].
    /// @param coefficients Raw contraction coefficients [n_prim]
    ///        (NOT normalized).
    /// @param n_occ Number of doubly-occupied orbitals [dimensionless].
    /// @param tol Convergence tolerance, Frobenius norm of P change
    ///        [dimensionless].
    /// @param max_iter Maximum number of SCF iterations [dimensionless].
    SCFSolver(int n_prim, const double* exponents, const double* coefficients,
              int n_occ, double tol, int max_iter);
#endif

    /// @brief One-time setup: allocate device workspaces and build the
    ///        owned `internal::*` step helpers.
    /// @param runtime CUDA execution context (streams).
    /// @param frame Device-resident system state, sized at this call.
    void build(cuda::Context& runtime, cuda::Frame& frame);

    /// @brief Execute the SCF self-consistency cycle on `frame`.
    ///
    /// Zeroes the density matrix internally, then loops the fixed-order
    /// SCF cycle (Fock/Hamiltonian build → diagonalize → density build →
    /// convergence check) up to `max_iter` iterations, breaking early when
    /// `frame.scf.converged` is set. Performs no device allocation.
    /// @param runtime CUDA execution context (streams).
    /// @param frame Device-resident system state, updated in place.
    void run(cuda::Context& runtime, cuda::Frame& frame);

    /// @brief Release the cuBLAS/cuSOLVER handles and per-helper scratch
    ///        buffers held by the owned `internal::*` step helpers.
    ///
    /// The shared SCF workspace in `frame.scf` is freed by the backend's
    /// `cuda::destroy(frame)`; this releases only the helper-private
    /// resources allocated in `build()`.
    /// @param runtime CUDA execution context (streams).
    /// @param frame Device-resident system state.
    void destroy(cuda::Context& runtime, cuda::Frame& frame);

    /// @brief Number of SCF iterations executed by the most recent `run()`.
    /// @return Iteration count [dimensionless].
    int iterations() const { return iterations_; }

  private:
    /// @brief Run-path selector — fixes which `internal::*` helpers drive
    ///        the Fock/Hamiltonian build and density rebuild.
    enum class Path {
        Synthetic,          ///< H(P) = H0 + alpha*P toy cycle.
        RHFHostIntegrals,   ///< RHF Fock build from a host IntegralPack.
        RHFDeviceIntegrals, ///< RHF Fock build from GPU integral kernels.
    };

    Path path_;     ///< Selected run path.
    double tol_;    ///< Convergence tolerance (Frobenius norm).
    int max_iter_;  ///< Maximum SCF iteration count.
    bool built_ = false;     ///< Whether build() has run.
    int iterations_ = 0;     ///< Iterations executed by the last run().

    // Step helpers shared by every path.
    std::unique_ptr<internal::EigenSolver> eigen_;
    std::unique_ptr<internal::ConvergenceCheck> convergence_;

    // Synthetic-path helpers.
    std::unique_ptr<internal::Hamiltonian> hamiltonian_;
    std::unique_ptr<internal::DensityBuilder> density_;

#ifdef ATV_HAS_RHF
    // RHF-path helpers (one Fock builder is active per solver).
    std::unique_ptr<internal::RHFFock> rhf_fock_;
    std::unique_ptr<internal::RHFFockDevice> rhf_fock_device_;
    std::unique_ptr<internal::RHFDensity> rhf_density_;
#endif
};

} // namespace scf
} // namespace cuda
} // namespace ATV
