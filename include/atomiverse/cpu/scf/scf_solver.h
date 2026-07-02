/// @file cpu/scf/scf_solver.h
/// @brief cpu::scf::SCFSolver — public low-level inner SCF solver.
///
/// **NOT a `Driver<...>` subclass** — `SCFSolver` is the *inner* solver
/// (the SCF self-consistency loop) that lives behind force-providing
/// evaluator components. The user-facing way to compute single-point
/// electronic structure is `cpu::Evaluate + cpu::scf::SCFEvaluator`;
/// `SCFSolver` is the internal machinery `SCFEvaluator` calls each step.
///
/// `SCFSolver` is also a public API for low-level use (debugging,
/// testing, density-matrix dumps, DIIS comparisons, initial-guess
/// experiments). It owns five `internal::*` step helpers, each with
/// build/step lifecycle and pre-allocated workspaces. The helpers are
/// **not** registered into any outer `Driver` Plan.
///
/// Usage (low-level direct SCF):
/// @code
///   cpu::scf::SCFSolver solver{cpu::scf::default_config(integrals)};
///   solver.build(integrals);
///   const auto& result = solver.run(integrals);
///   // result.energy, result.P_alpha, result.eigenvalues_alpha, ...
/// @endcode
///
/// Usage (single-point through the outer facade):
/// @code
///   cpu::Evaluate eval;
///   eval.add(cpu::scf::SCFEvaluator{config});  // owns an SCFSolver internally
///   eval.build();
///   eval.run(frame);
/// @endcode
#pragma once

#include <memory>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h" // SCFConfig, SCFInitialGuess, SCFResult
#include "atomiverse/cpu/scf/dense_matrix.h" // RealMatrix

namespace ATV {
namespace cpu {
namespace scf {

// Facade-private step helpers — forward-declared only. The definitions
// live in include/atomiverse/cpu/scf/internal/ and are included from
// scf_solver.cpp; this public header must not include internal/*.h
// (enforced by test_cpu_scf_solver_internal_isolation).
namespace internal {
class HamiltonianBuilder;
class EigenSolver;
class DensityBuilder;
class DIISMixer;
class DualConvergenceCheck;
} // namespace internal

/// @brief Inner SCF self-consistency loop. Public low-level surface.
///
/// Construction → build(reference) → run(integrals[, guess]).
/// `last_result()` is valid after each `run()` call.
class SCFSolver {
  public:
    /// @brief Construct from an SCF configuration.
    /// @param config  Method, spin, convergence, grid parameters.
    explicit SCFSolver(SCFConfig config);

    ~SCFSolver();
    SCFSolver(SCFSolver&&) noexcept;
    SCFSolver& operator=(SCFSolver&&) noexcept;
    SCFSolver(const SCFSolver&) = delete;
    SCFSolver& operator=(const SCFSolver&) = delete;

    /// @brief One-time setup: pre-allocates the n_basis-sized eigen /
    /// density workspaces. Geometry-dependent structures (X = S^{-1/2},
    /// shell-pair lists, KS grid) are deliberately rebuilt inside each
    /// `run()` — atoms move during AIMD / geometry optimization, so they
    /// must track the current positions. Subsequent `run()` calls assume
    /// the same n_basis (fixed by the basis set).
    /// @param reference IntegralPack used to size internal workspaces.
    void build(const IntegralPack& reference);

    /// @brief Execute the SCF cycle on `integrals`.
    /// @param integrals  Molecular integrals at the current geometry.
    /// @param guess      Optional density-matrix initial guess.
    /// @return Reference to the owned `SCFResult` (valid until the next
    ///         `run()` or `SCFSolver` destruction).
    SCFResult& run(const IntegralPack& integrals,
                   const SCFInitialGuess* guess = nullptr);

    /// @brief Access the most recent SCF result.
    const SCFResult& last_result() const { return result_; }

    /// @brief Access the current SCF configuration (for inspection).
    const SCFConfig& config() const { return config_; }

  private:
    SCFConfig config_;
    bool built_ = false;

    // Per-run state (rebuilt each run(); geometry-dependent).
    RealMatrix X_; ///< Orthogonalization matrix. Löwdin S^{-1/2} (n × n)
                   ///< when the basis is healthy or lindep filtering is
                   ///< disabled; rectangular canonical X (n × m, m < n)
                   ///< when overlap eigen-directions were dropped. Owned.
    int n_dropped_ = 0; ///< Overlap eigen-directions dropped this run.

    // Internal SCF step helpers. Held by pointer so this public header
    // only needs the forward declarations above (the internal/*.h
    // definitions stay facade-private).
    std::unique_ptr<internal::HamiltonianBuilder> hamiltonian_;
    std::unique_ptr<internal::EigenSolver> eigen_;
    std::unique_ptr<internal::DensityBuilder> density_;
    std::unique_ptr<internal::DIISMixer> diis_;
    std::unique_ptr<internal::DualConvergenceCheck> convergence_;

    // Owned result buffer (returned by reference from run()).
    SCFResult result_;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
