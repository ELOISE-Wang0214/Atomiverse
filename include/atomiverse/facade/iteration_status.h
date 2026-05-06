/// @file facade/iteration_status.h
/// @brief IterationStatus — outer-driver-loop control channel.
///
/// Replaces an earlier draft's `std::atomic<bool> Context::converged` with
/// a structured object that distinguishes:
///
///   - **converged** — Optimizer-style normal successful termination.
///     `OptimizeDriver` breaks the loop when this is set.
///     `DynamicsDriver` ignores this bit (time evolution doesn't "converge").
///   - **stop_requested** — explicit early termination from any source
///     (NaN detector, wall-time guard, future StopCondition components,
///     user-injected interrupts). All looping drivers honour this bit.
///   - **reason** — human-readable detail surfaced in error messages
///     and logs.
///
/// Ownership rule: only the **outer** category driver holds and resets the
/// status (per iteration). `ConvergenceCheck`-tagged components write to
/// it; SCF internals never touch it (they manage their own internal
/// convergence inside `SCFSolver`).
#pragma once

#include <string>

namespace ATV {

/// @brief Outer-driver loop status, owned by the category driver, written by
///        `ConvergenceCheck` components, polled at the end of each iteration.
struct IterationStatus {
    bool converged = false;        ///< Optimizer normal success.
    bool stop_requested = false;   ///< Explicit early termination.
    std::string reason;            ///< Human-readable detail.

    /// @brief Reset to "fresh iteration" state.
    void reset() noexcept {
        converged = false;
        stop_requested = false;
        reason.clear();
    }
};

} // namespace ATV
