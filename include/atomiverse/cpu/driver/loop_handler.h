/// @file cpu/driver/loop_handler.h
/// @brief LoopHandler — termination as a reusable, user-definable policy.
///
/// ## Why this is not a template
///
/// `LoopHandler` sees a `LoopState` and **never a phase**. That is the whole point: a
/// handler that knew about `md::Phase` could only ever serve MD. Knowing nothing about
/// stages, one handler serves every workflow — `StopOnNaN` is as meaningful to a band of
/// NEB images as to an MD trajectory, and neither of them had any way to say so before.
///
/// Type erasure in this architecture happens at the **node** boundary (a driver is a
/// node; a component is a node) and at the **policy** boundary (here). It never happens
/// at the stage boundary — stages stay typed enums, all the way down.
///
/// ## Writing one
///
/// @code
///   class StopWhenTooHot : public cpu::LoopHandler {
///     public:
///       explicit StopWhenTooHot(double t_max) : t_max_(t_max) {}
///       const char* name() const override { return "stop_when_too_hot"; }
///       Verdict after_iteration(const cpu::LoopState& s) override {
///           return (s.frame.scalars.temperature > t_max_) ? Verdict::Stop
///                                                         : Verdict::Continue;
///       }
///     private:
///       double t_max_;
///   };
///
///   sim.add_handler(StopWhenTooHot{2000.0});
/// @endcode
///
/// This is LAMMPS's `fix halt`, and it was **impossible** to write before: "stop or not"
/// was a private `if` inside a driver class.
///
/// ## MPI
///
/// **A handler's verdict is reduced across ranks by the driver, not by the handler.**
/// Ranks that disagree about when to stop deadlock at the next collective, and that is
/// a bug no handler author should have to know about. Write `after_iteration` as if the
/// process were alone; the driver all-reduces the severity.
#pragma once

#include <memory>
#include <string>

#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/facade/loop.h"

namespace ATV {
namespace cpu {

/// @brief Everything a handler is allowed to look at.
struct LoopState {
    int iteration;               ///< 0-based index of the iteration just completed.
    const cpu::Frame& frame;     ///< The system state as it now stands.
    const cpu::Context& runtime; ///< Execution context (threads, MPI, service slices).
    const LoopMetrics& metrics;  ///< What this loop's components measured. @see facade/loop.h
};

/// @brief A termination policy, asked once after every iteration of the loop it is on.
class LoopHandler {
  public:
    virtual ~LoopHandler() = default;

    /// @brief Identifier, used in the stop reason the driver reports.
    virtual const char* name() const = 0;

    /// @brief Judge the iteration just completed.
    /// @param state Iteration index, frame, context, metrics.
    /// @return `Continue`, `Converged` or `Stop`. @see ATV::Verdict
    virtual Verdict after_iteration(const LoopState& state) = 0;
};

// ─── Built-in handlers ──────────────────────────────────────────────────────
//
// Each of these used to BE a driver category. They are now ~10 lines apiece, and they
// compose: a Minimize run that must also give up after an hour is
// `ForceConverged{} + StopOnWallClock{3600}`, which needed a new class before.

/// @brief The five-criteria geometry-optimization convergence test.
///
/// ALL five must hold. It reads `LoopMetrics`, which the Optimizer component publishes —
/// so the tolerances live here, with the *policy*, and not (as they did) as five of the
/// eight positional constructor arguments of `LBFGS`, where changing a tolerance meant
/// touching the algorithm and where `CG` had to duplicate the identical test.
class ForceConverged : public LoopHandler {
  public:
    /// @param ftol Max force [eV/Å].
    /// @param rms_ftol RMS force [eV/Å].
    /// @param max_dtol Max displacement [Å].
    /// @param rms_dtol RMS displacement [Å].
    /// @param etol Energy change [eV].
    explicit ForceConverged(double ftol = 0.05, double rms_ftol = 0.03,
                            double max_dtol = 0.002, double rms_dtol = 0.001,
                            double etol = 1e-5)
        : ftol_(ftol), rms_ftol_(rms_ftol), max_dtol_(max_dtol), rms_dtol_(rms_dtol),
          etol_(etol) {}

    const char* name() const override { return "force_converged"; }

    Verdict after_iteration(const LoopState& s) override {
        const LoopMetrics& m = s.metrics;
        // Judge the METRICS — never the iteration counter. Whether a displacement is
        // meaningful is a fact only the component knows (it has the previous geometry),
        // and it says so in `has_displacement`. A handler that second-guesses this with
        // `iteration == 0` gets it wrong, and I know because I wrote that first: in the
        // minimizer phase order the optimizer has ALREADY taken a step by the time
        // `convergence_check` runs, so driver-iteration 0 is not the optimizer's first
        // step. Skipping the displacement test there declared an untouched lattice
        // converged, having moved nothing.
        const bool disp_ok = !m.has_displacement
                             || (m.max_displacement < max_dtol_
                                 && m.rms_displacement < rms_dtol_);
        const bool ok = m.max_force < ftol_ && m.rms_force < rms_ftol_ && disp_ok
                        && m.energy_change < etol_;
        (void)s.iteration;
        return ok ? Verdict::Converged : Verdict::Continue;
    }

  private:
    double ftol_, rms_ftol_, max_dtol_, rms_dtol_, etol_;
};

/// @brief A coupled chain is converged when EVERY replica is — never when the first is.
///
/// The AND that `PathDriver` used to spell out inline. Reads
/// `LoopMetrics::replicas_converged`, which the band optimizer publishes as a fraction.
class AllReplicasConverged : public LoopHandler {
  public:
    const char* name() const override { return "all_replicas_converged"; }
    Verdict after_iteration(const LoopState& s) override {
        return (s.metrics.replicas_converged >= 1.0) ? Verdict::Converged
                                                     : Verdict::Continue;
    }
};

/// @brief Abandon a run whose energy has gone non-finite.
///
/// `Stop`, never `Converged` — which is exactly why the severity order matters. A blown-up
/// system has zero force on every atom (they are all NaN, and NaN fails every `<`
/// comparison... but a *diverged* one can genuinely park at a force minimum of garbage).
/// Reporting success for it would be worse than crashing.
class StopOnNaN : public LoopHandler {
  public:
    const char* name() const override { return "stop_on_nan"; }
    Verdict after_iteration(const LoopState& s) override {
        const double e = s.metrics.energy;
        return (e != e || e > 1e300 || e < -1e300) ? Verdict::Stop : Verdict::Continue;
    }
};

} // namespace cpu
} // namespace ATV
