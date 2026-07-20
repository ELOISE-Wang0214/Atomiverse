/// @file cpu/component/criterion.h
/// @brief cpu::Criterion — termination as a Component kind.
///
/// A `Criterion` is a Component that also answers ONE question after every
/// iteration of the loop it is registered on: *should the loop stop?* It returns
/// a `Verdict` — `Continue`, `Converged` or `Stop` — and nothing else. Like a
/// `LoopHandler` it only ever RETURNS a verdict; it has no field with which to
/// write the loop's status, so a component still cannot declare its own loop
/// converged. The driver PULLS every registered criterion's `verdict()`, composes
/// them by severity (Stop > Converged > Continue), all-reduces once across MPI
/// ranks, and — the driver alone — writes `converged` / `stop_requested`.
///
/// Being a Component (not the free-standing `LoopHandler` concept) is what lets a
/// criterion also bind a phase and MEASURE: e.g. a later SCF `ScfConverged` binds
/// `scf::Phase::check` to read the energy Δ / DIIS residual it then judges. A
/// verdict-only criterion declares no `default_bindings()` and simply enters the
/// driver's polling list.
///
/// This is additive: the existing `LoopHandler` / `add_handler` API is untouched
/// in this sub-spec; the tree-wide migration onto `Criterion` comes later.
///
/// One kind per header (Concrete → Kind → Component). `loop_handler.h` supplies
/// `LoopState` and does not include `component.h`, so including both here is
/// cycle-free.
#pragma once

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/driver/loop_handler.h"
#include "atomiverse/facade/loop.h"

namespace ATV {
namespace cpu {

/// @brief A termination criterion: a Component the driver polls once per iteration.
class Criterion : public Component {
  public:
    /// @brief Judge the iteration just completed.
    ///
    /// Mirrors `LoopHandler::after_iteration`: non-const (a criterion may update
    /// its own measurement state), takes the loop's read-only `LoopState`. Only
    /// RETURNS a verdict — it cannot write the loop's status; the driver composes
    /// the returned verdicts and writes convergence itself.
    /// @param state Iteration index, frame, context, metrics.
    /// @return `Continue`, `Converged` or `Stop`. @see ATV::Verdict
    virtual ATV::Verdict verdict(const cpu::LoopState& state) = 0;
};

} // namespace cpu
} // namespace ATV
