/// @file cpu/single_point.h
/// @brief cpu::single_point — free-function syntactic sugar for one-shot
///        evaluator runs.
///
/// **NOT a class.** The spec's binding rules permit a helper function as
/// sugar for the `Evaluate + Evaluator` pattern, but explicitly forbid a
/// `SinglePoint` driver class — single-point energy is one outer
/// `Evaluate` task, not its own driver category.
///
/// Usage:
/// @code
///   cpu::scf::SCFEvaluator scf{cfg};
///   scf.with_gradient(false);
///   cpu::single_point(frame, scf);
///   // frame.energy populated; scf.last_result().energy in Hartree.
/// @endcode
///
/// The evaluator is borrowed (lvalue ref) so the caller retains
/// ownership and can inspect its `last_result()`-style state after the
/// call. To register multiple components, construct a `cpu::Evaluate`
/// directly and add them yourself.
#pragma once

#include "atomiverse/cpu/evaluate.h"
#include "atomiverse/cpu/frame.h"

namespace ATV {
namespace cpu {

/// @brief Run a single-pass evaluation with one evaluator component.
/// @param frame      System state, updated in place.
/// @param evaluator  Evaluator component to run (lvalue — caller retains).
/// @param phase      Phase slot to bind the evaluator on. Defaults to
///                   `evaluate::Phase::evaluate`. Override only if the
///                   evaluator needs to run during `prepare` or
///                   `finalize` (rare).
template <typename Evaluator>
void single_point(cpu::Frame& frame, Evaluator& evaluator,
                  ATV::evaluate::Phase phase = ATV::evaluate::Phase::evaluate) {
    cpu::Evaluate eval;
    eval.on(phase).add(evaluator);
    eval.build();
    eval.run(frame);
}

} // namespace cpu
} // namespace ATV
