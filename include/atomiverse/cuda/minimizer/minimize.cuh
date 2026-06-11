/// @file cuda/minimizer/minimize.cuh
/// @brief cuda::minimizer::Minimize — concrete CUDA energy-minimizer driver.
///
/// `Minimize` is the single user-facing CUDA minimizer. It is a
/// `cuda::OptimizeDriver` — the convergence-aware outer category — and so
/// inherits that base's `run(frame, max_iter)` loop and its build-time
/// `validate_components()` contract (≥ 1 `Optimizer`-role component and
/// ≥ 1 `Evaluator` carrying `Capability::Energy` or `Capability::Force`).
///
/// Algorithms are hot-swappable components, not driver subclasses:
/// `cuda::minimizer::SteepestDescent`, `cuda::minimizer::CG` and
/// `cuda::minimizer::LBFGS` are Optimizer components plugged in via
/// `add()`. This mirrors `cpu::minimizer::Minimize`, which accepts
/// `cpu::minimizer::LBFGS` / `CG` the same way.
///
/// The constructor adds only infrastructure components — it contributes
/// no runtime behaviour beyond the base. The phase-agnostic `Zero*`
/// initializers (`cuda/md/init.cuh`) are bound explicitly to
/// `minimizer::Phase` via `on(phase).add()`.
///
/// Early stop: `CG` and `LBFGS` bind the `convergence_check` phase and
/// set `outer_status.converged` under the shared five-criteria rule
/// (`cuda::minimizer::ConvCheck`), so `run()` returns as soon as the
/// system is converged. `SteepestDescent` carries no convergence check
/// and runs the full `max_iter` budget.
///
/// @see cuda::OptimizeDriver for the inherited loop + validation
/// @see cuda::minimizer::SteepestDescent for the steepest-descent Optimizer
/// @see cuda::minimizer::CG / cuda::minimizer::LBFGS for converging Optimizers
/// @see minimizer::Phase for the minimizer phase sequence
#pragma once

#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/optimize_driver.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace minimizer {

/// @brief CUDA energy-minimization driver — a `cuda::OptimizeDriver`.
///
/// Binds the phase-agnostic `Zero*` infrastructure components to
/// `minimizer::Phase::force`:
/// - ZeroForces at force phase (order -1000) — zeros fx, fy, fz.
/// - ZeroEnergy at force phase (order -999) — zeros potential energy.
///
/// Validation rules inherited from `cuda::OptimizeDriver`: requires ≥ 1
/// `Optimizer` component (e.g. `SteepestDescent`) and ≥ 1 `Evaluator`
/// carrying `Capability::Energy` or `Capability::Force` (e.g. `LJ`).
/// Users register the Optimizer component and force evaluators through
/// the inherited public `add()` / `on(phase).add()`.
///
/// Usage:
/// @code{.cpp}
/// cuda::minimizer::Minimize min;
/// min.add(cuda::minimizer::SteepestDescent{0.01});
/// min.on(minimizer::Phase::force)
///    .add(cuda::NeighborList{3.0}, -10)
///    .add(cuda::PairLJCut{1.0, 1.0, 3.0});
/// min.build();          // validate_components() then compile the Plan
/// min.run(frame, 500);
/// @endcode
///
/// Inheritance: Minimize → cuda::OptimizeDriver →
/// cuda::CUDADriver<minimizer::Phase> →
/// Driver<minimizer::Phase, cuda::device_tag>.
class Minimize : public OptimizeDriver {
  public:
    /// @brief Construct the minimizer, binding the phase-agnostic `Zero*`
    ///        infrastructure components to `minimizer::Phase::force`.
    Minimize() {
        on(ATV::minimizer::Phase::force).add(md::ZeroForces{}, -1000);
        on(ATV::minimizer::Phase::force).add(md::ZeroEnergy{}, -999);
    }
};

} // namespace minimizer
} // namespace cuda
} // namespace ATV
