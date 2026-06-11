/// @file cuda/minimizer/sd.cuh
/// @brief cuda::minimizer::SteepestDescent — steepest-descent Optimizer
///        component.
///
/// `SteepestDescent` is an `Optimizer`-role component that moves atoms
/// along the force (negative gradient) direction each iteration:
/// `x += alpha * f / m * kFtm2v`. It is hot-swappable: it plugs into the
/// single `cuda::minimizer::Minimize` driver exactly as the CPU
/// `cpu::minimizer::LBFGS` / `CG` Optimizer components plug into
/// `cpu::minimizer::Minimize`.
///
/// Phase sequence per iteration (driven by `cuda::minimizer::Minimize`):
///   search_direction (SD step) -> line_search -> force ->
///   convergence_check -> end_of_step
///
/// At step 0 forces are zero (from initialization), so no movement
/// occurs. From step 1 onward, atoms move downhill using forces from the
/// previous iteration.
///
/// @see cuda::minimizer::Minimize for the driver that owns the loop
/// @see minimizer::Phase for the minimizer phase sequence
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component/component.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace minimizer {

/// @brief Steepest-descent Optimizer component: `x += alpha * f / m * kFtm2v`.
///
/// Auto-binds to the `search_direction` phase. The `alpha` parameter
/// controls the step size (a pseudo-timestep in fs) — same unit
/// conversion as `VelocityVerlet`. Tagged `role()==Optimizer` so
/// `cuda::OptimizeDriver` build-time validation accepts it as the
/// required optimization step.
///
/// @see cuda::minimizer::Minimize for the driver this component plugs into
class SteepestDescent : public ComponentBase {
  public:
    /// @brief Construct a steepest-descent step.
    /// @param alpha Step size [fs] (pseudo-timestep, same units as
    ///        VelocityVerlet dt).
    explicit SteepestDescent(double alpha) : alpha_(alpha) {}

    /// @brief Component name.
    /// @return "steepest_descent".
    std::string name() const override { return "steepest_descent"; }

    /// @brief Structural role of this component.
    /// @return ComponentRole::Optimizer.
    ComponentRole role() const override { return ComponentRole::Optimizer; }

    /// @brief Default bindings: search_direction at order 0.
    /// @return Single binding to search_direction.
    std::vector<DefaultBinding> default_bindings() const override {
        return {
            {static_cast<int>(ATV::minimizer::Phase::search_direction), 0}};
    }

    /// @brief Launch the steepest-descent position-update kernel.
    /// @param ctx Launch context (frame, stream, runtime).
    void launch(const LaunchContext& ctx) override;

    /// @brief Step-size parameter.
    /// @return alpha [fs].
    double alpha() const { return alpha_; }

  private:
    double alpha_; ///< Step size (pseudo-timestep) [fs].
};

} // namespace minimizer
} // namespace cuda
} // namespace ATV
