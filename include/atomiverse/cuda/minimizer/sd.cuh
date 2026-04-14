/// @file cuda/minimizer/sd.cuh
/// @brief minimizer::SD -- steepest descent energy minimizer driver.
///
/// SteepestDescentStep moves atoms along the force (negative gradient) direction
/// each iteration: x += alpha * f / m * kFtm2v.
///
/// Phase sequence per iteration:
///   search_direction (SD step) -> line_search -> force -> convergence_check -> end_of_step
///
/// At step 0 forces are zero (from initialization), so no movement occurs.
/// From step 1 onward, atoms move downhill using forces from the previous iteration.
///
/// @see minimizer::Phase for the minimizer phase sequence
/// @see minimizer::CG for the conjugate gradient driver
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/cuda/minimizer/driver.cuh"
#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace minimizer {

/// @brief Steepest descent step: x += alpha * f / m * kFtm2v.
///
/// Auto-binds to search_direction phase. The alpha parameter controls the
/// step size (pseudo-timestep in fs). Same unit conversion as VelocityVerlet.
///
/// @see minimizer::SD for the driver that auto-adds zeroing infrastructure
class SteepestDescentStep : public ComponentBase {
public:
    /// @brief Construct steepest descent step.
    /// @param alpha Step size [fs] (pseudo-timestep, same units as VelocityVerlet dt).
    explicit SteepestDescentStep(double alpha) : alpha_(alpha) {}

    /// @brief Component name.
    /// @return "steepest_descent".
    std::string name() const override { return "steepest_descent"; }

    /// @brief Default bindings: search_direction at order 0.
    /// @return Single binding to search_direction.
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(Phase::search_direction), 0}};
    }

    /// @brief Launch SD position update kernel.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;

    /// @brief Step size parameter.
    /// @return alpha [fs].
    double alpha() const { return alpha_; }

private:
    double alpha_; ///< Step size (pseudo-timestep) [fs].
};

/// @brief Transitional steepest descent minimizer facade.
///
/// Usage:
/// @code{.cpp}
/// minimizer::SD min;
/// min.add(minimizer::SteepestDescentStep{0.01});
/// min.on(minimizer::Phase::force)
///    .add(md::NeighborList{3.0}, -10)
///    .add(md::LJ{1.0, 1.0, 3.0});
/// min.build();
/// min.run(frame, 500);
/// @endcode
///
/// @see minimizer::CG for the conjugate gradient driver
class SD : public ATV::detail::LegacyDriverFacade<Spec, Driver> {};

} // namespace minimizer
} // namespace ATV
