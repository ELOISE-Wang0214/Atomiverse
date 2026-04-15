/// @file cpu/minimizer/final_geometry.h
/// @brief FinalGeometry -- single-frame XYZ snapshot of the optimizer output.
///
/// Writes the current geometry to an XYZ file each `end_of_step`, truncating
/// any previous content. After `driver.run()` returns, the file holds the
/// last geometry visited by the optimizer -- whether convergence was reached
/// or the iteration budget was exhausted.
///
/// @note Single-frame by design. For a per-step trajectory, use
///       `cpu::minimizer::TrajectoryWriter`.
///
/// Usage:
/// @code
///   opt.on(minimizer::Phase::end_of_step)
///      .add(cpu::minimizer::FinalGeometry{"aspirin_final.xyz"});
/// @endcode
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Writer that keeps a single-frame XYZ reflecting the latest geometry.
class FinalGeometry : public cpu::ComponentBase {
public:
    /// @brief Construct.
    /// @param filename Output XYZ path (overwritten each step).
    explicit FinalGeometry(std::string filename);

    std::string name() const override { return "final_geometry"; }

    std::vector<DefaultBinding> default_bindings() const override {
        // Run late in end_of_step so that converged-step forces/energy have
        // already been written to the frame by earlier components.
        return {{static_cast<int>(ATV::minimizer::Phase::end_of_step), 10}};
    }

    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Number of times the file has been rewritten.
    int n_writes() const { return n_writes_; }

private:
    std::string filename_;
    int n_writes_ = 0;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
