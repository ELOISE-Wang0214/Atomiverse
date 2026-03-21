/// @file cpu/md/nose_hoover.h
/// @brief CPU Nose-Hoover thermostat.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Single-chain Nose-Hoover thermostat (Martyna et al., 1992).
///
/// Rescales velocities to maintain target temperature.
/// Auto-binds to initial_integrate and final_integrate at order +100.
class NoseHoover : public cpu::ComponentBase {
public:
    /// @param T_target Target temperature [K].
    /// @param tau Coupling time constant [fs].
    /// @param dt Timestep [fs].
    NoseHoover(Real T_target, Real tau, Real dt);

    std::string name() const override { return "nose_hoover"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), +100},
                {static_cast<int>(ATV::md::Phase::final_integrate), +100}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

private:
    Real T_target_, tau_, dt_;
    Real xi_ = 0;    ///< Thermostat variable.
    Real Q_ = 0;     ///< Thermostat mass.
};

} // namespace md
} // namespace cpu
} // namespace ATV
