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
    NoseHoover(double T_target, double tau, double dt);

    std::string name() const override { return "nose_hoover"; }
    ComponentRole role() const override { return ComponentRole::Thermostat; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), +100},
                {static_cast<int>(ATV::md::Phase::final_integrate), +100}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Retune the thermostat setpoint between `Driver::run()` calls.
    ///
    /// Updates `T_target_` and rescales the chain mass `Q_` to the new
    /// setpoint using the global DOF count cached at `build()`. Chain
    /// coordinate `xi_` is preserved by default — a continuous ramp keeps
    /// the integrator state, where a hard reset (e.g. between independent
    /// equilibration runs) would zero it.
    ///
    /// Invariant: must be called between `Driver::run(...)` invocations on
    /// a built driver, never during a launch. The setter does not touch
    /// `n_dof_global_`, so the post-build atom count is reused; rebuild
    /// the driver if the system size changes.
    ///
    /// @param T_target New target temperature [K].
    /// @param reset_chain If true, also resets the integrator variable
    ///        `xi_` to 0 (use only when starting a fresh equilibration).
    void set_T_target(double T_target, bool reset_chain = false);

    /// @return Current setpoint [K] — for diagnostics / observers.
    double T_target() const noexcept { return T_target_; }

private:
    double T_target_, tau_, dt_;
    double xi_ = 0;    ///< Thermostat variable.
    double Q_ = 0;     ///< Thermostat mass.
    int    n_dof_global_ = 0;  ///< Cached at build() — used by set_T_target.
};

} // namespace md
} // namespace cpu
} // namespace ATV
