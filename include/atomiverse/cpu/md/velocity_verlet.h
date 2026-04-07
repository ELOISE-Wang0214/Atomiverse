/// @file cpu/md/velocity_verlet.h
/// @brief CPU VelocityVerlet integrator.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Velocity Verlet symplectic integrator (Swope et al., 1982).
///
/// initial_integrate: v(t+dt/2) = v(t) + (dt/2)*f(t)/m, x(t+dt) = x(t) + dt*v(t+dt/2)
/// final_integrate:   v(t+dt)   = v(t+dt/2) + (dt/2)*f(t+dt)/m
class VelocityVerlet : public cpu::ComponentBase {
public:
    /// @param dt Timestep [fs].
    explicit VelocityVerlet(double dt);

    std::string name() const override { return "velocity_verlet"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), 0},
                {static_cast<int>(ATV::md::Phase::final_integrate), 0}};
    }
    void launch(const cpu::LaunchContext& ctx) override;

private:
    double dt_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
