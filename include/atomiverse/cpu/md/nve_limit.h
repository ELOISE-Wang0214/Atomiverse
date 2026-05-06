/// @file cpu/md/nve_limit.h
/// @brief CPU velocity-Verlet integrator with per-step velocity cap.
///        LAMMPS-equivalent `fix nve/limit`.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Velocity Verlet with per-atom velocity magnitude cap.
///
/// Mirrors LAMMPS `fix nve/limit xmax`
/// (`lammps/src/fix_nve_limit.cpp`): standard velocity-Verlet, but after
/// every `v += dtfm * f` kick a per-atom check enforces
/// `|v|_2 <= xmax / dt` -- if exceeded, `v` is rescaled to that
/// magnitude. The cap runs in BOTH halves of the integrator (after the
/// initial half-kick AND after the final half-kick), guaranteeing that
/// the displacement of any atom in any single step never exceeds
/// `xmax` in magnitude.
///
/// Use case: integrate from far-from-equilibrium initial conditions
/// (e.g. compressed amorphous polymer melts) without atoms flying apart
/// from transient force spikes.  Once the system relaxes, the cap is
/// rarely hit and the integrator behaves as plain velocity-Verlet.
///
/// Bound to `initial_integrate @ 0` and `final_integrate @ 0` --
/// replaces (do not combine with) `VelocityVerlet`.
class NVELimit : public cpu::ComponentBase {
public:
    /// @param dt   Timestep [fs].
    /// @param xmax Maximum per-step displacement magnitude [Å]. Implies
    ///             `vmax = xmax / dt`.
    NVELimit(double dt, double xmax);

    std::string name() const override { return "nve_limit"; }
    ComponentRole role() const override { return ComponentRole::Integrator; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), 0},
                {static_cast<int>(ATV::md::Phase::final_integrate), 0}};
    }
    void launch(const cpu::LaunchContext& ctx) override;

private:
    double dt_;
    double xmax_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
