/// @file cpu/md/virial_stress.h
/// @brief CPU kinetic-stress contribution to `frame.stress`.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Computes the kinetic-stress contribution `Σ m·v_α·v_β` and
/// writes it to `frame.stress_kinetic` (Voigt order: xx, yy, zz, xy,
/// xz, yz).
///
/// Pair virial is accumulated separately by the pair potentials (LJ
/// writes `r_ij · F_ij` into `frame.stress`). The split lets each
/// piece carry its own time slice -- the virial is fixed once forces
/// are computed at `post_force`, but the kinetic stress changes
/// every time velocities are updated. Integrators like
/// `NoseHooverNPT` that need the pressure *after* a velocity-rescale
/// step recompute the kinetic part fresh from current velocities
/// (the trace equals `2·KE`); they read only `frame.stress` (virial)
/// from this component's output. The reduce across MPI ranks happens
/// separately at `end_of_step @ -50` (see `ReduceEnergy`), so any
/// component scheduled after that point reads the global values.
class VirialStress : public cpu::ComponentBase {
public:
    std::string name() const override { return "virial_stress"; }
    ComponentRole role() const override { return ComponentRole::Observer; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::post_force), -50}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace md
} // namespace cpu
} // namespace ATV
