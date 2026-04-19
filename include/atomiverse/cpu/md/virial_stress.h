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

/// @brief Adds the kinetic-stress contribution `Σ m·v_α·v_β` to
/// `frame.stress` (Voigt order: xx, yy, zz, xy, xz, yz).
///
/// Pair virial is accumulated separately by the pair potentials (LJ
/// writes `r_ij · F_ij` into `frame.stress` in its launch). When this
/// component runs at `post_force @ -50`, both contributions are summed
/// in `frame.stress` -- still in **unnormalised virial units**: callers
/// (Barostat, observers) divide by box volume V to obtain pressure.
///
/// The reduce across MPI ranks happens separately at `end_of_step @
/// -50` (see `ReduceEnergy`), so a Barostat scheduled after that point
/// reads the global virial and computes the global pressure.
class VirialStress : public cpu::ComponentBase {
public:
    std::string name() const override { return "virial_stress"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::post_force), -50}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace md
} // namespace cpu
} // namespace ATV
