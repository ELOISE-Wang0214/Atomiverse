/// @file cpu/vib/displace_step.h
/// @brief cpu::vib::DisplaceStep — applies the +/-h Cartesian displacement.
///
/// Bound at `vib::Phase::displace`. For displacement index `s = ctx.step`
/// (`0 <= s < 6N`) it decodes the Cartesian coordinate and sign and sets
/// the frame coordinates to `R0 +/- h*e_{i.alpha}`:
///
/// @code
///   coordinate c = s / 2        (0 <= c < 3N)   atom i = c/3, axis = c%3
///   sign       = (s % 2 == 0) ? +1 : -1
/// @endcode
///
/// On the final finalize pass (`s == 6N`) it restores the frame to the
/// cached equilibrium geometry `R0`.
///
/// Displacement is applied directly in the frame's Angstrom coordinates
/// with `h = kHessianStepAngstrom` (1.0e-3 Angstrom).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/vib/hessian_workspace.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace vib {

/// @brief Structural component that displaces one Cartesian coordinate.
///
/// `role() == None`; it mutates positions (an implicit frame field) but
/// declares no `capabilities()` (Energy/Force/... are unaffected).
class DisplaceStep : public cpu::ComponentBase {
  public:
    /// @brief Construct bound to a driver-owned workspace.
    /// @param ws Non-owning pointer to the shared Hessian workspace.
    explicit DisplaceStep(HessianWorkspace* ws) : ws_(ws) {}

    std::string name() const override { return "vib_displace_step"; }
    ComponentRole role() const override { return ComponentRole::None; }

    /// @brief Default binding: `vib::Phase::displace` @ 0.
    /// @return Single binding placing this at the displace phase.
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::vib::Phase::displace), 0}};
    }

    /// @brief Set frame coordinates to `R0 +/- h*e_{i.alpha}` (or restore `R0`).
    /// @param ctx Launch context; `ctx.step` is the displacement index `s`.
    void launch(const cpu::LaunchContext& ctx) override {
        const int n = ws_->n_atoms;
        const int dim = 3 * n;

        // Restore every coordinate to R0 first — this both clears the
        // previous displacement and handles the final (s == 6N) restore pass.
        for (int i = 0; i < n; ++i) {
            ctx.frame.atoms.x(i) = ws_->r0[3 * i + 0];
            ctx.frame.atoms.y(i) = ws_->r0[3 * i + 1];
            ctx.frame.atoms.z(i) = ws_->r0[3 * i + 2];
        }

        const int s = ctx.step;
        if (s >= 2 * dim)
            return; // finalize pass: leave the frame at R0.

        const int coord = s / 2;
        const double sign = (s % 2 == 0) ? +1.0 : -1.0;
        const int atom = coord / 3;
        const int axis = coord % 3;
        const double shift = sign * kHessianStepAngstrom;

        if (axis == 0)
            ctx.frame.atoms.x(atom) += shift;
        else if (axis == 1)
            ctx.frame.atoms.y(atom) += shift;
        else
            ctx.frame.atoms.z(atom) += shift;
    }

  private:
    HessianWorkspace* ws_; ///< Non-owning; driver-owned.
};

} // namespace vib
} // namespace cpu
} // namespace ATV
