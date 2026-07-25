/// @file cpu/dynamics/fix_dfield.h
/// @brief Displacement-field finite-field MD foothold (MetalWalls D-field).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/bias.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Apply \(\mathbf{D}\) as an external field coupling \(U=-V\mathbf{D}\cdot\mathbf{P}\).
///
/// Foothold (vacuum / metal units): treat \(\mathbf{E}\equiv\mathbf{D}\) for
/// force \(F=qE\) on free charges, matching the finite-field electrolyte
/// path when \(\varepsilon=1\). For const-pot walls compose with
/// `ConstPotBO::Ex,Ey,Ez = Dx,Dy,Dz` so electrode \(Q\) sees the same field.
///
/// Full dielectric / Berry-phase D-field machinery is out of scope.
class FixDfield : public cpu::Bias {
  public:
    /// @param Dx Displacement field x [V/Å] (foothold unit match to E)
    /// @param Dy Displacement field y [V/Å]
    /// @param Dz Displacement field z [V/Å]
    FixDfield(double Dx, double Dy, double Dz) noexcept
        : Dx(Dx), Dy(Dy), Dz(Dz) {}

    double Dx = 0.0, Dy = 0.0, Dz = 0.0;

    std::string name() const override { return "fix_dfield"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::md::Phase::post_force, +100)};
    }

    void launch(const cpu::LaunchContext& ctx) override;

    double last_energy() const noexcept { return last_u_; }
    /// @brief Cell polarization proxy \(\sum q_i r_i / V\) z-component [e/Å²]
    double last_Pz() const noexcept { return last_Pz_; }

  private:
    double last_u_ = 0.0;
    double last_Pz_ = 0.0;
};

} // namespace md
} // namespace cpu
} // namespace ATV
