/// @file cuda/md/deform.cuh
/// @brief Deform -- uniaxial strain deformation component.
///
/// Rescales box and atom positions along a given axis at a constant engineering strain rate.
/// CUDA-free header.
#pragma once

#include <string>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cuda {
namespace md {

/// @brief Uniaxial constant-strain-rate deformation.
///
/// Each step: box[axis] *= (1 + rate*dt). Atom positions remapped proportionally.
/// Binds to initial_integrate at order -100 (before VelocityVerlet).
class Deform : public ComponentBase {
  public:
    /// @brief Construct deformation component.
    /// @param axis        Deformation axis (0=x, 1=y, 2=z).
    /// @param strain_rate Engineering strain rate [1/fs].
    /// @param dt          Time step [fs] (must match integrator).
    Deform(int axis, double strain_rate, double dt);

    std::string name() const override { return "deform"; }

    /// @brief Structural role — a box/position time-step integrator.
    ///
    /// Mirrors `cpu::md::Deform`: a strain-rate deformation advances the
    /// box and positions each step, so it occupies the `Integrator`
    /// slot.
    /// @return ComponentRole::Integrator.
    ComponentRole role() const override { return ComponentRole::Integrator; }

    std::vector<DefaultBinding> default_bindings() const override;
    void launch(const LaunchContext& ctx) override;

  private:
    int axis_;
    double strain_rate_, dt_;
};

} // namespace md
} // namespace cuda
} // namespace ATV
