/// @file cuda/md/kinetic_energy.cuh
/// @brief md::KineticEnergy -- computes system kinetic energy from velocities.
///
/// KE [eV] = 0.5 * kMvv2e * sum_i( mass_i * (vx_i^2 + vy_i^2 + vz_i^2) ).
///
/// Auto-binds to end_of_step (order -100, before user observers).
/// Accumulates into DeviceFrame::kinetic via shared-memory block reduction
/// + atomicAdd (same pattern as LJ energy accumulation).
///
/// Requires ZeroKinetic to run first (order -1000 at end_of_step) to zero
/// the accumulator before KineticEnergy writes.
///
/// @see ZeroKinetic for the zeroing infrastructure
/// @see Thermo for the component that reads frame.kinetic
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace md {

/// @brief Computes kinetic energy: KE = 0.5 * kMvv2e * sum(m_i * v_i^2) [eV].
///
/// Auto-binds to end_of_step at order -100 (after ZeroKinetic at -1000,
/// before user observers at order >= 0).
///
/// `role()` stays `None` and `capabilities()` stays `None` (the
/// `ComponentBase` defaults), mirroring `cpu::md::KineticEnergy`:
/// kinetic energy is written to `frame.kinetic`, which is **not** one of
/// the `Capability` bits (those cover potential `Energy`, `Force`,
/// `Stress`, charge / density / orbital fields). It is a structural
/// diagnostic, not a force/energy producer, so `DynamicsDriver`
/// validation does not count it toward the `Capability::Force`
/// requirement.
class KineticEnergy : public ComponentBase {
  public:
    /// @brief Component name.
    /// @return "kinetic_energy".
    std::string name() const override { return "kinetic_energy"; }

    /// @brief Default bindings: end_of_step at order -100.
    /// @return Single binding to end_of_step.
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), -100}};
    }

    /// @brief Launch KE reduction kernel.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;
};

} // namespace md
} // namespace cuda
} // namespace ATV
