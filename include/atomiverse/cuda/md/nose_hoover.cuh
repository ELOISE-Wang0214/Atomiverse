/// @file cuda/md/nose_hoover.cuh
/// @brief NoseHoover -- single-chain Nosé-Hoover thermostat for NVT.
///
/// Modifies velocities at initial_integrate and final_integrate (after VelocityVerlet).
/// Reference: Nosé (1984), Hoover (1985).
/// CUDA-free header.
#pragma once

#include <memory>
#include <string>

#include "atomiverse/cuda/component.cuh"

namespace ATV {
namespace md {

/// @brief Single-chain Nosé-Hoover thermostat.
///
/// Composes with VelocityVerlet (runs after it at same phases via order +100).
/// Maintains thermostat variable ξ on device. Rescales velocities to target temperature.
class NoseHoover : public ComponentBase {
public:
    /// @brief Construct thermostat.
    /// @param T_target Target temperature [K].
    /// @param tau      Coupling time constant [fs].
    /// @param dt       Time step [fs] (must match integrator).
    NoseHoover(Real T_target, Real tau, Real dt);
    ~NoseHoover() override;
    NoseHoover(NoseHoover&&) noexcept;
    NoseHoover& operator=(NoseHoover&&) noexcept;

    std::string name() const override { return "nose_hoover"; }
    std::vector<DefaultBinding> default_bindings() const override;
    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

private:
    Real T_target_, tau_, dt_;
    struct State;
    std::unique_ptr<State> state_;
};

} // namespace md
} // namespace ATV
