/// @file cuda/md/barostat.cuh
/// @brief Barostat -- Berendsen-style isotropic pressure coupling.
///
/// Rescales box and positions: μ = (1 - β*dt/τ_p * (P_target - P))^(1/3).
/// Reads stress tensor from DeviceFrame::stress (requires VirialStress to run first).
/// CUDA-free header.
#pragma once

#include <memory>
#include <string>

#include "atomiverse/cuda/component.cuh"

namespace ATV {
namespace md {

/// @brief Berendsen isotropic barostat.
///
/// Binds to final_integrate at order +200 (after VelocityVerlet + NoseHoover).
/// Requires VirialStress at post_force to populate DeviceFrame::stress.
class Barostat : public ComponentBase {
public:
    /// @brief Construct barostat.
    /// @param P_target Target pressure [eV/Angstrom^3].
    /// @param tau_p    Pressure coupling time [fs].
    /// @param dt       Time step [fs].
    Barostat(double P_target, double tau_p, double dt);
    ~Barostat() override;
    Barostat(Barostat&&) noexcept;
    Barostat& operator=(Barostat&&) noexcept;

    std::string name() const override { return "barostat"; }
    std::vector<DefaultBinding> default_bindings() const override;
    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

private:
    double P_target_, tau_p_, dt_;
    struct State;
    std::unique_ptr<State> state_;
};

} // namespace md
} // namespace ATV
