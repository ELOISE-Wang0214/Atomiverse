/// @file cpu/md/barostat.h
/// @brief CPU Berendsen isotropic barostat.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Berendsen isotropic pressure coupling (Berendsen et al. 1984).
///
/// At each call computes pressure `P = -(σ_xx + σ_yy + σ_zz) / (3·V)` from
/// `frame.stress` (already in unnormalised virial units; this component
/// applies the V division). Then rescales the box edges and every owned
/// atom's positions by `μ = (1 - dt/τ_p · (P - P_target))^(1/3)`.
///
/// Bound at `end_of_step @ -45`, after `ReduceEnergy @ -50`. In multi-rank
/// runs `frame.stress` is the globally-reduced virial after ReduceEnergy,
/// so every rank computes the SAME `μ` and applies the SAME box update --
/// the box stays consistent across ranks without any extra collective.
class Barostat : public cpu::ComponentBase {
public:
    /// @param P_target Target pressure [eV/Å³]. (1 bar ≈ 6.24e-7 eV/Å³.)
    /// @param tau_p    Coupling time constant [fs].
    /// @param dt       MD step size [fs].
    Barostat(double P_target, double tau_p, double dt);

    std::string name() const override { return "barostat"; }
    ComponentRole role() const override { return ComponentRole::Barostat; }
    /// @brief Reads `frame.stress` to derive internal pressure for
    /// the rescaling factor `μ`. Triggers auto-injection of
    /// `VirialStress` in `MD::install_infra()` when no canonical
    /// stress provider is otherwise registered.
    Capability reads() const override { return Capability::Stress; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), -45}};
    }
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Toggle pressure coupling at runtime (default: enabled).
    /// When disabled `launch()` is a no-op — useful for multi-stage
    /// protocols that mix NVT and NPT segments under one driver.
    void set_enabled(bool enabled) noexcept { enabled_ = enabled; }
    bool enabled() const noexcept { return enabled_; }

    /// @brief Retune the pressure setpoint between `Driver::run()` segments.
    void set_P_target(double P_target) noexcept { P_target_ = P_target; }
    double P_target() const noexcept { return P_target_; }

private:
    double P_target_, tau_p_, dt_;
    bool   enabled_ = true;
};

} // namespace md
} // namespace cpu
} // namespace ATV
