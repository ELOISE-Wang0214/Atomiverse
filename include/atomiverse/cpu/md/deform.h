/// @file cpu/md/deform.h
/// @brief CPU affine box-deformation driver (single-axis uniaxial strain).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Uniaxial affine deformation.
///
/// At each `initial_integrate` call multiplies the chosen box edge and
/// every owned atom's coordinate along that axis by `1 + ε̇·dt` -- so the
/// engineering strain on that axis ramps linearly at the user-supplied
/// rate `ε̇` per fs. Multi-rank: every rank applies the same deterministic
/// scaling, so the box stays consistent across ranks without any
/// explicit collective.
///
/// Bound at `initial_integrate @ +200`, after VelocityVerlet (+0) and
/// NoseHoover (+100).
class Deform : public cpu::ComponentBase {
public:
    /// @param strain_rate Engineering-strain rate `ε̇` along the chosen
    ///        axis [1 / fs]. Positive = stretch, negative = compress.
    /// @param dt          MD step size [fs] -- needed to convert ε̇ to a
    ///        per-step scale.
    /// @param axis        Box axis to deform: 0 = x, 1 = y, 2 = z.
    Deform(double strain_rate, double dt, int axis = 0);

    std::string name() const override { return "deform"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), 200}};
    }
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Retune the engineering-strain rate between `Driver::run()`
    /// segments. Setting `strain_rate` to 0 disables deformation while
    /// keeping the component registered (no driver rebuild required).
    /// @param strain_rate New ε̇ [1/fs].
    void set_strain_rate(double strain_rate) noexcept { strain_rate_ = strain_rate; }

    /// @return Current ε̇ [1/fs].
    double strain_rate() const noexcept { return strain_rate_; }

private:
    double strain_rate_;   // ε̇ [1/fs]
    double dt_;            // [fs]
    int    axis_;          // 0=x, 1=y, 2=z
};

} // namespace md
} // namespace cpu
} // namespace ATV
