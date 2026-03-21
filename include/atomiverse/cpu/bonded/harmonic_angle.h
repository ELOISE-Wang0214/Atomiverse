/// @file cpu/bonded/harmonic_angle.h
/// @brief CPU harmonic angle bending component.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"

namespace ATV {
namespace cpu {
namespace bonded {

/// @brief Harmonic angle bending: V = K(theta - theta0)^2.
///
/// Thread-parallel with per-thread force reduction buffers.
class HarmonicAngle : public cpu::ComponentBase {
public:
    /// @param K Force constants per angle type [eV/rad^2].
    /// @param theta0 Equilibrium angles per angle type [radians].
    HarmonicAngle(std::vector<Real> K, std::vector<Real> theta0);

    std::string name() const override { return "harmonic_angle"; }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
    void destroy(const cpu::DestroyContext& ctx) override;

private:
    std::vector<Real> K_, theta0_;

    int n_threads_ = 1;
    std::vector<std::vector<Real>> fx_thr_, fy_thr_, fz_thr_;
};

} // namespace bonded
} // namespace cpu
} // namespace ATV
