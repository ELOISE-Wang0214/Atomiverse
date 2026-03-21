/// @file cpu/bonded/harmonic_bond.h
/// @brief CPU harmonic bond stretching component.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"

namespace ATV {
namespace cpu {
namespace bonded {

/// @brief Harmonic bond stretching: V = K(r - r0)^2.
///
/// Thread-parallel with per-thread force reduction buffers (Newton's 3rd law).
class HarmonicBond : public cpu::ComponentBase {
public:
    /// @param K Force constants per bond type [eV/Angstrom^2].
    /// @param r0 Equilibrium distances per bond type [Angstrom].
    HarmonicBond(std::vector<Real> K, std::vector<Real> r0);

    std::string name() const override { return "harmonic_bond"; }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
    void destroy(const cpu::DestroyContext& ctx) override;

private:
    std::vector<Real> K_, r0_;

    int n_threads_ = 1;
    std::vector<std::vector<Real>> fx_thr_, fy_thr_, fz_thr_;
};

} // namespace bonded
} // namespace cpu
} // namespace ATV
