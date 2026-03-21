/// @file cpu/bonded/multi_harmonic_dihedral.h
/// @brief CPU multi-harmonic dihedral component.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"

namespace ATV {
namespace cpu {
namespace bonded {

/// @brief Multi-harmonic dihedral: V = sum(Cn * cos^n(phi)).
///
/// Thread-parallel with per-thread force reduction buffers.
class MultiHarmonicDihedral : public cpu::ComponentBase {
public:
    /// @param Cn Coefficients per dihedral type. Cn[type] = {C0, C1, C2, ...}.
    explicit MultiHarmonicDihedral(std::vector<std::vector<Real>> Cn);

    std::string name() const override { return "multi_harmonic_dihedral"; }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
    void destroy(const cpu::DestroyContext& ctx) override;

private:
    std::vector<std::vector<Real>> Cn_;

    int n_threads_ = 1;
    std::vector<std::vector<Real>> fx_thr_, fy_thr_, fz_thr_;
};

} // namespace bonded
} // namespace cpu
} // namespace ATV
