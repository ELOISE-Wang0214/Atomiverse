/// @file cpu/bonded/harmonic_bond.h
/// @brief CPU harmonic bond stretching component.
#pragma once

#include <cstdint>
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
    HarmonicBond(std::vector<double> K, std::vector<double> r0);

    std::string name() const override { return "harmonic_bond"; }
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force;
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
    void destroy(const cpu::DestroyContext& ctx) override;

private:
    std::vector<double> K_, r0_;

    int n_threads_ = 1;
    std::vector<std::vector<double>> fx_thr_, fy_thr_, fz_thr_;

    std::vector<int> local_bonds_;
    std::uint64_t local_bonds_generation_ = ~std::uint64_t{0};
};

} // namespace bonded
} // namespace cpu
} // namespace ATV
