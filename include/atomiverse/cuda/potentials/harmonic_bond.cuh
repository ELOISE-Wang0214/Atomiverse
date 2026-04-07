/// @file cuda/potentials/harmonic_bond.cuh
/// @brief HarmonicBond -- harmonic bond stretching potential.
///
/// V(r) = K (r - r0)^2,  F = -2K(r - r0) r_hat
///
/// Cross-driver component: returns empty default_bindings(), bind at force phase.
/// Reference: standard harmonic bond, Dreiding FF (Mayo et al. 1990).
/// CUDA-free header.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"

namespace ATV {

/// @brief Harmonic bond stretching potential.
///
/// One thread per bond. Reads BondList from DeviceFrame.
/// Forces accumulated via atomicAdd (Newton's 3rd law).
///
/// @see BondList for the topology data structure
class HarmonicBond : public ComponentBase {
public:
    /// @brief Construct with per-type parameters.
    /// @param K  Spring constants per bond type [eV/Angstrom^2].
    /// @param r0 Equilibrium distances per bond type [Angstrom].
    HarmonicBond(std::vector<double> K, std::vector<double> r0);
    ~HarmonicBond() override;
    HarmonicBond(HarmonicBond&&) noexcept;
    HarmonicBond& operator=(HarmonicBond&&) noexcept;

    std::string name() const override { return "harmonic_bond"; }
    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

private:
    std::vector<double> K_, r0_;
    struct State;
    std::unique_ptr<State> state_;
};

} // namespace ATV
