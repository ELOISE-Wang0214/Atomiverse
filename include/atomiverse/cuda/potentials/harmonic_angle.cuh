/// @file cuda/potentials/harmonic_angle.cuh
/// @brief HarmonicAngle -- harmonic angle bending potential.
///
/// V(θ) = K (θ - θ0)^2,  angle i-j-k where j is center atom.
///
/// Cross-driver component: returns empty default_bindings(), bind at force phase.
/// Reference: standard harmonic angle, Dreiding FF (Mayo et al. 1990).
/// CUDA-free header.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"

namespace ATV {

/// @brief Harmonic angle bending potential.
///
/// One thread per angle. Forces via chain rule dV/dθ · dθ/dr.
class HarmonicAngle : public ComponentBase {
public:
    /// @brief Construct with per-type parameters.
    /// @param K      Spring constants per angle type [eV/rad^2].
    /// @param theta0 Equilibrium angles per type [degrees] (converted to radians internally).
    HarmonicAngle(std::vector<Real> K, std::vector<Real> theta0);
    ~HarmonicAngle() override;
    HarmonicAngle(HarmonicAngle&&) noexcept;
    HarmonicAngle& operator=(HarmonicAngle&&) noexcept;

    std::string name() const override { return "harmonic_angle"; }
    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

private:
    std::vector<Real> K_, theta0_;
    struct State;
    std::unique_ptr<State> state_;
};

} // namespace ATV
