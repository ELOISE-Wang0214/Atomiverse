/// @file cuda/potentials/multi_harmonic_dihedral.cuh
/// @brief MultiHarmonicDihedral -- multi-harmonic (OPLS-style) dihedral potential.
///
/// V(φ) = Σ_{n=0}^{4} C_n cos^n(φ)
///
/// Cross-driver component: returns empty default_bindings(), bind at force phase.
/// Reference: Dreiding FF (Mayo et al. 1990), OPLS-style multi-harmonic.
/// CUDA-free header.
#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"

namespace ATV {

/// @brief Multi-harmonic dihedral torsion potential.
///
/// One thread per dihedral. Computes φ from 4 atoms, V = Σ C_n cos^n(φ).
class MultiHarmonicDihedral : public ComponentBase {
public:
    /// @brief Construct with per-type coefficients.
    /// @param C Per-type arrays of 5 coefficients {C0..C4} [eV].
    MultiHarmonicDihedral(std::vector<std::array<Real, 5>> C);
    ~MultiHarmonicDihedral() override;
    MultiHarmonicDihedral(MultiHarmonicDihedral&&) noexcept;
    MultiHarmonicDihedral& operator=(MultiHarmonicDihedral&&) noexcept;

    std::string name() const override { return "multi_harmonic_dihedral"; }
    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

private:
    std::vector<std::array<Real, 5>> C_;
    struct State;
    std::unique_ptr<State> state_;
};

} // namespace ATV
