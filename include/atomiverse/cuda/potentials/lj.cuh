/// @file cuda/potentials/lj.cuh
/// @brief LJ -- Lennard-Jones 12-6 pair potential evaluator.
///
/// Computes pairwise forces and energy from the Lennard-Jones potential:
/// @f[ V(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] - V_{\text{shift}} @f]
///
/// Forces are accumulated (atomicAdd) onto the force arrays. Energy is
/// reduced per-block using shared memory, then atomicAdd to the global scalar.
/// Potential is shifted to zero at the cutoff.
///
/// Cross-driver component: returns empty default_bindings(), must be
/// bound explicitly via `sim.on(phase).add(LJ{...})`.
/// CUDA-free header.
///
/// Reference: Jones (1924). On the Determination of Molecular Fields.
/// Proceedings of the Royal Society A, 106(738), 463-477.
#pragma once

#include <string>

#include "atomiverse/cuda/component.cuh"

namespace ATV {

/// @brief Lennard-Jones 12-6 pair potential evaluator.
///
/// Must be explicitly bound to a force phase (cross-driver component).
/// Requires a populated PairList in DeviceFrame (from NeighborList).
class LJ : public ComponentBase {
public:
    /// @brief Construct with LJ parameters.
    /// @param epsilon Well depth [eV].
    /// @param sigma Collision diameter [Angstrom].
    /// @param cutoff Interaction cutoff [Angstrom].
    LJ(double epsilon, double sigma, double cutoff);

    /// @brief Component name.
    /// @return "lj".
    std::string name() const override { return "lj"; }

    /// @brief Per-step force evaluation kernel.
    /// @param ctx Launch context with phase and step information.
    void launch(const LaunchContext& ctx) override;

    /// @brief Get well depth.
    /// @return Well depth [eV].
    double epsilon() const { return epsilon_; }

    /// @brief Get collision diameter.
    /// @return Collision diameter [Angstrom].
    double sigma() const { return sigma_; }

    /// @brief Get interaction cutoff.
    /// @return Interaction cutoff [Angstrom].
    double cutoff() const { return cutoff_; }

private:
    double epsilon_;  ///< Well depth [eV].
    double sigma_;    ///< Collision diameter [Angstrom].
    double cutoff_;   ///< Interaction cutoff [Angstrom].
    double v_shift_;  ///< Potential shift at cutoff [eV].
};

} // namespace ATV
