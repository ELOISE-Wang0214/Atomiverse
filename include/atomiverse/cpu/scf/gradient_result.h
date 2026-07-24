/// @file cpu/scf/gradient_result.h
/// @brief `GradientResult` — the public forces-API return type (rehomed from the
///        since-deleted scf_result.h; the public gradient.h stays a
///        tombstone).
#pragma once

#include <utility>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/units.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Constructor tag selecting the Ha/Bohr → eV/Å converting overload of
///        `GradientResult` (disambiguates the gradient-in-atomic-units input
///        from the eV/Å member-wise ctor; cf. `SlaterDeterminant::Generalized`).
struct FromHartreeBohr {};

/// @brief Analytic / FD nuclear gradient result (public forces API).
struct GradientResult {
    xt::xtensor<double, 2> forces; ///< [n_atoms, 3] eV/Å
    double energy = 0;             ///< eV

    /// @brief Empty result (zero energy [eV], empty forces).
    GradientResult() = default;

    /// @brief Member-wise construct from values already in the public currency.
    /// @param forces_ev Forces, shape [n_atoms, 3] [eV/Å]. Stored verbatim.
    /// @param energy_ev Energy [eV]. Stored verbatim.
    GradientResult(xt::xtensor<double, 2> forces_ev, double energy_ev)
        : forces(std::move(forces_ev)), energy(energy_ev) {}

    /// @brief Convert an energy gradient in atomic units to the public forces
    ///        currency: flip the sign (F = −∂E/∂R) and change units
    ///        (Ha/Bohr → eV/Å, Hartree → eV) via `units::convert` as the sole
    ///        primitive. `n_atoms` is read from `grad_ha_bohr.shape()[0]`.
    /// @param grad_ha_bohr Energy gradient ∂E/∂R, shape [n_atoms, 3] [Ha/Bohr].
    /// @param energy_ha    Total energy [Hartree].
    GradientResult(FromHartreeBohr, const xt::xtensor<double, 2>& grad_ha_bohr,
                   double energy_ha) {
        const int natm = static_cast<int>(grad_ha_bohr.shape()[0]);
        forces = xt::zeros<double>({natm, 3});
        for (int A = 0; A < natm; ++A)
            for (int d = 0; d < 3; ++d)
                forces(A, d) = -units::convert(grad_ha_bohr(A, d),
                                               units::Force::HartreePerBohr,
                                               units::Force::eVPerAngstrom);
        energy =
            units::convert(energy_ha, units::Energy::Hartree, units::Energy::eV);
    }
};

} // namespace scf
} // namespace cpu
} // namespace ATV
