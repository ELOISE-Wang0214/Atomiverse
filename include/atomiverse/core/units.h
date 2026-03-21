/// @file core/units.h
/// @brief Physical unit conversion constants for the ATV unit system.
///
/// Unit system: eV (energy), Angstrom (distance), fs (time), amu (mass).
///
/// Derived: eV/Angstrom (force), Angstrom/fs (velocity), Angstrom/fs^2 (accel).
#pragma once

#include "atomiverse/core/real.h"

namespace ATV {

/// @brief Force-to-motion: multiply (force/mass) by this to get acceleration in Angstrom/fs^2.
/// @see Allen & Tildesley (2017), Table 3.1.
/// @see LAMMPS metal units: force->ftm2v.
constexpr Real kFtm2v = Real(1) / Real(103.6427);

/// @brief Mass-velocity-velocity to energy conversion.
/// KE [eV] = 0.5 * kMvv2e * mass [amu] * v^2 [Angstrom^2/fs^2].
constexpr Real kMvv2e = Real(103.6427);

/// @brief Boltzmann constant in eV/K (CODATA 2018).
constexpr Real kBoltzmann = Real(8.617333262e-5);

/// @brief Hartree to eV conversion (CODATA 2018).
constexpr Real kHartree2eV = Real(27.211386245988);

/// @brief Bohr to Angstrom conversion (CODATA 2018).
constexpr Real kBohr2Ang = Real(0.529177210903);

/// @brief Force conversion: Hartree/Bohr to eV/Angstrom.
constexpr Real kForceAU2Metal = kHartree2eV / kBohr2Ang;

} // namespace ATV
