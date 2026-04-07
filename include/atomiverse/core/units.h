/// @file core/units.h
/// @brief Physical unit conversion constants for the ATV unit system.
///
/// Unit system: eV (energy), Angstrom (distance), fs (time), amu (mass).
///
/// Derived: eV/Angstrom (force), Angstrom/fs (velocity), Angstrom/fs^2 (accel).
#pragma once

namespace ATV {

/// @brief Force-to-motion: multiply (force/mass) by this to get acceleration in Angstrom/fs^2.
/// @see Allen & Tildesley (2017), Table 3.1.
/// @see LAMMPS metal units: force->ftm2v.
constexpr double kFtm2v = 1.0 / 103.6427;

/// @brief Mass-velocity-velocity to energy conversion.
/// KE [eV] = 0.5 * kMvv2e * mass [amu] * v^2 [Angstrom^2/fs^2].
constexpr double kMvv2e = 103.6427;

/// @brief Boltzmann constant in eV/K (CODATA 2018).
constexpr double kBoltzmann = 8.617333262e-5;

/// @brief Hartree to eV conversion (CODATA 2018).
constexpr double kHartree2eV = 27.211386245988;

/// @brief Bohr to Angstrom conversion (CODATA 2018).
constexpr double kBohr2Ang = 0.529177210903;

/// @brief Force conversion: Hartree/Bohr to eV/Angstrom.
constexpr double kForceAU2Metal = kHartree2eV / kBohr2Ang;

} // namespace ATV
