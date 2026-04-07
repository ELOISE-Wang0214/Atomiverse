/// @file core/constants.h
/// @brief Mathematical and unit-conversion constants (non-physical).
///
/// Physical constants (kBoltzmann, kFtm2v, kMvv2e) live in core/units.h.
/// This header covers pure-math and unit-conversion factors that are not
/// part of the ATV physical unit system.
#pragma once

namespace ATV {

/// @brief Pi (double-precision literal).
inline constexpr double kPI = 3.14159265358979323846;

/// @brief Degrees-to-radians conversion factor.
inline constexpr double kDegToRad = kPI / 180.0;

/// @brief Angstrom-to-Bohr-radius conversion factor (CODATA 2018).
inline constexpr double kAngstromToBohr = 1.8897259886;

} // namespace ATV
