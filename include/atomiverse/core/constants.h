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

/// @brief Angstrom-to-Bohr-radius conversion factor.
///
/// Pinned to PySCF's value `1 / BOHR` with `BOHR = 0.52917721092` Å
/// (CODATA 2006; see `pyscf/data/nist.py`). PySCF is our ground-truth
/// reference for every integral / SCF number, so the Angstrom→Bohr
/// constant that feeds geometry conversion must agree bit-for-bit.
inline constexpr double kAngstromToBohr = 1.8897261245650618;

} // namespace ATV
