/// @file cpu/potentials/harmonic_angle.h
/// @brief Backwards-compatibility alias header.
///
/// `ATV::HarmonicAngle` is now an alias for `ATV::cpu::AngleHarmonic`
/// after the 2026-05-19 LAMMPS-naming refactor.
#pragma once

#include "atomiverse/cpu/potentials/angle_harmonic.h"

namespace ATV {

using HarmonicAngle = ATV::cpu::AngleHarmonic;

} // namespace ATV
