/// @file cpu/potentials/multi_harmonic_dihedral.h
/// @brief Backwards-compatibility alias header.
///
/// `ATV::MultiHarmonicDihedral` is now an alias for
/// `ATV::cpu::DihedralMultiHarmonic` after the 2026-05-19 LAMMPS-naming
/// refactor.
#pragma once

#include "atomiverse/cpu/potentials/dihedral_multi_harmonic.h"

namespace ATV {

using MultiHarmonicDihedral = ATV::cpu::DihedralMultiHarmonic;

} // namespace ATV
