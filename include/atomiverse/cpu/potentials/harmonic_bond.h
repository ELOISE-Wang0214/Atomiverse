/// @file cpu/potentials/harmonic_bond.h
/// @brief Backwards-compatibility alias header.
///
/// `ATV::HarmonicBond` is now an alias for `ATV::cpu::BondHarmonic` after
/// the 2026-05-19 LAMMPS-naming refactor (`pe-potential-abstraction.md`).
/// New code should use `cpu::BondHarmonic` from
/// `atomiverse/cpu/potentials/bond_harmonic.h` directly. This header is
/// retained so existing call sites (`sim.add(HarmonicBond{...})`) compile
/// unchanged.
#pragma once

#include "atomiverse/cpu/potentials/bond_harmonic.h"

namespace ATV {

using HarmonicBond = ATV::cpu::BondHarmonic;

} // namespace ATV
