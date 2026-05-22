/// @file cuda/potentials/harmonic_bond.cuh
/// @brief Backwards-compatibility alias header.
///
/// `ATV::HarmonicBond` is now an alias for `ATV::cuda::BondHarmonic`
/// after the 2026-05-19 LAMMPS-naming refactor.
#pragma once

#include "atomiverse/cuda/potentials/bond_harmonic.cuh"

namespace ATV {
namespace cuda {

using HarmonicBond = ATV::cuda::BondHarmonic;

} // namespace cuda
} // namespace ATV
