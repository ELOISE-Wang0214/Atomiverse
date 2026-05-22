/// @file cuda/potentials/multi_harmonic_dihedral.cuh
/// @brief Backwards-compatibility alias header.
#pragma once

#include "atomiverse/cuda/potentials/dihedral_multi_harmonic.cuh"

namespace ATV {
namespace cuda {

using MultiHarmonicDihedral = ATV::cuda::DihedralMultiHarmonic;

} // namespace cuda
} // namespace ATV
