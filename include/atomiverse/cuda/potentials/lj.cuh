/// @file cuda/potentials/lj.cuh
/// @brief Backwards-compatibility alias header.
///
/// `ATV::LJ` is now an alias for `ATV::cuda::PairLJCut` after the
/// 2026-05-19 LAMMPS-naming refactor (`pe-potential-abstraction.md`).
/// New code should use `cuda::PairLJCut` from
/// `atomiverse/cuda/potentials/pair_lj_cut.cuh` directly.
#pragma once

#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"

namespace ATV {
namespace cuda {

using LJ = ATV::cuda::PairLJCut;

} // namespace cuda
} // namespace ATV
