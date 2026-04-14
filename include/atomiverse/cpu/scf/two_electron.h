/// @file cpu/scf/two_electron.h
/// @brief Two-electron (ERI) tensor builder — dense-storage path.
///
/// Dense ERI storage is O(N⁴) in memory and only appropriate for small
/// systems, validation, or code paths that must hand the full tensor to
/// another backend (e.g. CUDA RHFFock). The production CPU direct-SCF
/// path uses `build_shell_block_pair_list` in shell_pair.h instead.
///
/// Dispatch order inside the builder mirrors Phase-2 best practice:
///   1. Group AO shells into l-shared ShellBlocks (p-shell batching).
///   2. For each canonical (bra ≤ bra_max, ket ≤ bra) shell-block quartet,
///      call `native::rys_eri_block`.
///   3. Scatter with 8-fold symmetry into the dense tensor.
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/basis_set.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Dense two-electron repulsion tensor ERI(i,j,k,l) = (ij|kl).
///
/// Uses the Rys-quadrature engine with shell-block batching. Output
/// ordering is the physicists' convention (ij|kl) = ∫∫ φ_i(1)φ_j(1)
/// r₁₂⁻¹ φ_k(2)φ_l(2) dτ₁ dτ₂. Symmetry-complete: every permutation
/// of the 8-fold symmetry group is filled.
///
/// @param shells Resolved shell list.
/// @return (n × n × n × n) dense tensor in Hartree.
[[nodiscard]] xt::xtensor<double, 4> build_eri_dense(const ShellList& shells);

} // namespace scf
} // namespace cpu
} // namespace ATV
