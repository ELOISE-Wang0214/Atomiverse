/// @file core/pair_list.h
/// @brief Canonical PairList SoA field contract.
///
/// Defines the field names, types, and semantics that all backend PairList
/// implementations must provide. This is a compile-time contract, not a base
/// class — each backend implements PairList with its own storage type:
///   - cpu::PairList  — xtensor arrays
///   - cuda::PairList — device pointers (in cuda/frame.cuh)
///   - torch::Frame   — at::Tensor fields
///
/// @par Canonical SoA fields
/// | Field    | Type | Shape      | Description                         |
/// |----------|------|------------|-------------------------------------|
/// | i        | Int  | [n_pairs]  | First atom index                    |
/// | j        | Int  | [n_pairs]  | Second atom index                   |
/// | dx       | Real | [n_pairs]  | x displacement j-i [Angstrom]       |
/// | dy       | Real | [n_pairs]  | y displacement j-i [Angstrom]       |
/// | dz       | Real | [n_pairs]  | z displacement j-i [Angstrom]       |
/// | dist     | Real | [n_pairs]  | Pair distance [Angstrom]            |
/// | n_pairs  | Int  | scalar     | Number of valid pairs               |
///
/// @par Conventions
/// - Half-list: each pair (i,j) appears once with i < j.
/// - Displacement: d = pos[j] - pos[i], minimum image.
/// - Units: Angstrom for all distances.
#pragma once

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

namespace ATV {

/// @brief Number of SoA array fields in the canonical PairList layout.
inline constexpr int kPairListFieldCount = 6;  // i, j, dx, dy, dz, dist

} // namespace ATV
