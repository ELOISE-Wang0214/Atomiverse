/// @file cpu/semiempirical/internal/bcc_graph.h
/// @brief **PRIVATE** — the molecular graph `AM1BCCCharges` hands to molrs.
#pragma once

#include <array>
#include <string_view>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/am1_bcc_charges.h"

namespace ATV {
namespace cpu {
namespace semiempirical {
namespace internal {

/// @brief The Frame's bonds as local-row index pairs.
///
/// Per-atom topology stores a bond under one or both of its atoms and names
/// the partner by global tag. Each bond comes back once, as `{i, j}` with
/// `i < j`, sorted.
/// @param frame Frame whose local atoms form whole molecules.
/// @return Unique bonds over local rows; empty when the Frame has none.
/// @throws std::invalid_argument when a bond partner tag is not a local atom.
std::vector<std::array<int, 2>> local_bonds(const Frame& frame);

/// @brief molrs name of a correction table.
/// @param set Correction table.
/// @return `"bcc"` or `"abcg2"`, the names `molrs::bcc_model_new` accepts.
std::string_view molrs_name(BccParameterSet set);

} // namespace internal
} // namespace semiempirical
} // namespace cpu
} // namespace ATV
