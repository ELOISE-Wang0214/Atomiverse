/// @file cpu/md/owned_range.h
/// @brief Compatibility re-export of the backend-root `owned_count()`
///        helpers (moved to `cpu/owned_range.h` — shared infrastructure
///        consumed by md/, minimizer/, and observer/ components).
#pragma once

#include "atomiverse/cpu/owned_range.h"

namespace ATV {
namespace cpu {
namespace md {

using cpu::owned_count;

} // namespace md
} // namespace cpu
} // namespace ATV
