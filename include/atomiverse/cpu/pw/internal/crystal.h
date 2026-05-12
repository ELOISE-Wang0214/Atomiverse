/// @file cpu/pw/internal/crystal.h
/// @brief Si diamond crystal builder.
#pragma once

#include "atomiverse/cpu/frame.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

cpu::Frame build_si_diamond(Real a = Real(5.43));

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
