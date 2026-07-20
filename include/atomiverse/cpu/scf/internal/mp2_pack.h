/// @file cpu/scf/internal/mp2_pack.h
/// @brief Private pack-based MP2 entry (unit tests / method implementation).
#pragma once

#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/scf_result.h" // SCFResult

namespace ATV {
namespace cpu {
namespace scf {

MP2Result mp2_energy(const IntegralPack& ints, const SCFResult& scf);

} // namespace scf
} // namespace cpu
} // namespace ATV
