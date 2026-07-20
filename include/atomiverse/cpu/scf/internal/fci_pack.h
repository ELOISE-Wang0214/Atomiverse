/// @file cpu/scf/internal/fci_pack.h
/// @brief Private pack-based FCI entry.
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/fci.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/scf_result.h" // SCFResult

namespace ATV {
namespace cpu {
namespace scf {

xt::xtensor<double, 2> fci_hamiltonian(const IntegralPack& ints,
                                       const SCFResult& scf);

FCIResult fci_energy(const IntegralPack& ints, const SCFResult& scf);

} // namespace scf
} // namespace cpu
} // namespace ATV
