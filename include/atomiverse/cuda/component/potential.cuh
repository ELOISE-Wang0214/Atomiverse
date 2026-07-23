/// @file cuda/component/potential.cuh
/// @brief cuda::Potential — classical force term (pair, bond, coul, kspace, …).
///        Kind tag mirroring cpu::Potential; presets recognise it via
///        dynamic_cast.
#pragma once

#include "atomiverse/cuda/component/component.cuh"

namespace ATV {
namespace cuda {

/// @brief Classical potential energy / force term.
class Potential : public Component {};

} // namespace cuda
} // namespace ATV
