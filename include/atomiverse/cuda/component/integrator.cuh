/// @file cuda/component/integrator.cuh
/// @brief cuda::Integrator — time integration (x, v). Kind tag mirroring
///        cpu::Integrator; presets recognise it via dynamic_cast.
#pragma once

#include "atomiverse/cuda/component/component.cuh"

namespace ATV {
namespace cuda {

/// @brief Advances positions and/or velocities.
class Integrator : public Component {};

} // namespace cuda
} // namespace ATV
