/// @file cuda/backend_traits.cuh
/// @brief CUDA backend_traits specialization.
#pragma once

#include "atomiverse/core/device.h"
#include "atomiverse/core/backend_traits.h"
#include "atomiverse/cuda/component.cuh"
#include "atomiverse/cuda/backend.cuh"

namespace ATV {

/// @brief CUDA backend trait specialization.
template <>
struct backend_traits<cuda::device_tag> {
    using component_type = ComponentBase;
    using backend_type   = Backend;
};

} // namespace ATV
