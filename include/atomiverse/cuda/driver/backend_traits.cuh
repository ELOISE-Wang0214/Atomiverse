/// @file cuda/driver/backend_traits.cuh
/// @brief CUDA backend_traits specialization.
#pragma once

#include "atomiverse/core/backend_traits.h"
#include "atomiverse/core/device.h"
#include "atomiverse/cuda/driver/backend.cuh"
#include "atomiverse/cuda/component/component.cuh"

namespace ATV {

/// @brief CUDA backend trait specialization.
template <>
struct backend_traits<cuda::device_tag> {
    using component_type = cuda::ComponentBase;
    using backend_type = cuda::Backend;
};

} // namespace ATV
