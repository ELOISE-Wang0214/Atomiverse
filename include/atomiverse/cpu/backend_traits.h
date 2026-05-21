/// @file cpu/backend_traits.h
/// @brief CPU backend_traits specialization. CUDA-free header.
#pragma once

#include "atomiverse/core/backend_traits.h"
#include "atomiverse/core/device.h"
#include "atomiverse/cpu/backend.h"
#include "atomiverse/cpu/component.h"

namespace ATV {

/// @brief CPU backend trait specialization.
template <>
struct backend_traits<cpu::device_tag> {
    using component_type = cpu::ComponentBase;
    using backend_type = cpu::Backend;
};

} // namespace ATV
