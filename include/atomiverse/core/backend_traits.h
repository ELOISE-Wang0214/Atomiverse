/// @file core/backend_traits.h
/// @brief Forward-declared backend trait template. CUDA-free header.
///
/// Each device tag specializes this to export its component_type and backend_type.
/// @see cpu/backend_traits.h, cuda/backend_traits.cuh
#pragma once

namespace ATV {

/// @brief Trait template mapping a device tag to its component and backend types.
/// @tparam DeviceTag Device tag type (cpu::device_tag or cuda::device_tag).
template <typename DeviceTag>
struct backend_traits;

} // namespace ATV
