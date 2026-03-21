/// @file core/device.h
/// @brief Device tags for compile-time backend selection. CUDA-free header.
#pragma once

namespace ATV {
namespace cpu  { struct device_tag {}; }
namespace cuda { struct device_tag {}; }
} // namespace ATV
