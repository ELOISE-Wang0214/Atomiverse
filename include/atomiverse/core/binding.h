/// @file core/binding.h
/// @brief DefaultBinding -- shared between CPU and CUDA component bases. CUDA-free header.
#pragma once

namespace ATV {

/// @brief Default stage binding declared by a component.
struct DefaultBinding {
    int phase; ///< Target phase index.
    int order; ///< Ordering within phase (lower = earlier).
};

} // namespace ATV
