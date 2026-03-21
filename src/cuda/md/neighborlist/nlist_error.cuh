#pragma once

/// @file nlist_error.cuh
/// @brief CUDA error checking for neighborlist internals.

#include "atomiverse/core/real.h"

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

namespace ATV::md::nlist {

inline void cuda_check(cudaError_t err, const char* file, int line) {
    if (err != cudaSuccess) {
        throw std::runtime_error(
            std::string("CUDA error: ") + cudaGetErrorString(err) +
            " at " + file + ":" + std::to_string(line));
    }
}

#define ATV_CUDA_CHECK(err) ::ATV::md::nlist::cuda_check((err), __FILE__, __LINE__)

} // namespace ATV::md::nlist
