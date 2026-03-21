/// @file cuda/launch.cuh
/// @brief GPU kernel launch configuration helpers. CUDA-free header.
///
/// Centralises the thread-block size so every kernel in the codebase can
/// be re-tuned from a single constant. Also provides a grid_size() helper
/// to eliminate the repeated `(n + block - 1) / block` boilerplate.
#pragma once

namespace ATV {
namespace cuda {

/// @brief Default thread block size for 1D CUDA kernels.
///
/// 256 threads per block balances occupancy and register pressure on
/// Volta / Ampere / Ada / Hopper architectures. Changing this single
/// constant rescales every kernel in the codebase.
inline constexpr int kDefaultBlock = 256;

/// @brief Compute 1D grid dimension for @p n_elements.
/// @return At least 1 block, even when @p n_elements is 0.
inline int grid_size(int n_elements, int block = kDefaultBlock) {
    int g = (n_elements + block - 1) / block;
    return g < 1 ? 1 : g;
}

} // namespace cuda
} // namespace ATV
