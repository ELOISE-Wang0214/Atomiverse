/// @file cuda/transfer.cuh
/// @brief cuda::to_cpu() / cuda::destroy() -- CUDA device transfer functions.
///
/// CUDA-free declarations (implementations in src/cuda/transfer.cu).
///
/// @see cuda::Frame for the CUDA device data
/// @see cpu::Frame for the CPU device data
/// @see cpu::to_cuda() for the reverse direction
#pragma once

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/frame.cuh"

namespace ATV {
namespace cuda {

/// @brief Transfer CUDA Frame back to CPU Frame (D2H copy, no deallocation).
///
/// Copies positions, velocities, forces, energy, kinetic energy, stress,
/// and box from device to host.
///
/// @param src CUDA Frame to read from.
/// @param dst CPU Frame to write into (must be pre-sized).
void to_cpu(const Frame& src, cpu::Frame& dst);

/// @brief Free all device memory in a CUDA Frame.
///
/// Frees all device pointers including PairList and topology arrays.
/// cudaFree(nullptr) is safe, so unallocated fields are handled.
/// Resets the Frame to zero after freeing.
///
/// @param frame CUDA Frame to free (zeroed after this call).
void destroy(Frame& frame);

} // namespace cuda
} // namespace ATV
