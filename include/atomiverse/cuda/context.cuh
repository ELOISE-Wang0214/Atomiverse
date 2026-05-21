/// @file cuda/context.cuh
/// @brief cuda::Context -- CUDA execution context (stream handles).
///
/// Holds stream handles as void* to keep this header CUDA-free.
/// Components cast to cudaStream_t in their .cu implementations.
///
/// @see cuda::Frame for the device-side data
#pragma once

namespace ATV {
namespace cuda {

/// @brief CUDA execution context passed to components via build/launch/destroy contexts.
///
/// Components must launch kernels onto these streams, never the default stream.
struct Context {
    void* stream = nullptr; ///< Main CUDA stream (cast to cudaStream_t in .cu).
    void* defer_stream =
        nullptr; ///< Secondary stream for deferrable operations.
};

} // namespace cuda
} // namespace ATV
