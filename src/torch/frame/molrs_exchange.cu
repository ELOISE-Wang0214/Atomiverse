/// @file torch/frame/molrs_exchange.cu
/// @brief Implementation of torch::from_molrs / torch::to_molrs.
///
/// Both directions stage through a transient `cpu::Frame`: the CPU
/// materializer (`cpu::from_molrs` / `cpu::to_molrs`) handles the
/// molrs.Frame introspection, and `torch::from_cpu` / `torch::to_cpu`
/// perform the explicit host<->device hop.

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/frame/molrs_exchange.cuh"

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/frame/molrs_exchange.h"
#include "atomiverse/torch/frame/transfer.cuh"

namespace ATV {
namespace torch {

Frame from_molrs(const molrs::FrameRef& src) {
    // Stage 1: molrs.Frame -> cpu::Frame (introspection + required-field
    // validation; throws std::runtime_error naming any missing field).
    cpu::Frame cpu_frame = cpu::from_molrs(src);

    // Stage 2: cpu::Frame -> torch::Frame (explicit H2D transfer).
    return from_cpu(cpu_frame);
}

void to_molrs(const Frame& src, molrs::FrameRef& dst) {
    // Stage 1: torch::Frame -> cpu::Frame (explicit D2H transfer). The
    // cpu::Frame must be pre-sized to the torch frame's atom count.
    cpu::Frame cpu_frame(src.n_atoms());
    to_cpu(src, cpu_frame);

    // Stage 2: cpu::Frame -> molrs.Frame (writes per-atom result columns
    // and the system scalars; topology and box untouched).
    cpu::to_molrs(cpu_frame, dst);
}

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
