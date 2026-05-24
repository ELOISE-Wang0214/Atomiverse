/// @file cuda/context.cuh
/// @brief cuda::Context -- CUDA execution context (stream handles + host-side box).
///
/// Holds stream handles as void* to keep this header CUDA-free.
/// Components cast to cudaStream_t in their .cu implementations.
///
/// @see cuda::Frame for the device-side data
#pragma once

#include "atomiverse/cuda/box.cuh"
#include "atomiverse/facade/iteration_status.h"

namespace ATV {
namespace cuda {

struct PairList;

/// @brief CUDA execution context passed to components via build/launch/destroy contexts.
///
/// Components must launch kernels onto these streams, never the default stream.
struct Context {
    void* stream = nullptr; ///< Main CUDA stream (cast to cudaStream_t in .cu).
    void* defer_stream =
        nullptr; ///< Secondary stream for deferrable operations.

    /// @brief Non-owning handle to the active neighbor pair list.
    ///
    /// Published by `NeighborList::build()` (and re-published on each
    /// `launch()`) and consumed by downstream force components (LJ) and
    /// observers (RDF). Kept here -- rather than on `cuda::Frame` -- so the
    /// frame stays confined to fundamental physical quantities. Null on
    /// drivers / configurations that carry no neighbor list. Mirrors
    /// `cpu::Context::pair_list`.
    /// @see cuda::PairList
    const PairList* pair_list = nullptr;

    /// @brief Outer-driver loop status, owned by the category driver.
    ///
    /// Reset before each outer iteration; written by
    /// `ConvergenceCheck`-tagged components; polled by `cuda::CUDADriver`
    /// subclasses (`OptimizeDriver` for `converged`, `DynamicsDriver`
    /// for `stop_requested`). Mirrors `cpu::Context::outer_status`.
    /// @see facade/iteration_status.h
    ATV::IterationStatus outer_status;

    /// @brief Host-resident snapshot of the simulation cell, passed by
    ///        value into integrator + bonded + pair kernels so they
    ///        avoid per-thread FP64 division and uncoalesced
    ///        widely-spaced `frame.box[0]/[4]/[8]` loads.
    ///
    /// Populated by `Backend::materialize(cpu::Frame)` at build time
    /// from the canonical `cpu::Box` in `cpu::Frame`. Mutated by
    /// integrators that change the cell (`cuda::md::NoseHooverNPT`,
    /// `cuda::md::Barostat`, `cuda::md::Deform`) after they write the
    /// new diagonal to device-side `frame.box`. All consumers read
    /// from this field at launch — no per-step D2H of box[0,4,8]
    /// (cuda-md-review-cleanup-04-box-by-value-pbc).
    cuda::Box host_box;
};

} // namespace cuda
} // namespace ATV
