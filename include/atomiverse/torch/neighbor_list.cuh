/// @file torch/neighbor_list.cuh
/// @brief NeighborList for the torch backend (CUDA path).
///
/// Holds ATV::NeighborList internally. build() allocates once; update()
/// rebuilds pair list from current positions. Results are written to
/// frame.pair_* as Tensor (requires_grad=false).
///
/// @see ATV::NeighborList for the CUDA backend component
/// @see torch::Frame for the Tensor-based system state
#pragma once

#ifdef ATV_HAS_TORCH

#include <memory>

#include "atomiverse/torch/frame.cuh"
#include "atomiverse/cuda/context.cuh"
#include "atomiverse/cuda/frame.cuh"

namespace ATV {

class NeighborList;

namespace torch {

/// @brief Neighbor list builder for torch backend (CUDA tensors).
///
/// Wraps ATV::NeighborList. The underlying component is created once
/// in build() and reused across update() calls.
///
/// Usage:
/// @code{.cpp}
/// torch::NeighborList nl(cutoff);
/// nl.build(frame);    // allocate once
/// nl.update(frame);   // rebuild pair list (call after positions change)
/// // frame.pair_* are now populated Tensors
/// @endcode
class NeighborList {
public:
    /// @brief Construct with cutoff and optional skin width.
    /// @param cutoff Interaction cutoff [Angstrom].
    /// @param skin   Verlet skin width [Angstrom]. Negative = auto (0.3 * cutoff).
    explicit NeighborList(double cutoff, double skin = -1.0);
    ~NeighborList();

    NeighborList(NeighborList&&) noexcept;
    NeighborList& operator=(NeighborList&&) noexcept;

    /// @brief Allocate internal resources from frame dimensions.
    /// @param frame torch::Frame with positions, box, and metadata set.
    void build(Frame& frame);

    /// @brief Rebuild pair list from current positions.
    /// @param frame torch::Frame with up-to-date positions.
    void update(Frame& frame);

private:
    double cutoff_;
    double skin_;
    std::unique_ptr<ATV::NeighborList> nl_;
    cuda::Context ctx_;
    cuda::Frame cframe_;
    bool built_ = false;
};

/// @brief Convenience: build + update in one call.
/// @param frame   torch::Frame with positions and box.
/// @param cutoff  Interaction cutoff [Angstrom].
/// @param skin    Verlet skin width [Angstrom]. Negative = auto (0.3 * cutoff).
void neighbor_list(Frame& frame, double cutoff, double skin = -1.0);

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
