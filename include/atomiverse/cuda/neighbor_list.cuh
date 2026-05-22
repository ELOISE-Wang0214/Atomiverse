/// @file cuda/neighbor_list.cuh
/// @brief NeighborList -- cluster-based pair-building component with full lifecycle.
///
/// Demonstrates the build/launch/destroy lifecycle:
/// - build(): allocates PairList in DeviceFrame + all internal GPU buffers.
/// - launch(): pure GPU pair search + export. Zero D2H/sync/malloc.
/// - destroy(): frees PairList and all internal buffers.
///
/// Cross-driver component: returns empty default_bindings(), must be
/// bound explicitly via `sim.on(phase).add(NeighborList{...})`.
/// CUDA-free header.
#pragma once

#include <memory>
#include <string>

#include "atomiverse/cuda/component.cuh"

namespace ATV {
namespace cuda {

/// @brief Cluster-based neighbor list builder.
///
/// Builds atom-level pair lists using a cluster-based algorithm:
/// morton sort, AABB computation, grid CSR, coarse/fine pair search,
/// and atom-level pair export.
///
/// @see PairList for the output data structure
class NeighborList : public ComponentBase {
  public:
    /// @brief Construct with cutoff and optional skin width.
    /// @param cutoff Interaction cutoff [Angstrom].
    /// @param skin   Verlet skin width [Angstrom]. Negative = auto (0.3 * cutoff).
    explicit NeighborList(double cutoff, double skin = -1.0);

    ~NeighborList() override;
    NeighborList(NeighborList&&) noexcept;
    NeighborList& operator=(NeighborList&&) noexcept;

    /// @brief Component name.
    /// @return "neighbor_list".
    std::string name() const override { return "neighbor_list"; }

    /// @brief Allocate PairList in DeviceFrame + all internal GPU buffers.
    ///
    /// Reads grid dimensions from device (build-time D2H allowed).
    /// Allocates pair capacity based on atom count heuristic (n * 50).
    ///
    /// @param ctx Build context with Runtime and DeviceFrame.
    void build(const BuildContext& ctx) override;

    /// @brief Rebuild pair list from current positions. Pure GPU -- zero D2H.
    ///
    /// Full pipeline: morton sort, cluster AABB, grid CSR, pair search, export.
    ///
    /// @param ctx Launch context with Runtime and DeviceFrame.
    void launch(const LaunchContext& ctx) override;

    /// @brief Free PairList in DeviceFrame + all internal buffers.
    /// @param ctx Destroy context with Runtime and DeviceFrame.
    void destroy(const DestroyContext& ctx) override;

    /// @brief Get interaction cutoff.
    /// @return Interaction cutoff [Angstrom].
    double cutoff() const { return cutoff_; }

    /// @brief Get Verlet skin width.
    /// @return Skin width [Angstrom].
    double skin() const { return skin_; }

  private:
    double cutoff_; ///< Interaction cutoff [Angstrom].
    double skin_;   ///< Verlet skin width [Angstrom].
    struct State;   ///< Opaque GPU state (pimpl, defined in .cu).
    std::unique_ptr<State> state_; ///< Internal GPU buffers.
};

} // namespace cuda
} // namespace ATV
