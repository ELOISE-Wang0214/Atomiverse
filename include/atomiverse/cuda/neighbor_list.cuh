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
#include "atomiverse/cuda/types.cuh"

namespace ATV {
namespace cuda {

/// @brief Neighbor pair list -- all arrays on device.
///
/// Computed method-specific scratch, **not** part of `cuda::Frame` (see the
/// CLAUDE.md Frame purity rule). Owned by `NeighborList` -- allocated in
/// `build()`, populated in `launch()`, freed in `destroy()` -- and published
/// to consumers as a non-owning `const PairList*` on `cuda::Context`,
/// mirroring `cpu::md::NeighborList` / `cpu::Context::pair_list`.
struct PairList {
    int* i = nullptr;        ///< First atom index [device].
    int* j = nullptr;        ///< Second atom index [device].
    real_t* dx = nullptr;    ///< x-displacement (j - i) [Angstrom, device].
    real_t* dy = nullptr;    ///< y-displacement (j - i) [Angstrom, device].
    real_t* dz = nullptr;    ///< z-displacement (j - i) [Angstrom, device].
    real_t* dist = nullptr;  ///< Pair distance [Angstrom, device].
    int* n_pairs = nullptr;  ///< Current pair count [device, 1 element].
    int* capacity = nullptr; ///< Pre-allocated slot count [device, 1 element].
};

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

    /// @brief Identifies this component as the neighbor-list supplier
    ///        for `cuda::md::MD::install_infra_`'s typed query.
    bool provides_neighbor_list() const override { return true; }

    /// @brief Allocate the owned PairList + all internal GPU buffers, and
    ///        publish the pair list to `ctx.runtime.pair_list`.
    ///
    /// Reads grid dimensions from device (build-time D2H allowed).
    /// Allocates pair capacity based on atom count heuristic (n * 50).
    ///
    /// @param ctx Build context with Runtime and DeviceFrame.
    void build(const BuildContext& ctx) override;

    /// @brief Rebuild pair list from current positions. Pure GPU -- zero D2H.
    ///
    /// Full pipeline: morton sort, cluster AABB, grid CSR, pair search, export.
    /// Re-publishes the owned pair list to `ctx.runtime.pair_list`.
    ///
    /// @param ctx Launch context with Runtime and DeviceFrame.
    void launch(const LaunchContext& ctx) override;

    /// @brief Free the owned PairList + all internal buffers.
    /// @param ctx Destroy context with Runtime and DeviceFrame.
    void destroy(const DestroyContext& ctx) override;

    /// @brief Access the owned neighbor pair list.
    /// @return Device-resident pair list (valid after `build()`).
    const PairList& pairs() const { return pairs_; }

    /// @brief Get interaction cutoff.
    /// @return Interaction cutoff [Angstrom].
    double cutoff() const { return cutoff_; }

    /// @brief Get Verlet skin width.
    /// @return Skin width [Angstrom].
    double skin() const { return skin_; }

  private:
    double cutoff_; ///< Interaction cutoff [Angstrom].
    double skin_;   ///< Verlet skin width [Angstrom].
    PairList pairs_; ///< Owned device pair list (allocated in build()).
    struct State;   ///< Opaque GPU state (pimpl, defined in .cu).
    std::unique_ptr<State> state_; ///< Internal GPU buffers.
};

} // namespace cuda
} // namespace ATV
