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
    int* i = nullptr;     ///< First atom index [device].
    int* j = nullptr;     ///< Second atom index [device].
    real_t* dx = nullptr; ///< x-displacement (j - i) [Angstrom, device].
    real_t* dy = nullptr; ///< y-displacement (j - i) [Angstrom, device].
    real_t* dz = nullptr; ///< z-displacement (j - i) [Angstrom, device].
    real_t* dist_sq =
        nullptr; ///< Squared pair distance dx²+dy²+dz² [Å², device].
    ///<
    ///< Per cuda-md-review-cleanup-06-launch-config-streams, the
    ///< neighbor-list export stores r² rather than r — pair styles
    ///< need inv_r² for the LJ force calculation and can skip the
    ///< sqrt entirely. Consumers that need r (e.g. cuda::RDF for
    ///< histogram binning) take sqrt at read time.
    int* n_pairs = nullptr;  ///< Current pair count [device, 1 element].
    int* capacity = nullptr; ///< Pre-allocated slot count [device, 1 element].

    // ── Cluster neighbor list (atom-centric, atomic-free force path) ──────
    // The builder internally produces a per-cluster CSR. These expose it so a
    // force kernel can own each atom (one thread per atom slot), traverse its
    // cluster's FULL neighbor-cluster list, compute distances from CURRENT
    // atom positions, accumulate the atom's force in registers, and write it
    // ONCE — no atomics, no flat pair list, no per-step export. row_ptr/col_idx
    // are FULL (every neighbor cluster, both directions). Addresses are stable
    // across the run; rebuilt in place on each skin-gated rebuild. orig_index
    // maps a (cluster*kClusterSize + slot) to the global atom index; pad slots
    // hold a negative sentinel. Null until build() publishes them.
    const int* cl_row_ptr = nullptr;    ///< [num_clusters+1] CSR row offsets.
    const int* cl_col_idx = nullptr;    ///< neighbor cluster indices (full).
    const int* cl_orig_index = nullptr; ///< cluster slot -> global atom index.
    int num_clusters = 0;               ///< cluster count (host scalar).
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
    /// Pair capacity is sized from physical density × cutoff sphere volume
    /// (`max(100, 3 × ρ × (4/3)π(r_cut+skin)³) per atom`); falls back to
    /// `n × 50` only for atomless / zero-volume frames. See `neighbor_list.cu`
    /// `NeighborList::build` for the density-based rationale.
    ///
    /// @param ctx Build context with Runtime and DeviceFrame.
    void build(const BuildContext& ctx) override;

    /// @brief Rebuild pair list from current positions, gated by the
    ///        Verlet-skin rebuild check.
    ///
    /// Pipeline (run only when `rebuild_flag` from the previous step is
    /// set): morton sort → CUB radix sort → fused scatter+AABB → grid
    /// CSR → coarse count → fill → compact → csr_to_flat → snapshot
    /// ref_x/y/z. Every step (rebuild or skip): mask refresh via
    /// `smart_rolling_prune` → atom-pair export → `atom_displacement_max`
    /// for the next step's rebuild gate.
    ///
    /// First launch always rebuilds (no prior flag exists). Subsequent
    /// launches use a 1-step-lagged async D2H of `rebuild_flag` —
    /// `skin/2` width absorbs the extra step of motion (standard
    /// LAMMPS / GROMACS pattern).
    ///
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

    /// @brief Number of full neighbor-list rebuilds since `build()`.
    ///
    /// Test introspection only — used by `test_neighborlist_skin_gating`
    /// to assert that a quiescent system rebuilds at most O(1) times
    /// over many steps. Production code should not depend on this
    /// counter (the rebuild cadence is a perf-tuning detail, not a
    /// public contract).
    /// @return Rebuild count over the lifetime of this component.
    int rebuild_count() const noexcept;

  private:
    double cutoff_; ///< Interaction cutoff [Angstrom].
    double skin_;   ///< Verlet skin width [Angstrom].
    PairList pairs_; ///< Owned device pair list (allocated in build()).
    struct State;   ///< Opaque GPU state (pimpl, defined in .cu).
    std::unique_ptr<State> state_; ///< Internal GPU buffers.
};

} // namespace cuda
} // namespace ATV
