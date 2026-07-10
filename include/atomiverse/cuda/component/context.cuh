/// @file cuda/component/context.cuh
/// @brief cuda::Context -- CUDA execution context: flat core state (streams
/// + host-side box + outer status) plus the classical force-pipeline slice
/// (arch-refactor-05-context).
///
/// Holds stream handles as void* to keep this header CUDA-free.
/// Components cast to cudaStream_t in their .cu implementations.
///
/// Homomorphic with `cpu::Context`: same-named slices carry the same
/// semantics. CUDA has no Sampling / Optimize driver and no analysis
/// executor (async observation goes through `AsyncBridge`), so only the
/// `md` slice exists on this backend — homomorphism means "same slice,
/// same meaning", not "same field set".
///
/// @see cuda::Frame for the device-side data
#pragma once

#include "atomiverse/cuda/frame/box.cuh"
#include "atomiverse/facade/iteration_status.h"

namespace ATV {
namespace cuda {

struct PairList;

/// @brief Classical force-pipeline services (publishers: `NeighborList`,
///        `PairLJCut::build`, the multi-GPU comm layer). Mirrors
///        `cpu::DynamicsServices`; the exclusion CSR and the force-row
///        split live here because they belong to the same
///        publish-consume chain as `pair_list`.
struct DynamicsServices {
    /// @brief Non-owning handle to the active neighbor pair list.
    ///
    /// Published by `NeighborList::build()` (and re-published on each
    /// `launch()`) and consumed by downstream force components (LJ) and
    /// observers. Kept here -- rather than on `cuda::Frame` -- so the
    /// frame stays confined to fundamental physical quantities. Null on
    /// drivers / configurations that carry no neighbor list. Mirrors
    /// `cpu::DynamicsServices::pair_list`.
    /// @see cuda::PairList
    const PairList* pair_list = nullptr;

    /// @brief Opt-in request for the flat (i,j,dx,dy,dz,r²) atom pair list.
    ///
    /// The atom-centric LJ force kernel reads the cluster CSR directly and
    /// does NOT need the flat list, so NeighborList skips the per-step
    /// prune+export (~60% of step time) by default. Consumers that
    /// genuinely need the flat list set this true in their build();
    /// NeighborList::launch() then emits it every step. Set before the
    /// first launch (all build()s precede any launch()).
    bool flat_pairs_requested = false;

    /// @brief Special-bonds exclusion CSR (per-atom, LOCAL indices),
    ///        published by `PairLJCut::build` and consumed by
    ///        `NeighborList::launch` to bake exclusions directly into the
    ///        per-atom Verlet list (LAMMPS style — find_special at build
    ///        time, not per pair at force time).
    ///
    /// `excl_offsets[i]..excl_offsets[i+1]` index `excl_partners` for atom
    /// i. `s_excl` is the single special-bonds scale (set_special_bonds
    /// requires all tiers equal). When `s_excl == 0` and the neighbor list
    /// bakes the exclusions (sets `exclusions_baked`), the force kernel
    /// skips the per-pair exclusion scan entirely. Null/false on configs
    /// without special bonds.
    const int* excl_offsets = nullptr;
    const int* excl_partners = nullptr;
    double s_excl = 1.0;
    bool exclusions_baked = false;

    /// @brief Pair-force row restriction (multi-GPU owned/ghost split).
    ///
    /// A pair style computes one force row per atom. On a multi-GPU
    /// subdomain the frame's atom range is [owned | ghost]; ghost rows
    /// must NOT be computed: each owned atom's force is already complete
    /// from its own full neighbor list, so a ghost row's pair force is (a)
    /// wasted work and (b) WRONG when `reverse_comm` later ships ghost
    /// rows home — the pair part would double onto the owner (only bonded
    /// ghost deposits may travel). It also double-counts cross-boundary
    /// pair energy (each half-pair contribution appears on both GPUs).
    ///
    /// `n_force_rows < 0` (default): compute every atom row (single-GPU /
    /// driver path — unchanged). `n_force_rows >= 0` with `force_rows ==
    /// nullptr`: compute the owned prefix [0, n_force_rows). With
    /// `force_rows != nullptr`: compute the explicit row subset
    /// `force_rows[0..n_force_rows)` — the multi-GPU interior/boundary
    /// split that overlaps interior force with the halo exchange. Neighbor
    /// POSITIONS are always read from the full owned+ghost set regardless.
    int n_force_rows = -1;
    const int* force_rows = nullptr;
};

/// @brief CUDA execution context passed to components via
///        build/launch/destroy contexts.
///
/// Components must launch kernels onto these streams, never the default
/// stream.
struct Context {
    void* stream = nullptr; ///< Main CUDA stream (cast to cudaStream_t in .cu).
    void* defer_stream =
        nullptr; ///< Secondary stream for deferrable operations.

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
    /// (cuda-md-review-cleanup-04-box-by-value-pbc). Core state, not a
    /// nullable service — stays flat.
    cuda::Box host_box;

    DynamicsServices md; ///< Classical force-pipeline slice.
};

} // namespace cuda
} // namespace ATV
