/// @file cpu/md/mpi/halo_exchange.h
/// @brief Ghost-atom halo exchange components for MPI domain decomposition.
///
/// Phase 2 implementation:
///   - `HaloExchangeForward::launch` rebuilds per-face send lists from
///     current positions, packs owned positions (with PBC shift across
///     periodic-edge faces), issues 6-face `MPI_Isend/Irecv`, and unpacks
///     received ghosts into `frame.atoms[n_local : n_local + n_ghost)`.
///   - `HaloExchangeReverse::launch` symmetrically reduces ghost-row
///     `fx/fy/fz` contributions back onto the owners.
///
/// Both classes share a `HaloExchangeState` record (per-face send/recv
/// indices, packed-position buffers, MPI request handles) held by
/// `DomainDecomposition`. `HaloExchangeForward` is the producer (writes
/// `state->recv_counts_per_face` and the ghost rows); `HaloExchangeReverse`
/// is the consumer (uses the same per-face counts to know how many ghost
/// forces to pack and where they go on the owner's side).
///
/// Default constructor: `cutoff = 0` makes both classes structural no-ops
/// (auto-added in `cpu::md::MD`'s ctor). For real multi-rank runs the user
/// adds an explicit `HaloExchangeForward{cutoff_with_skin}` via
/// `sim.on(pre_force).add(..., -20)`; the no-op auto-instance is harmless
/// (same order, no atoms to send when cutoff is zero).
///
/// **Bonded skin auto-extension** (spec A.5): `build()` scans
/// `frame.bonds/angles/dihedrals` (tag-aware: multi-rank scans resolve
/// atom tags through `dd->tag_to_local` and `MPI_Allreduce(MAX)` across
/// ranks) to find the widest intra-tuple span. `launch()` uses
/// `max(cutoff, bonded_extension)` as the effective halo width -- the
/// bonded extension only *raises* the halo when the user-supplied pair
/// cutoff is smaller than a bonded tuple's reach. Callers no longer have
/// to hand-roll the bonded skin width on top of the pair cutoff.
///
/// @see .claude/specs/parallel-md-engine.md  § A.4, A.5
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

namespace detail {
/// @brief Per-face MPI buffer storage (PIMPL). Definition lives in
/// `halo_exchange.cpp` so the header stays free of MPI-adjacent state.
/// Held per-instance by `HaloExchangeForward` / `HaloExchangeReverse` so
/// multiple components in the same process do not clobber each other's
/// send/recv buffers.
struct FaceBufs;
} // namespace detail

/// @brief Forward pass: owned positions → ghost rows on face neighbors.
///
/// **Real implementation (Phase 2).** Per-face send lists are rebuilt every
/// step from current owned positions; packed buffers are allocated once in
/// `build()` based on a generous initial estimate (no `launch()` realloc on
/// the steady-state path). Single-rank (`runtime.mpi_size == 1`) and
/// `cutoff == 0` both early-return.
///
/// Atoms within `cutoff` of a face plane are flagged as ghosts to the
/// corresponding face neighbor. PBC: atoms sent across the periodic edge of
/// the global box are shifted by `±L_d` on the way out so the receiver sees
/// them in its own coordinate frame.
class HaloExchangeForward : public cpu::ComponentBase {
public:
    /// @brief No-op default constructor (cutoff = 0).
    HaloExchangeForward();

    /// @brief Real-comm constructor with the pair cutoff.
    /// @param cutoff Pair-interaction cutoff [Angstrom]. Typically
    ///        `r_cut + skin` where `r_cut` matches the neighbor list /
    ///        pair potential cutoff. The effective halo width used at
    ///        `launch()` is `max(cutoff, bonded_extension)`, where
    ///        `bonded_extension` is the widest intra-tuple span in
    ///        `frame.bonds/angles/dihedrals` computed once at `build()`.
    ///        When `cutoff >= bonded_extension` the bonded scan has no
    ///        effect (pair cutoff already dominates).
    /// @param every Interval (in steps) between **full** rebuilds of
    ///        per-face send lists and ghost-slot assignments. On steps
    ///        between rebuilds the halo only refreshes ghost positions
    ///        for the existing send lists -- slot indices stay stable
    ///        so the downstream neighbor list remains reusable. Matches
    ///        LAMMPS's reneighboring cadence. Default 10 lines up with
    ///        `NeighborList`'s default `every=10`.
    explicit HaloExchangeForward(double cutoff, int every = 10);

    // unique_ptr<incomplete type> needs user-defined special members,
    // defined in the .cpp where detail::FaceBufs is complete.
    ~HaloExchangeForward() override;
    HaloExchangeForward(HaloExchangeForward&&) noexcept;
    HaloExchangeForward& operator=(HaloExchangeForward&&) noexcept;
    HaloExchangeForward(const HaloExchangeForward&) = delete;
    HaloExchangeForward& operator=(const HaloExchangeForward&) = delete;

    std::string name() const override { return "halo_forward"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::pre_force), -20}};
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Pair cutoff this exchange was constructed with (no bonded skin).
    double cutoff() const { return cutoff_; }

    /// @brief Bonded-skin extension computed by `build()` from the frame's
    /// bonds / angles / dihedrals. Zero before `build()` or when the frame
    /// has no topology.
    double bonded_extension() const { return bonded_extension_; }

    /// @brief Effective halo cutoff: the larger of the pair cutoff and the
    /// bonded skin extension. An atom near the boundary needs its farthest
    /// pair OR bonded partner ghosted, not both summed. `launch()` uses
    /// this value to decide which owned atoms to ship as ghosts. Returns 0
    /// before `build()`.
    double halo_cutoff() const {
        return (cutoff_ > bonded_extension_) ? cutoff_ : bonded_extension_;
    }

    /// @brief Compute the bonded-skin extension from a frame's topology.
    ///
    /// Scans `frame.bonds`, `frame.angles`, `frame.dihedrals` and returns
    /// the largest pairwise atom-to-atom distance found across any tuple.
    /// For an empty topology the result is `0.0`.
    ///
    /// Intended for unit testing and for callers that want the same number
    /// the component will use internally. The component invokes this
    /// helper once, inside `build()`, and caches the result.
    ///
    /// @param frame Frame whose positions and topology to scan.
    /// @return Maximum intra-tuple pairwise distance in Angstroms.
    static double compute_bonded_extension(const cpu::Frame& frame);

    /// @brief Test-only: per-face send-buffer capacity (doubles) after
    /// `build()` has preallocated the instance's own `FaceBufs`. Returns 0
    /// before `build()` or when MPI is not compiled in. Two instances
    /// built from different `dd->capacity` values are expected to report
    /// different capacities; used as a regression check that the PIMPL
    /// buffers are per-instance rather than a process-global static.
    std::size_t send_buf_capacity_for_test() const;

private:
    double cutoff_ = 0.0;
    double bonded_extension_ = 0.0;
    std::unique_ptr<detail::FaceBufs> bufs_;

    // Reneighboring cadence. On `ago_ >= every_` or when migration bumps
    // `dd->neighbor_generation`, the next `launch()` does a full rebuild
    // of per-face send lists and ghost-slot assignments. Otherwise the
    // launch reuses existing send lists and only refreshes ghost x/y/z.
    int every_ = 10;
    int ago_ = 0;
    std::uint64_t last_rebuild_generation_ = 0;
};

/// @brief Reverse pass: ghost-row forces → owner's force arrays.
///
/// Phase 2: pairs with the most recent `HaloExchangeForward` exchange via
/// the shared per-face state on `DomainDecomposition`. Ghost forces in
/// `frame.atoms.fx/fy/fz[n_local : n_local + n_ghost)` are packed by face,
/// `MPI_Isend/Irecv` to the originating ranks, and added into the owners'
/// local force arrays at the indices in their own `send_lists`.
///
/// @see HaloExchangeForward for the matching producer of `dd->send_lists` /
///      `dd->recv_counts`; Reverse consumes those counts to decide how
///      many ghost forces to pack and where to add them on the owner side.
class HaloExchangeReverse : public cpu::ComponentBase {
public:
    HaloExchangeReverse();
    ~HaloExchangeReverse() override;
    HaloExchangeReverse(HaloExchangeReverse&&) noexcept;
    HaloExchangeReverse& operator=(HaloExchangeReverse&&) noexcept;
    HaloExchangeReverse(const HaloExchangeReverse&) = delete;
    HaloExchangeReverse& operator=(const HaloExchangeReverse&) = delete;

    std::string name() const override { return "halo_reverse"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::post_force), -1000}};
    }

    /// @brief Preallocate per-face buffers from `dd->capacity`.
    /// @param ctx Build context.
    void build(const cpu::BuildContext& ctx) override;

    /// @brief Pack ghost-row forces, `MPI_Isend/Irecv` back to owners,
    /// accumulate into owner force arrays.
    /// @param ctx Launch context.
    void launch(const cpu::LaunchContext& ctx) override;

private:
    std::unique_ptr<detail::FaceBufs> bufs_;
};

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
