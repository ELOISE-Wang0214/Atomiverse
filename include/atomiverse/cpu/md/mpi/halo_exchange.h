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
/// @see .claude/specs/parallel-md-engine.md  § A.4
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

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
    HaloExchangeForward() = default;

    /// @brief Real-comm constructor with the ghost cutoff.
    /// @param cutoff Halo cutoff [Angstrom]. Typically `r_cut + skin` where
    ///        `r_cut` matches the neighbor list / pair potential cutoff.
    explicit HaloExchangeForward(double cutoff) : cutoff_(cutoff) {}

    std::string name() const override { return "halo_forward"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::pre_force), -20}};
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Cutoff this exchange was constructed with (for diagnostics).
    double cutoff() const { return cutoff_; }

private:
    double cutoff_ = 0.0;
};

/// @brief Reverse pass: ghost-row forces → owner's force arrays.
///
/// Phase 2: pairs with the most recent `HaloExchangeForward` exchange via
/// the shared per-face state on `DomainDecomposition`. Ghost forces in
/// `frame.atoms.fx/fy/fz[n_local : n_local + n_ghost)` are packed by face,
/// `MPI_Isend/Irecv` to the originating ranks, and added into the owners'
/// local force arrays at the indices in their own `send_lists`.
class HaloExchangeReverse : public cpu::ComponentBase {
public:
    std::string name() const override { return "halo_reverse"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::post_force), -1000}};
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
