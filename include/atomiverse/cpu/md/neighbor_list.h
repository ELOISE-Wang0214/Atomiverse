/// @file cpu/md/neighbor_list.h
/// @brief NeighborList -- component wrapper for spatial index strategies.
///
/// @tparam Strategy Spatial indexing strategy (CellList, NSQ, etc.).
#pragma once

#include <string>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/md/pair_list.h"
#include "atomiverse/cpu/md/cell_list.h"
#include "atomiverse/cpu/md/mpi/domain.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Neighbor list component -- builds and maintains spatial index.
/// @tparam Strategy Spatial indexing strategy (CellList by default).
///
/// Rebuild cadence mirrors LAMMPS `neigh_modify every N delay 0 check yes`:
/// the Verlet list is re-enumerated only when `ago >= every` *and* the
/// max-displacement skin check fires. Between rebuilds the stored CSR is
/// reused verbatim — force kernels continue reading the same `pair_list`.
/// Cost of `launch()` on non-rebuild steps collapses to a single O(n) skin
/// scan (or zero when `ago < every`). Default `every = 1` preserves the
/// LAMMPS default.
template <typename Strategy = CellList>
class NeighborList : public cpu::ComponentBase {
public:
    /// @param cutoff Interaction cutoff [Angstrom].
    /// @param skin Verlet skin width [Angstrom]. Negative = auto (0.3 * cutoff).
    /// @param every Attempt rebuild at most every N steps (LAMMPS `neigh_modify
    ///        every`; default 1 reproduces LAMMPS's built-in default).
    explicit NeighborList(double cutoff, double skin = -1.0, int every = 1)
        : strategy_(cutoff, skin), cutoff_(cutoff),
          every_(every < 1 ? 1 : every) {}

    std::string name() const override { return "neighbor_list"; }

    void build(const cpu::BuildContext& ctx) override {
        reinit_strategy(ctx.runtime, ctx.frame);
        const int n_active = active_count(ctx.runtime, ctx.frame);
        strategy_.build(ctx.frame, n_active);
        // One-time arena sizing. The xtensor CSR arrays grow (never shrink)
        // through `refresh_half`'s `ensure_at_least`; prime them here to
        // the expected upper bound so steady-state rebuilds don't enter
        // the allocator. `pairs_.neighbors` is a `std::vector<int>` and
        // amortises growth itself -- leave it empty, the first refresh
        // reaches the correct size within one resize.
        strategy_.refresh_half(ctx.frame, n_active, pairs_);
        ++pairs_.version;
        trim_inum_to_local(ctx.runtime, pairs_);
        ctx.runtime.pair_list = &pairs_;
        ago_ = 0;
    }

    void launch(const cpu::LaunchContext& ctx) override {
        const int n_active = active_count(ctx.runtime, ctx.frame);
        ++ago_;

        // `ago < every`: reuse the list verbatim when the atom count is
        // stable. An upstream `AtomMigration` or halo-size change, however,
        // makes the current CSR size-incompatible with the trimmed `inum`
        // that force kernels will walk — so we must rebuild as soon as the
        // CSR size is stale, regardless of `ago`. The `every` cadence still
        // gates the displacement-driven rebuild path.
        const bool size_stale =
            n_active != strategy_.n_at_last_build();
        const bool disp_rebuild =
            ago_ >= every_ && strategy_.needs_rebuild(ctx.frame, n_active);
        if (size_stale || disp_rebuild) {
            strategy_.build(ctx.frame, n_active);
            strategy_.refresh_half(ctx.frame, n_active, pairs_);
            ++pairs_.version;
            ago_ = 0;
        }

        trim_inum_to_local(ctx.runtime, pairs_);
        // Re-publish defensively: the component that owns pairs_ is stable
        // across launches, but a subsequent Backend::build_components() cycle
        // (e.g. after a frame-resize) would reset the runtime pointer.
        ctx.runtime.pair_list = &pairs_;
    }

    /// @brief Access the current pair list (read-only).
    const PairList& pairs() const { return pairs_; }

    /// @brief Access the underlying strategy.
    const Strategy& strategy() const { return strategy_; }

private:
    /// @brief True when running with multiple MPI ranks (and the DD struct
    /// is available). Single-rank or null-`dd` returns false.
    static bool is_multi_rank(const cpu::Context& rt) {
        return rt.dd != nullptr && rt.mpi_size > 1;
    }

    /// @brief Active atom count: owners + current ghosts under multi-rank,
    /// otherwise the full atom array (single-rank legacy).
    static int active_count(const cpu::Context& rt, const cpu::Frame& frame) {
        if (is_multi_rank(rt)) return rt.dd->n_local + rt.dd->n_ghost;
        return frame.atoms.size();
    }

    /// @brief Multi-rank trim: `inum = dd->n_local` so force kernels only
    /// iterate owned atoms as `i`, while `neighbors[]` may still point to
    /// ghost rows. Single-rank: no-op.
    static void trim_inum_to_local(const cpu::Context& rt, PairList& pl) {
        if (is_multi_rank(rt)) {
            pl.inum = rt.dd->n_local;
        }
    }

    /// @brief (Re-)initialise the underlying spatial-index strategy from
    /// the current frame.box.  Called from `build()` (one-time setup)
    /// and from `launch()` whenever `frame.box` has drifted enough that
    /// the cell-grid sizing is stale.
    void reinit_strategy(const cpu::Context& rt, cpu::Frame& f) {
        if (is_multi_rank(rt)) {
            strategy_.init_extended(f, cutoff_, rt.dd->dims, rt.dd->capacity);
        } else {
            strategy_.init(f, cutoff_);
        }
    }

    Strategy strategy_;
    PairList pairs_;
    double cutoff_;
    int every_ = 1;
    /// @brief Steps since the last actual rebuild. Starts at 0 after `build()`
    /// (the initial enumeration), incremented once per `launch()`.
    int ago_ = 0;
};

} // namespace md
} // namespace cpu
} // namespace ATV
