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
template <typename Strategy = CellList>
class NeighborList : public cpu::ComponentBase {
public:
    /// @param cutoff Interaction cutoff [Angstrom].
    /// @param skin Verlet skin width [Angstrom]. Negative = auto (0.3 * cutoff).
    explicit NeighborList(double cutoff, double skin = -1.0)
        : strategy_(cutoff, skin), cutoff_(cutoff) {}

    std::string name() const override { return "neighbor_list"; }

    void build(const cpu::BuildContext& ctx) override {
        reinit_strategy(ctx.runtime, ctx.frame);
        const int n_active = active_count(ctx.runtime, ctx.frame);
        strategy_.build(ctx.frame, n_active);
        pairs_ = strategy_.build_half(ctx.frame, n_active);
        trim_inum_to_local(ctx.runtime, pairs_);
        ctx.frame.pair_list = &pairs_;
    }

    void launch(const cpu::LaunchContext& ctx) override {
        const int n_active = active_count(ctx.runtime, ctx.frame);
        strategy_.update(ctx.frame, n_active);
        pairs_ = strategy_.build_half(ctx.frame, n_active);
        trim_inum_to_local(ctx.runtime, pairs_);
        ctx.frame.pair_list = &pairs_;
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
};

} // namespace md
} // namespace cpu
} // namespace ATV
