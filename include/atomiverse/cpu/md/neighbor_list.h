/// @file cpu/md/neighbor_list.h
/// @brief NeighborList -- component wrapper for spatial index strategies.
///
/// @tparam Strategy Spatial indexing strategy (CellList, NSQ, etc.).
#pragma once

#include <string>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/md/pair_list.h"
#include "atomiverse/cpu/md/cell_list.h"

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
        strategy_.init(ctx.frame, cutoff_);
        strategy_.build(ctx.frame);
        pairs_ = strategy_.build_half(ctx.frame);
        ctx.frame.pair_list = &pairs_;
    }

    void launch(const cpu::LaunchContext& ctx) override {
        strategy_.update(ctx.frame);
        pairs_ = strategy_.build_half(ctx.frame);
        ctx.frame.pair_list = &pairs_;
    }

    /// @brief Access the current pair list (read-only).
    const PairList& pairs() const { return pairs_; }

    /// @brief Access the underlying strategy.
    const Strategy& strategy() const { return strategy_; }

private:
    Strategy strategy_;
    PairList pairs_;
    double cutoff_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
