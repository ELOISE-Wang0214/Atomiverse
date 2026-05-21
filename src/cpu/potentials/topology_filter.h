/// @file src/cpu/potentials/topology_filter.h
/// @brief Internal helper — multi-rank ownership filter for Bond/Angle/Dihedral.
///
/// Single template that the three bonded bases reuse. Caches the index
/// list against `dd->neighbor_generation` so the filter runs once per
/// migration, not once per step. NOT part of the public API; only the
/// .cpp files in this directory include it.
#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/md/mpi/domain.h"

namespace ATV::cpu::detail {

/// @brief Refresh `owned` to the index list of tuples this rank owns.
/// @param ctx              LaunchContext (provides runtime + DD pointer).
/// @param n_tuples         Total tuples in the frame's topology block.
/// @param owner_tag_of     Callback `idx → global tag of owner atom`.
/// @param owned            Cache vector owned by the calling base class.
/// @param cached_gen       Cache generation owned by the calling base.
/// @return Const reference to `owned`.
inline const std::vector<int>&
refresh_owned_tuples(const LaunchContext& ctx, int n_tuples,
                     const std::function<std::int64_t(int)>& owner_tag_of,
                     std::vector<int>& owned, std::uint64_t& cached_gen) {
    const auto* dd = ctx.runtime.dd;
    const bool multi_rank = (dd != nullptr) && (ctx.runtime.mpi_size > 1);
    if (!multi_rank) {
        if (static_cast<int>(owned.size()) != n_tuples) {
            owned.resize(static_cast<std::size_t>(n_tuples));
            for (int t = 0; t < n_tuples; ++t)
                owned[t] = t;
        }
        return owned;
    }
    if (cached_gen == dd->neighbor_generation)
        return owned;
    owned.clear();
    owned.reserve(static_cast<std::size_t>(n_tuples));
    for (int t = 0; t < n_tuples; ++t) {
        const std::int64_t owner_tag = owner_tag_of(t);
        auto it = dd->tag_to_local.find(owner_tag);
        if (it != dd->tag_to_local.end() && it->second < dd->n_local) {
            owned.push_back(t);
        }
    }
    cached_gen = dd->neighbor_generation;
    return owned;
}

} // namespace ATV::cpu::detail
