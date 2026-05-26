/// @file src/cpu/potentials/topology_filter.h
/// @brief Internal helper — multi-rank ownership filter for Bond/Angle/Dihedral.
///
/// Single template that the three bonded bases reuse. Caches the index
/// list against `dd->neighbor_generation` so the filter runs once per
/// migration, not once per step. NOT part of the public API; only the
/// .cpp files in this directory include it.
#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <type_traits>
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

/// @brief Per-step pre-resolved local indices for K-body bonded tuples.
///
/// `tuple_idx[0..n_work)` lists tuple indices that survive both the owner
/// filter (`refresh_owned_tuples`) and per-step local-index resolution
/// against `dd->tag_to_local`. `local[k][0..n_work)` is participant k's
/// local row in `frame.atoms` for every surviving tuple — usable by the
/// hot OMP loop with no `unordered_map` lookups per atom.
template <int K> struct OwnedLocalsView {
    static_assert(K == 2 || K == 3 || K == 4,
                  "K must be 2 (bond), 3 (angle), or 4 (dihedral)");
    const int* tuple_idx;
    std::array<const int*, K> local;
    int n_work;
};

/// @brief Refresh and return a per-step view of resolved local indices.
///
/// Two-phase: (1) migration-cached owner filter; (2) per-step tag→local
/// resolution that ignores tuples whose partners have drifted out of
/// halo skin (silently skipped — see project note
/// `cpu_mpi_halo_edge_corner_bug` for the larger context of why the
/// silent skip is dangerous).
template <int K, typename OwnerTagFn>
inline OwnedLocalsView<K> refresh_resolved_locals(
    const LaunchContext& ctx, int n_tuples,
    const std::array<const int*, K>& atom_arrs, OwnerTagFn owner_tag_of,
    std::vector<int>& owned, std::uint64_t& cached_gen,
    std::vector<int>& step_idx,
    std::array<std::vector<int>, K>& step_locals) {
    std::function<std::int64_t(int)> fn = owner_tag_of;
    refresh_owned_tuples(ctx, n_tuples, fn, owned, cached_gen);

    const int n_owned = static_cast<int>(owned.size());
    step_idx.resize(static_cast<std::size_t>(n_owned));
    for (int k = 0; k < K; ++k)
        step_locals[k].resize(static_cast<std::size_t>(n_owned));

    const auto* dd = ctx.runtime.dd;
    const bool multi_rank = (dd != nullptr) && (ctx.runtime.mpi_size > 1);
    int n_work = 0;
    if (!multi_rank) {
        for (int ib = 0; ib < n_owned; ++ib) {
            const int b = owned[static_cast<std::size_t>(ib)];
            step_idx[static_cast<std::size_t>(n_work)] = b;
            for (int k = 0; k < K; ++k)
                step_locals[k][static_cast<std::size_t>(n_work)]
                    = atom_arrs[k][b];
            ++n_work;
        }
    } else {
        for (int ib = 0; ib < n_owned; ++ib) {
            const int b = owned[static_cast<std::size_t>(ib)];
            std::array<int, K> locals{};
            bool ok = true;
            for (int k = 0; k < K; ++k) {
                const auto it = dd->tag_to_local.find(
                    static_cast<std::int64_t>(atom_arrs[k][b]));
                if (it == dd->tag_to_local.end()) {
                    ok = false;
                    break;
                }
                locals[k] = it->second;
            }
            if (!ok)
                continue;
            step_idx[static_cast<std::size_t>(n_work)] = b;
            for (int k = 0; k < K; ++k)
                step_locals[k][static_cast<std::size_t>(n_work)] = locals[k];
            ++n_work;
        }
    }
    OwnedLocalsView<K> view;
    view.tuple_idx = step_idx.data();
    for (int k = 0; k < K; ++k)
        view.local[k] = step_locals[k].data();
    view.n_work = n_work;
    return view;
}

} // namespace ATV::cpu::detail
