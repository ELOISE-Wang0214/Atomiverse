/// @file cpu/potentials/special_scale.h
/// @brief Shared LAMMPS special_bonds scaling table for pair styles.
///
/// Builds a tag→(partner, scale) map from `kAtom_nspecial` / `kAtom_special`.
/// Scales multiply pair energy and force (0 = skip, 0.5 = Amber 1–4, 1 = full).
#pragma once

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Per-owner special-bond weight lookup (LAMMPS special_bonds).
struct SpecialScaleTable {
    double s12 = 1.0;
    double s13 = 1.0;
    double s14 = 1.0;

    /// tag → sorted list of (partner_tag, scale) for scale != 1
    std::unordered_map<std::int64_t, std::vector<std::pair<std::int64_t, double>>>
        by_tag;

    void clear() { by_tag.clear(); }

    /// @return true if any of s12/s13/s14 is not 1 (specials active).
    bool active() const {
        return s12 != 1.0 || s13 != 1.0 || s14 != 1.0;
    }

    /// @brief Rebuild from current topology of owner rows [0, n_owners).
    void rebuild(const Frame& f, int n_owners) {
        by_tag.clear();
        if (!active())
            return;

        auto add = [&](std::int64_t a, std::int64_t b, double s) {
            if (a == b || s == 1.0)
                return;
            by_tag[a].push_back({b, s});
            by_tag[b].push_back({a, s});
        };

        const auto nspecial = f.span(kAtom_nspecial);
        const auto special = f.span(kAtom_special);
        const TagInt* tag = f.span(kAtom_tag);
        for (int i = 0; i < n_owners; ++i) {
            const std::int64_t ti = tag[i];
            const int n12 = nspecial(i, 0);
            const int n13 = nspecial(i, 1);
            const int n14 = nspecial(i, 2);
            for (int s = 0; s < n12; ++s)
                add(ti, special(i, s), s12);
            for (int s = n12; s < n13; ++s)
                add(ti, special(i, s), s13);
            for (int s = n13; s < n14; ++s)
                add(ti, special(i, s), s14);
        }
        for (auto& kv : by_tag) {
            auto& list = kv.second;
            std::sort(list.begin(), list.end(),
                      [](const auto& a, const auto& b) {
                          return a.first < b.first;
                      });
            // Prefer nearer tier if duplicate partners (keep first = smaller
            // tier already added first). Unique by partner keeping first scale.
            std::vector<std::pair<std::int64_t, double>> uniq;
            uniq.reserve(list.size());
            for (const auto& p : list) {
                if (uniq.empty() || uniq.back().first != p.first)
                    uniq.push_back(p);
            }
            list.swap(uniq);
        }
    }

    /// @brief Scale for pair (ti,tj); default 1.0 when not special.
    double scale(std::int64_t ti, std::int64_t tj) const {
        auto it = by_tag.find(ti);
        if (it == by_tag.end())
            return 1.0;
        const auto& list = it->second;
        auto lb = std::lower_bound(
            list.begin(), list.end(), tj,
            [](const std::pair<std::int64_t, double>& p, std::int64_t t) {
                return p.first < t;
            });
        if (lb != list.end() && lb->first == tj)
            return lb->second;
        return 1.0;
    }
};

} // namespace md
} // namespace cpu
} // namespace ATV
