/// @file cpu/potentials/bond.cpp
/// @brief `cpu::Bond::select_owned_tuples` — multi-rank ownership cache.

#include "atomiverse/cpu/potentials/bond.h"
#include "topology_filter.h"

namespace ATV {
namespace cpu {

const std::vector<int>& Bond::select_owned_tuples(
    const LaunchContext& ctx,
    int n_tuples,
    const std::function<std::int64_t(int)>& owner_tag_of) {
    return detail::refresh_owned_tuples(ctx, n_tuples, owner_tag_of,
                                        owned_, cached_generation_);
}

} // namespace cpu
} // namespace ATV
