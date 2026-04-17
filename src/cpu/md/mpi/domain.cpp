/// @file src/cpu/md/mpi/domain.cpp
/// @brief `DomainDecomposition::bind` -- single-rank population from a Frame.

#include "atomiverse/cpu/md/mpi/domain.h"

#include "atomiverse/cpu/frame.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

void DomainDecomposition::bind(const Frame& frame) {
    const int n = frame.atoms.size();
    n_local = n;
    n_ghost = 0;
    capacity = n;

    tag_to_local.clear();
    tag_to_local.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        tag_to_local.emplace(frame.atoms.tag(i), i);
    }
}

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
