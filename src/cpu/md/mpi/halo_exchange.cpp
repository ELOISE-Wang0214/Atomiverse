/// @file src/cpu/md/mpi/halo_exchange.cpp
/// @brief Phase 1 no-op implementations of `HaloExchangeForward`/`HaloExchangeReverse`.

#include "atomiverse/cpu/md/mpi/halo_exchange.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

void HaloExchangeForward::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size == 1) return;
    // Phase 2: non-blocking position send/recv over 6 face neighbors.
}

void HaloExchangeReverse::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size == 1) return;
    // Phase 2: reduce ghost-row forces back onto owners.
}

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
