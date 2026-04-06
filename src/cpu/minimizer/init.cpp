/// @file cpu/minimizer/init.cpp
/// @brief Minimizer infrastructure component implementations.

#include "atomiverse/cpu/minimizer/init.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace minimizer {

void ZeroForces::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    Int n = f.atoms.size();
    #pragma omp parallel for schedule(static)
    for (Int i = 0; i < n; ++i) {
        f.atoms.fx[i] = 0;
        f.atoms.fy[i] = 0;
        f.atoms.fz[i] = 0;
    }
}

void ZeroEnergy::launch(const cpu::LaunchContext& ctx) {
    ctx.frame.energy = 0;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
