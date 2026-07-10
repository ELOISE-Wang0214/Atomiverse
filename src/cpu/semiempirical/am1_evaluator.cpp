/// @file src/cpu/semiempirical/am1_evaluator.cpp
/// @brief AM1 evaluator implementation.

#include "atomiverse/cpu/semiempirical/am1_evaluator.h"

#include "atomiverse/core/units.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

void AM1Evaluator::build(const cpu::BuildContext& ctx) {
    solver_.build(ctx.frame);
}

void AM1Evaluator::launch(const cpu::LaunchContext& ctx) {
    auto& result = solver_.run(ctx.frame);
    ctx.frame.scalars.energy = result.electronic_energy_hartree * kHartree2eV;
    const int n = ctx.frame.atoms.size();
    for (int i = 0; i < n; ++i) {
        ctx.frame.atoms.charge(i) =
            result.mulliken_charges[static_cast<std::size_t>(i)];
    }
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
