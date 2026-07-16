/// @file cpu/pw/pw_evaluator.cpp
/// @brief PWEvaluator implementation — plane-wave energy component.
#include "atomiverse/cpu/pw/pw_evaluator.h"

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/pw/qe_units.h"

namespace ATV {
namespace cpu {
namespace pw {

void PWEvaluator::build(const cpu::BuildContext& ctx) {
    solver_.build(ctx.frame);
}

void PWEvaluator::launch(const cpu::LaunchContext& ctx) {
    const auto& result = solver_.run(ctx.frame);
    // Project units convention: frame-level energies are eV.
    ctx.frame.meta.get(kEnergy) = result.energy * qe_units::rydberg_to_ev;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
