/// @file cpu/pw/pw_evaluator.cpp
/// @brief PWEvaluator implementation — plane-wave energy (+ FD forces).

#include "atomiverse/cpu/pw/pw_evaluator.h"

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/pw/qe_units.h"

#include <stdexcept>

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

    if (!compute_forces_ || !ctx.forces_requested)
        return;

    const auto forces =
        pw_forces_fd(ctx.frame, solver_.config(), fd_step_); // eV/Å
    const int n = ctx.frame.nlocal();
    if (static_cast<int>(forces.size()) != n)
        throw std::runtime_error("PWEvaluator: force count mismatch");

    auto fx = ctx.frame.span(kAtom_fx);
    auto fy = ctx.frame.span(kAtom_fy);
    auto fz = ctx.frame.span(kAtom_fz);
    for (int i = 0; i < n; ++i) {
        fx[i] = static_cast<AccReal>(forces[static_cast<std::size_t>(i)][0]);
        fy[i] = static_cast<AccReal>(forces[static_cast<std::size_t>(i)][1]);
        fz[i] = static_cast<AccReal>(forces[static_cast<std::size_t>(i)][2]);
    }
}

} // namespace pw
} // namespace cpu
} // namespace ATV
