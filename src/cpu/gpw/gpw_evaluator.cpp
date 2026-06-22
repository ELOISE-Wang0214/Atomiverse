/// @file src/cpu/gpw/gpw_evaluator.cpp
/// @brief Implementation of the GPW energy evaluator component (see header).

#include "atomiverse/cpu/gpw/gpw_evaluator.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/gpw/internal/gpw_forces.h"

namespace ATV {
namespace cpu {
namespace gpw {

namespace {
constexpr double kAng2Bohr = 1.8897259885789;

const GpwElement& find_element(const std::vector<GpwElement>& elements, int z) {
    for (const auto& e : elements)
        if (e.z == z) return e;
    throw std::runtime_error("GPWEvaluator: no basis/potential configured for Z=" +
                             std::to_string(z));
}
} // namespace

void GPWEvaluator::build(const cpu::BuildContext& ctx) {
    // Validate that every atom in the frame has a configured element.
    const auto& a = ctx.frame.atoms;
    for (int i = 0; i < a.size(); ++i) (void)find_element(config_.elements, a.type(i));
}

void GPWEvaluator::launch(const cpu::LaunchContext& ctx) {
    cpu::Frame& frame = ctx.frame;
    const auto& a = frame.atoms;
    const int n = a.size();

    // Assemble atoms at the current geometry (Frame is in Angstrom -> Bohr).
    std::vector<GpwAtom> atoms;
    atoms.reserve(static_cast<std::size_t>(n));
    double z_val_total = 0.0;
    for (int i = 0; i < n; ++i) {
        const GpwElement& e = find_element(config_.elements, a.type(i));
        const double pos[3] = {a.x(i) * kAng2Bohr, a.y(i) * kAng2Bohr, a.z(i) * kAng2Bohr};
        atoms.push_back(make_gpw_atom(pos, e.basis, e.pot, a.type(i)));
        z_val_total += static_cast<double>(e.pot.zion);
    }

    GpwConfig cfg;
    cfg.cell[0] = frame.box.length(0) * kAng2Bohr;
    cfg.cell[1] = frame.box.length(1) * kAng2Bohr;
    cfg.cell[2] = frame.box.length(2) * kAng2Bohr;
    cfg.cutoff_ry = config_.cutoff_ry;
    cfg.n_occ = static_cast<int>(std::lround(z_val_total / 2.0));
    cfg.max_iter = config_.max_iter;
    cfg.tol = config_.tol;
    cfg.mix = config_.mix;
    cfg.xc = config_.xc;
    cfg.nonperiodic = config_.nonperiodic;
    cfg.d3 = config_.d3;
    cfg.d3_param_file = config_.d3_param_file;
    if (std::getenv("ATV_GPW_OT")) cfg.use_ot = true; // diagnostic: force OT path

    last_ = gpw_scf(atoms, cfg);
    frame.scalars.energy = last_.e_total * ATV::kHartree2eV;

    if (config_.compute_forces) {
        const auto F = gpw_forces_analytic(atoms, cfg); // [Hartree/Bohr]; one SCF + gradient (was 6N+1 FD)
        for (int i = 0; i < n; ++i) {
            frame.atoms.fx(i) = F[static_cast<std::size_t>(i)][0] * ATV::kForceAU2Metal;
            frame.atoms.fy(i) = F[static_cast<std::size_t>(i)][1] * ATV::kForceAU2Metal;
            frame.atoms.fz(i) = F[static_cast<std::size_t>(i)][2] * ATV::kForceAU2Metal;
        }
    }
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
