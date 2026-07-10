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
        if (e.z == z)
            return e;
    throw std::runtime_error(
        "GPWEvaluator: no basis/potential configured for Z="
        + std::to_string(z));
}
} // namespace

namespace {
/// Assemble the GpwAtom set at the frame's current geometry (Å → Bohr).
std::vector<GpwAtom> assemble_atoms(const cpu::Frame& frame,
                                    const GpwEvaluatorConfig& config) {
    const auto& a = frame.atoms;
    const int n = a.size();
    std::vector<GpwAtom> atoms;
    atoms.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const GpwElement& e = find_element(config.elements, a.type(i));
        const double pos[3] = {a.x(i) * kAng2Bohr, a.y(i) * kAng2Bohr,
                               a.z(i) * kAng2Bohr};
        atoms.push_back(make_gpw_atom(pos, e.basis, e.pot, a.type(i)));
    }
    return atoms;
}
} // namespace

void GPWEvaluator::build(const cpu::BuildContext& ctx) {
    // Validate that every atom in the frame has a configured element, and
    // count valence electrons for the closed-shell occupation.
    const auto& a = ctx.frame.atoms;
    double z_val_total = 0.0;
    for (int i = 0; i < a.size(); ++i)
        z_val_total += static_cast<double>(
            find_element(config_.elements, a.type(i)).pot.zion);

    // Assemble the solver configuration from the build-time frame. The
    // ATV_GPW_OT diagnostic switch is read once here — re-reading an
    // environment variable per launch had no meaning.
    GpwConfig cfg;
    cfg.cell[0] = ctx.frame.box.length(0) * kAng2Bohr;
    cfg.cell[1] = ctx.frame.box.length(1) * kAng2Bohr;
    cfg.cell[2] = ctx.frame.box.length(2) * kAng2Bohr;
    cfg.cutoff_ry = config_.cutoff_ry;
    cfg.n_occ = static_cast<int>(std::lround(z_val_total / 2.0));
    cfg.max_iter = config_.max_iter;
    cfg.tol = config_.tol;
    cfg.mix = config_.mix;
    cfg.xc = config_.xc;
    cfg.nonperiodic = config_.nonperiodic;
    cfg.d3 = config_.d3;
    cfg.d3_param_file = config_.d3_param_file;
    if (std::getenv("ATV_GPW_OT"))
        cfg.use_ot = true; // diagnostic: force OT path
    built_cell_ = {cfg.cell[0], cfg.cell[1], cfg.cell[2]};

    solver_ = std::make_unique<GPWSolver>(std::move(cfg));
    solver_->build(assemble_atoms(ctx.frame, config_));
}

void GPWEvaluator::launch(const cpu::LaunchContext& ctx) {
    cpu::Frame& frame = ctx.frame;
    const int n = frame.atoms.size();

    // GPW does not support a variable cell: the grids, lattice images and
    // Ewald-like electrostatics are sized from the build-time cell. The old
    // implementation silently solved in the new box — undefined behavior.
    for (int d = 0; d < 3; ++d) {
        const double now = frame.box.length(d) * kAng2Bohr;
        if (std::abs(now - built_cell_[static_cast<std::size_t>(d)])
            > 1e-9
                  * std::max(
                      1.0,
                      std::abs(built_cell_[static_cast<std::size_t>(d)]))) {
            throw std::runtime_error(
                "GPWEvaluator::launch: frame box changed since build() "
                "(GPW does not support a variable cell)");
        }
    }

    const auto atoms = assemble_atoms(frame, config_);
    last_ = solver_->run(atoms);
    frame.scalars.energy = last_.e_total * ATV::kHartree2eV;

    if (config_.compute_forces) {
        const auto F = gpw_forces_analytic(atoms, solver_->config());
        // [Hartree/Bohr]; one SCF + gradient (was 6N+1 FD)
        for (int i = 0; i < n; ++i) {
            frame.atoms.fx(i) =
                F[static_cast<std::size_t>(i)][0] * ATV::kForceAU2Metal;
            frame.atoms.fy(i) =
                F[static_cast<std::size_t>(i)][1] * ATV::kForceAU2Metal;
            frame.atoms.fz(i) =
                F[static_cast<std::size_t>(i)][2] * ATV::kForceAU2Metal;
        }
    }
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
