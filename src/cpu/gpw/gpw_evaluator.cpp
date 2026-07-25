/// @file src/cpu/gpw/gpw_evaluator.cpp
/// @brief Implementation of the GPW energy evaluator component (see header).

#include "atomiverse/cpu/gpw/gpw_evaluator.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/gpw/internal/gpw_forces.h"
#include "atomiverse/cpu/gpw/internal/gpw_stress.h"

namespace ATV {
namespace cpu {
namespace gpw {

namespace {
constexpr double kAng2Bohr = 1.8897259885789;

const Kind& find_kind(const std::vector<Kind>& kinds, int z) {
    for (const auto& k : kinds)
        if (static_cast<int>(k.element) == z)
            return k;
    throw std::runtime_error(
        "GPWEvaluator: no basis/potential configured for Z="
        + std::to_string(z));
}
} // namespace

namespace {
/// Assemble the GpwAtom set at the frame's current geometry (Å → Bohr).
std::vector<GpwAtom> assemble_atoms(const cpu::Frame& frame,
                                    const GPWEvaluator& self) {
    const int n = frame.nlocal();
    const int* __restrict sp_type = frame.span(kAtom_type);
    const StorageReal* __restrict sp_x = frame.span(kAtom_x);
    const StorageReal* __restrict sp_y = frame.span(kAtom_y);
    const StorageReal* __restrict sp_z = frame.span(kAtom_z);
    std::vector<GpwAtom> atoms;
    atoms.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const Kind& k = find_kind(self.kinds, sp_type[i]);
        const double pos[3] = {sp_x[i] * kAng2Bohr, sp_y[i] * kAng2Bohr,
                               sp_z[i] * kAng2Bohr};
        atoms.push_back(GpwAtom(pos, k.basis, k.pot, sp_type[i]));
    }
    return atoms;
}
} // namespace

void GPWEvaluator::build(const cpu::BuildContext& ctx) {
    // Validate that every atom in the frame has a configured element, and
    // count valence electrons for the closed-shell occupation.
    const int n_build = ctx.frame.nlocal();
    const int* __restrict sp_type_build = ctx.frame.span(kAtom_type);
    double z_val_total = 0.0;
    for (int i = 0; i < n_build; ++i)
        z_val_total += static_cast<double>(
            find_kind(kinds, sp_type_build[i]).pot.zion);

    // Assemble the solver from the build-time frame. Product path is the
    // `use_ot` member; env `ATV_GPW_OT` remains a diagnostic force-on.
    GpwDft cfg;
    cfg.cell[0] = ctx.frame.simbox.length(0) * kAng2Bohr;
    cfg.cell[1] = ctx.frame.simbox.length(1) * kAng2Bohr;
    cfg.cell[2] = ctx.frame.simbox.length(2) * kAng2Bohr;
    cfg.cutoff_ry = cutoff_ry;
    cfg.n_alpha = n_alpha;
    cfg.n_beta = n_beta;
    if (n_alpha + n_beta > 0) {
        // UKS: occupations are explicit; n_occ unused.
        cfg.n_occ = 0;
    } else {
        cfg.n_occ = static_cast<int>(std::lround(z_val_total / 2.0));
    }
    cfg.max_iter = max_iter;
    cfg.tol = tol;
    cfg.mix = mix;
    cfg.xc = xc;
    cfg.nonperiodic = nonperiodic;
    cfg.d3 = d3;
    cfg.d3bj = d3bj;
    cfg.d4 = d4;
    cfg.use_ot = use_ot || (std::getenv("ATV_GPW_OT") != nullptr);
    cfg.plus_u = plus_u;
    cfg.atom_spin_guess = atom_spin_guess;

    solver_ = std::make_unique<GpwDft>(std::move(cfg));
    solver_->build(assemble_atoms(ctx.frame, *this));
}

void GPWEvaluator::launch(const cpu::LaunchContext& ctx) {
    cpu::Frame& frame = ctx.frame;
    const int n = frame.nlocal();

    // Variable cell: re-target the solver to the frame's current orthorhombic
    // box. GpwDft::run rebuilds every geometry-dependent structure (grids,
    // lattice images, S/T/Vpp) from `cell` per call, so a box change between
    // launches is well-defined — cell-relax drivers rescale simbox freely.
    for (int d = 0; d < 3; ++d) {
        const double now = frame.simbox.length(d) * kAng2Bohr;
        if (!(now > 0.0))
            throw std::runtime_error(
                "GPWEvaluator::launch: non-positive box edge");
        solver_->cell[d] = now;
    }

    const auto atoms = assemble_atoms(frame, *this);
    // BOMD warm-start: reuse previous MO coefficients when shape matches.
    // GpwDft::run accepts optional guess_C; empty last_.C → core guess.
    const ATV::cpu::linalg::RealMatrix* guess = nullptr;
    if (last_.C.size() > 0) {
        guess = &last_.C;
    }
    last_ = solver_->run(atoms, guess);
    frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

    if (compute_forces) {
        const auto F = gpw_forces_analytic(atoms, *solver_);
        // [Hartree/Bohr]; one SCF + gradient (was 6N+1 FD)
        AccReal* __restrict sp_fx = frame.span(kAtom_fx);
        AccReal* __restrict sp_fy = frame.span(kAtom_fy);
        AccReal* __restrict sp_fz = frame.span(kAtom_fz);
        for (int i = 0; i < n; ++i) {
            sp_fx[i] = F[static_cast<std::size_t>(i)][0]
                       * ATV::units::kHartreePerBohrInEvPerAngstrom;
            sp_fy[i] = F[static_cast<std::size_t>(i)][1]
                       * ATV::units::kHartreePerBohrInEvPerAngstrom;
            sp_fz[i] = F[static_cast<std::size_t>(i)][2]
                       * ATV::units::kHartreePerBohrInEvPerAngstrom;
        }
    }

    if (compute_stress) {
        // Strain-FD rebuilds SCF per strain sample (independent of solver_).
        // Warm-start every strained SCF from the just-converged MOs; skip
        // masked axes entirely (variable-cell relax of constrained boxes).
        last_stress_ha_ =
            gpw_stress(atoms, *solver_, stress_strain_h, stress_axes,
                       last_.C.size() > 0 ? &last_.C : nullptr);
        // 1 Ha/Bohr³ = kHartreeInEv / (kBohrInAngstrom)³ eV/Å³
        const double ha_bohr3_to_ev_a3 =
            ATV::units::kHartreeInEv
            / (ATV::units::kBohrInAngstrom * ATV::units::kBohrInAngstrom
               * ATV::units::kBohrInAngstrom);
        auto& st = frame.meta.get(kStress);
        for (int k = 0; k < 6; ++k)
            st[static_cast<std::size_t>(k)] =
                last_stress_ha_[static_cast<std::size_t>(k)] * ha_bohr3_to_ev_a3;
    }
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
