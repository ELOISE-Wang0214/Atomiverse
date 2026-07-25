/// @file src/cpu/gpw/gpw_evaluator.cpp
/// @brief Implementation of the GPW energy evaluator component (see header).

#include "atomiverse/cpu/gpw/gpw_evaluator.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <xtensor/views/xview.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/gpw/internal/gpw_forces.h"
#include "atomiverse/cpu/gpw/internal/gpw_stress.h"
#include "atomiverse/cpu/linalg/density.h"

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
        z_val_total +=
            static_cast<double>(find_kind(kinds, sp_type_build[i]).pot.zion);

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

    // XL-BOMD guards (no-fallback): the auxiliary-density guess feeds the
    // closed-shell RKS diagonalization path only.
    if (xl_bomd) {
        if (n_alpha + n_beta > 0)
            throw std::invalid_argument(
                "GPWEvaluator: xl_bomd supports closed-shell RKS only "
                "(n_alpha/n_beta are set)");
        if (cfg.use_ot)
            throw std::invalid_argument(
                "GPWEvaluator: xl_bomd requires the diagonalization SCF path "
                "(use_ot parameterizes orbitals, not the density)");
        if (cpmd)
            throw std::invalid_argument(
                "GPWEvaluator: xl_bomd and cpmd are mutually exclusive "
                "(XL is SCF-converged auxiliary density; CP is fictitious "
                "orbital dynamics)");
        xl_.emplace(xl_bomd_k); // validates xl_bomd_k in [3, 9]
    } else {
        xl_.reset();
    }

    // Car–Parrinello guards: closed-shell RKS Γ diagonalization path only.
    if (cpmd) {
        if (n_alpha + n_beta > 0)
            throw std::invalid_argument(
                "GPWEvaluator: cpmd supports closed-shell RKS only "
                "(n_alpha/n_beta are set)");
        if (cfg.use_ot)
            throw std::invalid_argument(
                "GPWEvaluator: cpmd requires the diagonalization SCF path for "
                "the cold seed (use_ot is refused — no silent OT-CP hybrid)");
        if (!(cpmd_mu > 0.0) || !std::isfinite(cpmd_mu))
            throw std::invalid_argument(
                "GPWEvaluator: cpmd_mu must be finite and > 0 (a.u. mass)");
        if (!(cpmd_dt > 0.0) || !std::isfinite(cpmd_dt))
            throw std::invalid_argument(
                "GPWEvaluator: cpmd_dt must be finite and > 0 (fs); must match "
                "the nuclear VelocityVerlet timestep");
        if (!compute_forces)
            throw std::invalid_argument(
                "GPWEvaluator: cpmd requires compute_forces=true (nuclear "
                "forces from the fixed-density electronic state)");
        cpmd_.emplace(cpmd_mu);
        last_te_ev_ = 0.0;
    } else {
        cpmd_.reset();
        last_te_ev_ = 0.0;
    }

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

    // ── Car–Parrinello path ──────────────────────────────────────────────
    // Seed: one full SCF + BO analytic forces, then hand occupied C to Cpmd.
    // Later steps: fixed-density KS eval (max_iter=0) + forces at that
    // density + electronic Verlet. Nuclear PE is E_KS only; T_e is exposed
    // via last_electronic_ke_ev() for the CP conserved energy.
    if (cpmd_ && cpmd_->ready()) {
        auto fixed_density_eval = [&](const linalg::RealMatrix& C_occ) {
            const int nb = static_cast<int>(C_occ.shape(0));
            const int n_occ = static_cast<int>(C_occ.shape(1));
            const linalg::RealMatrix P =
                linalg::build_density_matrix(C_occ, n_occ, nb, 2.0);
            const int saved_max = solver_->max_iter;
            solver_->max_iter = 0;
            try {
                last_ = solver_->run(atoms, /*guess_C=*/nullptr, &P);
            } catch (...) {
                solver_->max_iter = saved_max;
                throw;
            }
            solver_->max_iter = saved_max;
        };

        // Pass 1: KS eval at transported C to obtain the new geometry's S.
        fixed_density_eval(cpmd_->C());
        // Nuclei moved ⇒ S changed; re-orthonormalize electronic history.
        cpmd_->align_to_overlap(last_.S);
        // Pass 2: energy / H / forces at the S-orthonormal density.
        fixed_density_eval(cpmd_->C());

        frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

        if (compute_forces) {
            const auto F = gpw_forces_at_density(atoms, *solver_, last_.P,
                                                 cpmd_->C(), last_.H_ks);
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

        // Electronic step for the *next* nuclear force evaluation.
        cpmd_->advance(last_.H_ks, last_.S, cpmd_dt);
        last_te_ev_ = cpmd_->kinetic_ha() * ATV::units::kHartreeInEv;
        return;
    }

    // BOMD warm-start: reuse previous MO coefficients when shape matches.
    // GpwDft::run accepts optional guess_C; empty last_.C → core guess.
    const ATV::cpu::linalg::RealMatrix* guess = nullptr;
    if (last_.C.size() > 0) {
        guess = &last_.C;
    }
    // XL-BOMD: the propagated auxiliary density is the SCF initial-density
    // guess (Niklasson dissipative recursion, seeded on the first step).
    const ATV::cpu::linalg::RealMatrix* guess_p =
        (xl_ && xl_->ready()) ? &xl_->aux() : nullptr;
    last_ = solver_->run(atoms, guess, guess_p);
    if (xl_) {
        if (last_.converged && last_.P.size() > 0)
            xl_->advance(last_.P);
        else
            xl_->reset(); // discontinuity: never propagate an unconverged D
    }
    frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

    if (compute_forces) {
        // XL mode warm-starts the force-pass SCF from the MO coefficients
        // just converged at this exact geometry (same fixed point, fewer
        // iterations). Default mode keeps the cold solve unchanged.
        const auto F = gpw_forces_analytic(
            atoms, *solver_, (xl_ && last_.C.size() > 0) ? &last_.C : nullptr);
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

    // CPMD cold seed after the first full SCF + BO forces.
    if (cpmd_ && !cpmd_->ready()) {
        if (!last_.converged || last_.C.size() == 0 || solver_->n_occ < 1)
            throw std::runtime_error(
                "GPWEvaluator: cpmd seed requires a converged RKS SCF with "
                "occupied MO coefficients");
        const int nb = static_cast<int>(last_.C.shape(0));
        const int n_occ = solver_->n_occ;
        linalg::RealMatrix C_occ = xt::zeros<double>(
            {static_cast<std::size_t>(nb), static_cast<std::size_t>(n_occ)});
        for (int mu = 0; mu < nb; ++mu)
            for (int i = 0; i < n_occ; ++i)
                C_occ(mu, i) = last_.C(mu, i);
        cpmd_->seed(C_occ);
        last_te_ev_ = 0.0;
    }

    if (compute_stress) {
        // Product default: Hellmann–Feynman fixed-P analytic stress (one SCF
        // already done above + fixed-P energy samples). Optional multi-SCF
        // strain-FD oracle via stress_use_fd. Warm-start from last MOs; skip
        // masked axes (variable-cell relax of constrained boxes).
        const ATV::cpu::linalg::RealMatrix* stress_guess =
            last_.C.size() > 0 ? &last_.C : nullptr;
        last_stress_ha_ = stress_use_fd
                              ? gpw_stress_fd(atoms, *solver_, stress_strain_h,
                                              stress_axes, stress_guess)
                              : gpw_stress(atoms, *solver_, stress_strain_h,
                                           stress_axes, stress_guess);
        // 1 Ha/Bohr³ = kHartreeInEv / (kBohrInAngstrom)³ eV/Å³
        const double ha_bohr3_to_ev_a3 =
            ATV::units::kHartreeInEv
            / (ATV::units::kBohrInAngstrom * ATV::units::kBohrInAngstrom
               * ATV::units::kBohrInAngstrom);
        auto& st = frame.meta.get(kStress);
        for (int k = 0; k < 6; ++k)
            st[static_cast<std::size_t>(k)] =
                last_stress_ha_[static_cast<std::size_t>(k)]
                * ha_bohr3_to_ev_a3;
    }
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
