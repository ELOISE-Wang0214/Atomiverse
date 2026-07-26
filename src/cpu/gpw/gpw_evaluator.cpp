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
    cfg.ts = ts;
    cfg.xdm = xdm;
    cfg.mbd = mbd;
    for (int d = 0; d < 3; ++d)
        cfg.mbd_k_mesh[d] = mbd_k_mesh[d];
    cfg.vdw_df = vdw_df;
    cfg.vdw_df2 = vdw_df2;
    cfg.rvv10 = rvv10;
    cfg.hirshfeld = hirshfeld;
    cfg.use_ot = use_ot || (std::getenv("ATV_GPW_OT") != nullptr);
    cfg.plus_u = plus_u;
    cfg.atom_spin_guess = atom_spin_guess;

    // XL-BOMD guards (no-fallback): auxiliary-density guess on the
    // diagonalization path only. RKS: one total-P propagator. UKS: dual
    // α/β propagators (same Table-I coefficients).
    if (xl_bomd) {
        if (cfg.use_ot)
            throw std::invalid_argument(
                "GPWEvaluator: xl_bomd requires the diagonalization SCF path "
                "(use_ot parameterizes orbitals, not the density)");
        if (cpmd)
            throw std::invalid_argument(
                "GPWEvaluator: xl_bomd and cpmd are mutually exclusive "
                "(XL is SCF-converged auxiliary density; CP is fictitious "
                "orbital dynamics)");
        if (second_gen_cp)
            throw std::invalid_argument(
                "GPWEvaluator: xl_bomd and second_gen_cp are mutually "
                "exclusive (XL = Niklasson aux density; SGCP = ASPC + "
                "partial SCF corrector)");
        xl_.emplace(xl_bomd_k); // validates xl_bomd_k in [3, 9]
        if (n_alpha + n_beta > 0)
            xl_beta_.emplace(xl_bomd_k);
        else
            xl_beta_.reset();
    } else {
        xl_.reset();
        xl_beta_.reset();
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
        if (second_gen_cp)
            throw std::invalid_argument(
                "GPWEvaluator: cpmd and second_gen_cp are mutually exclusive "
                "(CPMD = fictitious mass; SGCP = ASPC + partial SCF)");
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

    // Second-generation CP (Kühne PRL 2007): ASPC + partial SCF corrector.
    if (second_gen_cp) {
        if (n_alpha + n_beta > 0)
            throw std::invalid_argument(
                "GPWEvaluator: second_gen_cp supports closed-shell RKS only "
                "(n_alpha/n_beta are set)");
        if (cfg.use_ot)
            throw std::invalid_argument(
                "GPWEvaluator: second_gen_cp requires the diagonalization "
                "SCF path (OT corrector is a follow-on)");
        if (!(second_gen_cp_k >= 0))
            throw std::invalid_argument(
                "GPWEvaluator: second_gen_cp_k must be >= 0 (ASPC order K)");
        if (!(second_gen_cp_max_iter >= 1))
            throw std::invalid_argument(
                "GPWEvaluator: second_gen_cp_max_iter must be >= 1 "
                "(corrector SCF budget)");
        if (!compute_forces)
            throw std::invalid_argument(
                "GPWEvaluator: second_gen_cp requires compute_forces=true "
                "(forces at the post-corrector density)");
        sgcp_aspc_.emplace(second_gen_cp_k);
    } else {
        sgcp_aspc_.reset();
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

    // ── Second-generation Car–Parrinello (ASPC + partial SCF) ───────────
    // Cold steps (insufficient history): full max_iter SCF, record P.
    // Warm steps: ASPC predict → SCF with corrector budget → forces at that
    // density → Kolafa correct + record (Kühne PRL 98, 066401 (2007)).
    if (sgcp_aspc_) {
        xt::xtensor<double, 2> P_pred;
        const bool have_pred = sgcp_aspc_->predict(P_pred);
        const ATV::cpu::linalg::RealMatrix* guess_c =
            (last_.C.size() > 0) ? &last_.C : nullptr;
        const ATV::cpu::linalg::RealMatrix* guess_p = nullptr;
        const int saved_max = solver_->max_iter;
        if (have_pred) {
            // Layout convert: ASPC uses default-layout xtensor; GPW is
            // column-major RealMatrix.
            const auto n0 = P_pred.shape(0);
            const auto n1 = P_pred.shape(1);
            sgcp_guess_p_ = xt::zeros<double>({n0, n1});
            for (std::size_t i = 0; i < n0; ++i)
                for (std::size_t j = 0; j < n1; ++j)
                    sgcp_guess_p_(i, j) = P_pred(i, j);
            guess_p = &sgcp_guess_p_;
            solver_->max_iter = second_gen_cp_max_iter;
        }
        try {
            last_ = solver_->run(atoms, guess_c, guess_p);
        } catch (...) {
            solver_->max_iter = saved_max;
            throw;
        }
        solver_->max_iter = saved_max;

        if (last_.P.size() > 0) {
            // Ensure H_ks / S for force-at-density when SCF stored orbitals
            // without the fixed-density payload.
            if (last_.H_ks.size() == 0 || last_.S.size() == 0) {
                const int saved = solver_->max_iter;
                solver_->max_iter = 0;
                try {
                    auto tmp = solver_->run(atoms, /*guess_C=*/nullptr,
                                            &last_.P);
                    last_.H_ks = std::move(tmp.H_ks);
                    last_.S = std::move(tmp.S);
                    // Keep energy / P / C / iters from the corrector SCF.
                } catch (...) {
                    solver_->max_iter = saved;
                    throw;
                }
                solver_->max_iter = saved;
            }
            xt::xtensor<double, 2> P_scf = xt::zeros<double>(
                {last_.P.shape(0), last_.P.shape(1)});
            for (std::size_t i = 0; i < last_.P.shape(0); ++i)
                for (std::size_t j = 0; j < last_.P.shape(1); ++j)
                    P_scf(i, j) = last_.P(i, j);
            if (have_pred) {
                xt::xtensor<double, 2> P_corr;
                sgcp_aspc_->correct(P_scf, P_pred, P_corr);
                sgcp_aspc_->record(P_corr);
            } else {
                sgcp_aspc_->record(P_scf);
            }
        } else {
            sgcp_aspc_->reset();
        }

        frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

        if (compute_forces) {
            if (last_.C.size() == 0 || last_.P.size() == 0
                || last_.H_ks.size() == 0)
                throw std::runtime_error(
                    "GPWEvaluator: second_gen_cp force path needs P, C, H_ks "
                    "from the corrector SCF");
            const int nb = static_cast<int>(last_.C.shape(0));
            const int n_occ = solver_->n_occ;
            linalg::RealMatrix C_occ = xt::zeros<double>(
                {static_cast<std::size_t>(nb),
                 static_cast<std::size_t>(n_occ)});
            for (int mu = 0; mu < nb; ++mu)
                for (int i = 0; i < n_occ; ++i)
                    C_occ(mu, i) = last_.C(mu, i);
            const auto F = gpw_forces_at_density(atoms, *solver_, last_.P,
                                                 C_occ, last_.H_ks);
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
        return;
    }

    // BOMD warm-start: reuse previous MO coefficients when shape matches.
    // GpwDft::run accepts optional guess_C; empty last_.C → core guess.
    const ATV::cpu::linalg::RealMatrix* guess = nullptr;
    if (last_.C.size() > 0) {
        guess = &last_.C;
    }
    // XL-BOMD: propagated auxiliary density is the SCF initial-density
    // guess (Niklasson dissipative recursion, seeded on the first step).
    // RKS: total P on guess_P. UKS: P_α / P_β on guess_P / guess_P_beta.
    const ATV::cpu::linalg::RealMatrix* guess_p =
        (xl_ && xl_->ready()) ? &xl_->aux() : nullptr;
    const ATV::cpu::linalg::RealMatrix* guess_p_beta =
        (xl_beta_ && xl_beta_->ready()) ? &xl_beta_->aux() : nullptr;
    last_ = solver_->run(atoms, guess, guess_p, guess_p_beta);
    if (xl_) {
        if (last_.converged && last_.P.size() > 0) {
            xl_->advance(last_.P);
            if (xl_beta_) {
                if (last_.P_beta.size() > 0)
                    xl_beta_->advance(last_.P_beta);
                else
                    xl_beta_->reset();
            }
        } else {
            // Discontinuity: never propagate an unconverged D.
            xl_->reset();
            if (xl_beta_)
                xl_beta_->reset();
        }
    }
    frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

    if (compute_forces) {
        // Prefer forces at the density just converged (no second SCF). UKS has
        // no guess_C warm-start, so a re-SCF would double wall; RKS cold
        // re-SCF is kept only when the energy result lacks orbitals (should
        // not happen after a normal converged run).
        std::vector<std::array<double, 3>> F;
        if (last_.converged && last_.P.size() > 0 && last_.C.size() > 0
            && !last_.eps.empty()) {
            F = gpw_forces_from_scf(atoms, *solver_, last_);
        } else {
            // XL mode warm-starts the force-pass SCF from the MO coefficients
            // just converged at this exact geometry. Default mode keeps the
            // cold solve unchanged.
            F = gpw_forces_analytic(
                atoms, *solver_,
                (xl_ && last_.C.size() > 0) ? &last_.C : nullptr);
        }
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
