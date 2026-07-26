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

    // Car–Parrinello guards: RKS Γ or collinear Γ-UKS diagonalization path.
    if (cpmd) {
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
        // UKS CPMD: alpha is the primary channel (mirrors XL dual setup).
        if (n_alpha + n_beta > 0 && n_alpha < 1)
            throw std::invalid_argument(
                "GPWEvaluator: cpmd UKS requires n_alpha >= 1");
        cpmd_.emplace(cpmd_mu);
        if (n_alpha + n_beta > 0 && n_beta > 0)
            cpmd_beta_.emplace(cpmd_mu);
        else
            cpmd_beta_.reset();
        last_te_ev_ = 0.0;
    } else {
        cpmd_.reset();
        cpmd_beta_.reset();
        last_te_ev_ = 0.0;
    }

    // Second-generation CP (Kühne PRL 2007): ASPC + partial SCF corrector.
    // RKS: one total-P ASPC. UKS: dual α/β ASPC (same order K).
    if (second_gen_cp) {
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
        if (n_alpha + n_beta > 0)
            sgcp_aspc_beta_.emplace(second_gen_cp_k);
        else
            sgcp_aspc_beta_.reset();
    } else {
        sgcp_aspc_.reset();
        sgcp_aspc_beta_.reset();
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
    // RKS: one propagator (occ=2). UKS: dual α/β (occ=1 each).
    const bool cpmd_uks = cpmd_ && (n_alpha + n_beta > 0);
    const bool cpmd_ready =
        cpmd_
        && (cpmd_uks ? (cpmd_->ready()
                        && (n_beta < 1 || (cpmd_beta_ && cpmd_beta_->ready())))
                     : cpmd_->ready());
    if (cpmd_ready) {
        auto fixed_density_eval_rks = [&](const linalg::RealMatrix& C_occ) {
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
        auto fixed_density_eval_uks = [&](const linalg::RealMatrix& Ca,
                                          const linalg::RealMatrix& Cb) {
            const int nb = static_cast<int>(
                (Ca.size() > 0 ? Ca.shape(0) : Cb.shape(0)));
            const int na = Ca.size() > 0 ? static_cast<int>(Ca.shape(1)) : 0;
            const int nb_spin =
                Cb.size() > 0 ? static_cast<int>(Cb.shape(1)) : 0;
            const linalg::RealMatrix Pa =
                na > 0 ? linalg::build_density_matrix(Ca, na, nb, 1.0)
                       : linalg::RealMatrix(
                             xt::zeros<double>({static_cast<std::size_t>(nb),
                                                static_cast<std::size_t>(nb)}));
            const linalg::RealMatrix Pb =
                nb_spin > 0
                    ? linalg::build_density_matrix(Cb, nb_spin, nb, 1.0)
                    : linalg::RealMatrix(xt::zeros<double>(
                          {static_cast<std::size_t>(nb),
                           static_cast<std::size_t>(nb)}));
            const int saved_max = solver_->max_iter;
            solver_->max_iter = 0;
            try {
                last_ = solver_->run(atoms, /*guess_C=*/nullptr, &Pa, &Pb);
            } catch (...) {
                solver_->max_iter = saved_max;
                throw;
            }
            solver_->max_iter = saved_max;
        };

        if (cpmd_uks) {
            // Empty beta C when n_beta==0 (half-metal edge).
            const linalg::RealMatrix empty_c;
            const linalg::RealMatrix& Ca = cpmd_->C();
            const linalg::RealMatrix& Cb =
                (cpmd_beta_ && cpmd_beta_->ready()) ? cpmd_beta_->C()
                                                    : empty_c;
            fixed_density_eval_uks(Ca, Cb);
            cpmd_->align_to_overlap(last_.S);
            if (cpmd_beta_ && cpmd_beta_->ready())
                cpmd_beta_->align_to_overlap(last_.S);
            const linalg::RealMatrix& Ca2 = cpmd_->C();
            const linalg::RealMatrix& Cb2 =
                (cpmd_beta_ && cpmd_beta_->ready()) ? cpmd_beta_->C()
                                                    : empty_c;
            fixed_density_eval_uks(Ca2, Cb2);

            frame.meta.get(kEnergy) =
                last_.e_total * ATV::units::kHartreeInEv;

            if (compute_forces) {
                if (last_.H_ks.size() == 0 || last_.H_ks_beta.size() == 0)
                    throw std::runtime_error(
                        "GPWEvaluator: UKS CPMD force path needs H_ks and "
                        "H_ks_beta from fixed-density eval");
                const auto F = gpw_forces_at_density(
                    atoms, *solver_, last_.P, Ca2, last_.H_ks, last_.P_beta,
                    Cb2, last_.H_ks_beta);
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

            // Dual electronic steps (occ=1 per spin).
            double te_ha = 0.0;
            cpmd_->advance(last_.H_ks, last_.S, cpmd_dt, /*occ=*/1.0);
            te_ha += cpmd_->kinetic_ha();
            if (cpmd_beta_ && cpmd_beta_->ready()) {
                cpmd_beta_->advance(last_.H_ks_beta, last_.S, cpmd_dt,
                                    /*occ=*/1.0);
                te_ha += cpmd_beta_->kinetic_ha();
            }
            last_te_ev_ = te_ha * ATV::units::kHartreeInEv;
            return;
        }

        // Closed-shell RKS CPMD.
        fixed_density_eval_rks(cpmd_->C());
        cpmd_->align_to_overlap(last_.S);
        fixed_density_eval_rks(cpmd_->C());

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

        cpmd_->advance(last_.H_ks, last_.S, cpmd_dt, /*occ=*/2.0);
        last_te_ev_ = cpmd_->kinetic_ha() * ATV::units::kHartreeInEv;
        return;
    }

    // ── Second-generation Car–Parrinello (ASPC + partial SCF) ───────────
    // Cold steps (insufficient history): full max_iter SCF, record P.
    // Warm steps: ASPC predict → SCF with corrector budget → forces at that
    // density → Kolafa correct + record (Kühne PRL 98, 066401 (2007)).
    // RKS: total P. UKS: dual ASPC on P_α / P_β (same order K).
    if (sgcp_aspc_) {
        const bool uks_sgcp = static_cast<bool>(sgcp_aspc_beta_);
        xt::xtensor<double, 2> P_pred;
        xt::xtensor<double, 2> P_pred_beta;
        const bool have_pred_a = sgcp_aspc_->predict(P_pred);
        bool have_pred_b = false;
        if (uks_sgcp)
            have_pred_b = sgcp_aspc_beta_->predict(P_pred_beta);
        // Both channels must be warm together (shared MD step count).
        const bool have_pred =
            uks_sgcp ? (have_pred_a && have_pred_b) : have_pred_a;
        const ATV::cpu::linalg::RealMatrix* guess_c =
            (last_.C.size() > 0) ? &last_.C : nullptr;
        const ATV::cpu::linalg::RealMatrix* guess_p = nullptr;
        const ATV::cpu::linalg::RealMatrix* guess_p_beta = nullptr;
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
            if (uks_sgcp) {
                const auto m0 = P_pred_beta.shape(0);
                const auto m1 = P_pred_beta.shape(1);
                sgcp_guess_p_beta_ = xt::zeros<double>({m0, m1});
                for (std::size_t i = 0; i < m0; ++i)
                    for (std::size_t j = 0; j < m1; ++j)
                        sgcp_guess_p_beta_(i, j) = P_pred_beta(i, j);
                guess_p_beta = &sgcp_guess_p_beta_;
            }
            solver_->max_iter = second_gen_cp_max_iter;
        }
        try {
            last_ = solver_->run(atoms, guess_c, guess_p, guess_p_beta);
        } catch (...) {
            solver_->max_iter = saved_max;
            throw;
        }
        solver_->max_iter = saved_max;

        auto record_channel = [](scf::ASPCPredictor& aspc,
                                 const linalg::RealMatrix& Pmat,
                                 bool do_correct,
                                 const xt::xtensor<double, 2>& P_pred_in) {
            xt::xtensor<double, 2> P_scf =
                xt::zeros<double>({Pmat.shape(0), Pmat.shape(1)});
            for (std::size_t i = 0; i < Pmat.shape(0); ++i)
                for (std::size_t j = 0; j < Pmat.shape(1); ++j)
                    P_scf(i, j) = Pmat(i, j);
            if (do_correct) {
                xt::xtensor<double, 2> P_corr;
                aspc.correct(P_scf, P_pred_in, P_corr);
                aspc.record(P_corr);
            } else {
                aspc.record(P_scf);
            }
        };

        if (last_.P.size() > 0) {
            if (!uks_sgcp) {
                // RKS force-at-density needs H_ks / S when the SCF stored
                // orbitals without the fixed-density payload.
                if (last_.H_ks.size() == 0 || last_.S.size() == 0) {
                    const int saved = solver_->max_iter;
                    solver_->max_iter = 0;
                    try {
                        auto tmp = solver_->run(atoms, /*guess_C=*/nullptr,
                                                &last_.P);
                        last_.H_ks = std::move(tmp.H_ks);
                        last_.S = std::move(tmp.S);
                    } catch (...) {
                        solver_->max_iter = saved;
                        throw;
                    }
                    solver_->max_iter = saved;
                }
            }
            record_channel(*sgcp_aspc_, last_.P, have_pred, P_pred);
            if (uks_sgcp) {
                if (last_.P_beta.size() > 0)
                    record_channel(*sgcp_aspc_beta_, last_.P_beta, have_pred,
                                   P_pred_beta);
                else
                    sgcp_aspc_beta_->reset();
            }
        } else {
            sgcp_aspc_->reset();
            if (sgcp_aspc_beta_)
                sgcp_aspc_beta_->reset();
        }

        frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

        if (compute_forces) {
            std::vector<std::array<double, 3>> F;
            if (uks_sgcp) {
                // UKS: forces at the post-corrector density (partial SCF may
                // not hit tol). gpw_forces_from_scf requires the converged
                // flag; SGCP semantics are forces at the corrector state.
                if (last_.P.size() == 0 || last_.C.size() == 0
                    || last_.eps.empty() || last_.P_beta.size() == 0
                    || last_.C_beta.size() == 0 || last_.eps_beta.empty())
                    throw std::runtime_error(
                        "GPWEvaluator: second_gen_cp UKS force path needs "
                        "P/C/eps for both spins from the corrector SCF");
                GpwResult force_scf = last_;
                force_scf.converged = true;
                F = gpw_forces_from_scf(atoms, *solver_, force_scf);
            } else {
                if (last_.C.size() == 0 || last_.P.size() == 0
                    || last_.H_ks.size() == 0)
                    throw std::runtime_error(
                        "GPWEvaluator: second_gen_cp force path needs P, C, "
                        "H_ks from the corrector SCF");
                const int nb = static_cast<int>(last_.C.shape(0));
                const int n_occ = solver_->n_occ;
                linalg::RealMatrix C_occ = xt::zeros<double>(
                    {static_cast<std::size_t>(nb),
                     static_cast<std::size_t>(n_occ)});
                for (int mu = 0; mu < nb; ++mu)
                    for (int i = 0; i < n_occ; ++i)
                        C_occ(mu, i) = last_.C(mu, i);
                F = gpw_forces_at_density(atoms, *solver_, last_.P, C_occ,
                                          last_.H_ks);
            }
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
        if (!last_.converged || last_.C.size() == 0)
            throw std::runtime_error(
                "GPWEvaluator: cpmd seed requires a converged SCF with "
                "occupied MO coefficients");
        const int nb = static_cast<int>(last_.C.shape(0));
        if (n_alpha + n_beta > 0) {
            // Collinear UKS: seed α (and β when present) occupied blocks.
            if (n_alpha < 1)
                throw std::runtime_error(
                    "GPWEvaluator: UKS CPMD seed requires n_alpha >= 1 "
                    "(alpha propagator is primary)");
            if (static_cast<int>(last_.C.shape(1)) < n_alpha)
                throw std::runtime_error(
                    "GPWEvaluator: UKS CPMD seed: alpha MO count < n_alpha");
            linalg::RealMatrix Ca = xt::zeros<double>(
                {static_cast<std::size_t>(nb),
                 static_cast<std::size_t>(n_alpha)});
            for (int mu = 0; mu < nb; ++mu)
                for (int i = 0; i < n_alpha; ++i)
                    Ca(mu, i) = last_.C(mu, i);
            cpmd_->seed(Ca);
            if (n_beta > 0) {
                if (!cpmd_beta_)
                    throw std::runtime_error(
                        "GPWEvaluator: UKS CPMD missing beta propagator");
                if (last_.C_beta.size() == 0
                    || static_cast<int>(last_.C_beta.shape(1)) < n_beta)
                    throw std::runtime_error(
                        "GPWEvaluator: UKS CPMD seed requires beta MO "
                        "coefficients");
                linalg::RealMatrix Cb = xt::zeros<double>(
                    {static_cast<std::size_t>(nb),
                     static_cast<std::size_t>(n_beta)});
                for (int mu = 0; mu < nb; ++mu)
                    for (int i = 0; i < n_beta; ++i)
                        Cb(mu, i) = last_.C_beta(mu, i);
                cpmd_beta_->seed(Cb);
            }
        } else {
            if (solver_->n_occ < 1)
                throw std::runtime_error(
                    "GPWEvaluator: RKS CPMD seed requires n_occ >= 1");
            const int n_occ = solver_->n_occ;
            linalg::RealMatrix C_occ = xt::zeros<double>(
                {static_cast<std::size_t>(nb),
                 static_cast<std::size_t>(n_occ)});
            for (int mu = 0; mu < nb; ++mu)
                for (int i = 0; i < n_occ; ++i)
                    C_occ(mu, i) = last_.C(mu, i);
            cpmd_->seed(C_occ);
        }
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
