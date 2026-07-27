/// @file src/cpu/gpw/gpw_evaluator.cpp
/// @brief Implementation of the GPW energy evaluator component (see header).

#include "atomiverse/cpu/gpw/gpw_evaluator.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/views/xview.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/gpw/internal/gpw_forces.h"
#include "atomiverse/cpu/gpw/internal/gpw_stress.h"
#include "atomiverse/cpu/linalg/density.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace gpw {

namespace {
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
        const double pos[3] = {sp_x[i] * units::kAngstromInBohr,
                               sp_y[i] * units::kAngstromInBohr,
                               sp_z[i] * units::kAngstromInBohr};
        atoms.push_back(GpwAtom(pos, k.basis, k.pot, sp_type[i]));
    }
    return atoms;
}

/// Publish the two CPMD electronic energies [eV] on the frame, so Observers,
/// Thermo and `cpu::md::CpmdConserved` see the adiabaticity signal that used
/// to be reachable only through `last_electronic_ke_ev()`. Called from every
/// CPMD branch (Γ RKS, Γ UKS, multi-k RKS) plus the cold seed, which
/// publishes zeros — the electronic degrees of freedom have not moved yet.
/// The keys are open-schema (not pre-inserted by `Frame(int)`), so the first
/// launch inserts and later launches overwrite (the `kDipole` pattern).
/// @param frame Frame to write.
/// @param te_ev Fictitious electronic kinetic energy \(T_e\) [eV].
/// @param etherm_ev Electronic Nosé–Hoover thermostat energy [eV].
void publish_electronic_energies(cpu::Frame& frame, double te_ev,
                                 double etherm_ev) {
    if (!frame.meta.contains(kEnergyElectronKinetic.key))
        frame.meta.insert(kEnergyElectronKinetic, te_ev);
    else
        frame.meta.set(kEnergyElectronKinetic, te_ev);
    if (!frame.meta.contains(kEnergyElectronThermostat.key))
        frame.meta.insert(kEnergyElectronThermostat, etherm_ev);
    else
        frame.meta.set(kEnergyElectronThermostat, etherm_ev);
}
} // namespace

void GPWEvaluator::set_timestep(double dt_fs) {
    if (!(dt_fs > 0.0) || !std::isfinite(dt_fs))
        throw std::invalid_argument(
            "GPWEvaluator::set_timestep: nuclear timestep must be finite and "
            "> 0 (fs), got "
            + std::to_string(dt_fs));
    cpmd_dt_fs_ = dt_fs;
}

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
    cfg.cell[0] = ctx.frame.box.length(0) * units::kAngstromInBohr;
    cfg.cell[1] = ctx.frame.box.length(1) * units::kAngstromInBohr;
    cfg.cell[2] = ctx.frame.box.length(2) * units::kAngstromInBohr;
    cfg.cutoff_ry = cutoff_ry;
    cfg.ngrids = ngrids;
    cfg.rel_cutoff = rel_cutoff;
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
    cfg.d3_atm = d3_atm;
    cfg.d4 = d4;
    cfg.d4_atm = d4_atm;
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
    for (int d = 0; d < 3; ++d)
        cfg.k_mesh[d] = k_mesh[d];
    cfg.k_time_reversal = k_time_reversal;

    const bool multi_k =
        (k_mesh[0] != 1 || k_mesh[1] != 1 || k_mesh[2] != 1);

    // XL-BOMD guards (no-fallback): auxiliary-density guess on the
    // diagonalization path only.
    //   Γ RKS: one total-P real propagator
    //   Γ UKS: dual α/β real propagators
    //   multi-k RKS: one complex propagator per k (seeded on first SCF)
    //   multi-k UKS: dual complex propagators per k (P_k^α, P_k^β)
    if (xl_bomd_linearized && !xl_bomd)
        throw std::invalid_argument(
            "GPWEvaluator: xl_bomd_linearized requires xl_bomd=true "
            "(zero-SCF shadow forces ride on the XL auxiliary density)");
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
        if (xl_bomd_linearized) {
            // Zero-SCF linearized forces: Γ RKS diagonalization foothold.
            if (!compute_forces)
                throw std::invalid_argument(
                    "GPWEvaluator: xl_bomd_linearized requires "
                    "compute_forces=true (shadow forces are the product)");
            if (n_alpha + n_beta > 0)
                throw std::invalid_argument(
                    "GPWEvaluator: xl_bomd_linearized is Γ RKS only "
                    "(UKS linearized XL is a follow-on)");
            if (multi_k)
                throw std::invalid_argument(
                    "GPWEvaluator: xl_bomd_linearized is Γ-point only "
                    "(multi-k linearized XL is a follow-on)");
        }
        // Validates xl_bomd_k ∈ [3, 9] via the real constructor even when
        // multi-k uses XlBomdCplx (same Table I).
        if (multi_k) {
            xl_.reset();
            xl_beta_.reset();
            xl_k_.clear(); // sized after the first multi-k SCF (n_kpts known)
            xl_k_beta_.clear();
            (void)electronic::XlBomd(xl_bomd_k); // throw early on bad K
        } else {
            xl_.emplace(xl_bomd_k);
            if (n_alpha + n_beta > 0)
                xl_beta_.emplace(xl_bomd_k);
            else
                xl_beta_.reset();
            xl_k_.clear();
            xl_k_beta_.clear();
        }
    } else {
        xl_.reset();
        xl_beta_.reset();
        xl_k_.clear();
        xl_k_beta_.clear();
    }

    // Car–Parrinello guards: RKS Γ, collinear Γ-UKS, or multi-k RKS.
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
        // The electronic step rides the nuclear clock, so the dt is PUSHED by
        // the MD integrator (`cpu::TimestepProvider` → `set_timestep`), never
        // re-declared here. No push ⇒ either there is no MD loop around this
        // evaluator, or its integrator publishes nothing: both are broken CP.
        if (!(cpmd_dt_fs_ > 0.0))
            throw std::runtime_error(
                "GPWEvaluator: cpmd requires an MD integrator publishing the "
                "timestep (TimestepProvider); CPMD outside MD is meaningless");
        if (!compute_forces)
            throw std::invalid_argument(
                "GPWEvaluator: cpmd requires compute_forces=true (nuclear "
                "forces from the fixed-density electronic state)");
        // Every CPMD branch of `launch` returns before the stress block, so a
        // stress request would be silently dropped rather than refused.
        if (compute_stress)
            throw std::runtime_error(
                "GPWEvaluator: CPMD stress path is not productized; "
                "variable-cell CP is deferred");
        // Multi-k UKS CPMD is a residual (dual-spin Bloch C_k^σ).
        if (multi_k && (n_alpha + n_beta > 0))
            throw std::invalid_argument(
                "GPWEvaluator: multi-k UKS + cpmd is not productized "
                "(use Γ-UKS CPMD or multi-k RKS CPMD)");
        // UKS CPMD: alpha is the primary channel (mirrors XL dual setup).
        if (n_alpha + n_beta > 0 && n_alpha < 1)
            throw std::invalid_argument(
                "GPWEvaluator: cpmd UKS requires n_alpha >= 1");
        if (cpmd_thermostat) {
            if (multi_k)
                throw std::invalid_argument(
                    "GPWEvaluator: multi-k CPMD electronic thermostat is not "
                    "productized (NVE electrons only on the k path)");
            if (!(cpmd_ekin_target_ha > 0.0)
                || !std::isfinite(cpmd_ekin_target_ha))
                throw std::invalid_argument(
                    "GPWEvaluator: cpmd_thermostat requires "
                    "cpmd_ekin_target_ha finite and > 0 (Ha)");
            if (!(cpmd_tau_e > 0.0) || !std::isfinite(cpmd_tau_e))
                throw std::invalid_argument(
                    "GPWEvaluator: cpmd_thermostat requires cpmd_tau_e "
                    "finite and > 0 (fs)");
        }
        if (multi_k) {
            // Complex propagators sized after the first multi-k SCF (n_kpts).
            cpmd_.reset();
            cpmd_beta_.reset();
            cpmd_k_.clear();
            (void)electronic::CpmdCplx(cpmd_mu); // early mu validation
        } else {
            cpmd_.emplace(cpmd_mu);
            cpmd_->thermostat = cpmd_thermostat;
            cpmd_->ekin_target_ha = cpmd_ekin_target_ha;
            cpmd_->tau_e_fs = cpmd_tau_e;
            if (n_alpha + n_beta > 0 && n_beta > 0) {
                cpmd_beta_.emplace(cpmd_mu);
                cpmd_beta_->thermostat = cpmd_thermostat;
                cpmd_beta_->ekin_target_ha = cpmd_ekin_target_ha;
                cpmd_beta_->tau_e_fs = cpmd_tau_e;
            } else {
                cpmd_beta_.reset();
            }
            cpmd_k_.clear();
        }
        last_te_ev_ = 0.0;
        last_etherm_ev_ = 0.0;
    } else {
        if (cpmd_thermostat)
            throw std::invalid_argument(
                "GPWEvaluator: cpmd_thermostat requires cpmd=true");
        cpmd_.reset();
        cpmd_beta_.reset();
        cpmd_k_.clear();
        last_te_ev_ = 0.0;
        last_etherm_ev_ = 0.0;
    }

    // Second-generation CP (Kühne PRL 2007): ASPC + partial SCF corrector.
    // Γ RKS: one total-P real ASPC (DIIS or OT corrector).
    // Γ UKS: dual α/β real ASPC (DIIS or OT).
    // multi-k RKS: one complex ASPC per k (DIIS only — OT + P_k throws).
    // multi-k UKS: dual complex ASPC per k.
    if (second_gen_cp) {
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
        if (multi_k) {
            if (cfg.use_ot)
                throw std::invalid_argument(
                    "GPWEvaluator: multi-k second_gen_cp requires the "
                    "diagonalization SCF path (use_ot + P_k density guess is "
                    "not productized)");
            // Validate K early; k-histories sized after the first multi-k SCF
            // (n_kpts known).
            (void)electronic::ASPCPredictorCplx(second_gen_cp_k);
            sgcp_aspc_.reset();
            sgcp_aspc_beta_.reset();
            sgcp_aspc_k_.clear();
            sgcp_aspc_k_beta_.clear();
        } else {
            sgcp_aspc_.emplace(second_gen_cp_k);
            if (n_alpha + n_beta > 0)
                sgcp_aspc_beta_.emplace(second_gen_cp_k);
            else
                sgcp_aspc_beta_.reset();
            sgcp_aspc_k_.clear();
            sgcp_aspc_k_beta_.clear();
        }
    } else {
        sgcp_aspc_.reset();
        sgcp_aspc_beta_.reset();
        sgcp_aspc_k_.clear();
        sgcp_aspc_k_beta_.clear();
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
    // launches is well-defined — cell-relax drivers rescale box freely.
    for (int d = 0; d < 3; ++d) {
        const double now = frame.box.length(d) * units::kAngstromInBohr;
        if (!(now > 0.0))
            throw std::runtime_error(
                "GPWEvaluator::launch: non-positive box edge");
        solver_->cell[d] = now;
    }

    const auto atoms = assemble_atoms(frame, *this);

    // ── Car–Parrinello path ──────────────────────────────────────────────
    // Seed: one full SCF + BO analytic forces, then hand occupied C to Cpmd.
    // Later steps: fixed-density KS eval (max_iter=0) + forces at that
    // density + electronic Verlet. Nuclear PE is E_KS only; T_e is published
    // on frame.meta[kEnergyElectronKinetic] (and last_electronic_ke_ev())
    // for cpu::md::CpmdConserved to assemble the CP conserved energy.
    // RKS Γ: one real propagator (occ=2). UKS Γ: dual α/β (occ=1 each).
    // Multi-k RKS: one CpmdCplx per sampled k (occ=2).
    const bool cpmd_multi_k =
        !cpmd_k_.empty() && cpmd_k_[0].ready()
        && [&]() {
               for (const auto& ch : cpmd_k_)
                   if (!ch.ready())
                       return false;
               return true;
           }();
    const bool cpmd_uks = cpmd_ && (n_alpha + n_beta > 0);
    const bool cpmd_ready =
        cpmd_multi_k
        || (cpmd_
            && (cpmd_uks
                    ? (cpmd_->ready()
                       && (n_beta < 1
                           || (cpmd_beta_ && cpmd_beta_->ready())))
                    : cpmd_->ready()));
    if (cpmd_ready && cpmd_multi_k) {
        // Multi-k RKS CPMD: fixed-density {P_k} → forces → electronic VV.
        using Cplx = std::complex<double>;
        const int nk = static_cast<int>(cpmd_k_.size());
        auto density_from_cocc = [](const linalg::CplxMatrix& C,
                                    int n_occ) -> linalg::CplxMatrix {
            const int nb = static_cast<int>(C.shape(0));
            linalg::CplxMatrix P = xt::zeros<Cplx>(
                {static_cast<std::size_t>(nb), static_cast<std::size_t>(nb)});
            for (int o = 0; o < n_occ; ++o)
                for (int mu = 0; mu < nb; ++mu) {
                    const Cplx cmu = C(mu, o);
                    for (int nu = 0; nu < nb; ++nu)
                        P(mu, nu) += 2.0 * cmu * std::conj(C(nu, o));
                }
            return P;
        };
        auto fixed_density_eval_k = [&]() {
            const int n_occ = solver_->n_occ;
            cpmd_guess_p_k_.resize(static_cast<std::size_t>(nk));
            for (int ik = 0; ik < nk; ++ik)
                cpmd_guess_p_k_[static_cast<std::size_t>(ik)] =
                    density_from_cocc(cpmd_k_[static_cast<std::size_t>(ik)].C(),
                                      n_occ);
            const int saved_max = solver_->max_iter;
            solver_->max_iter = 0;
            try {
                last_ = solver_->run(atoms, /*guess_C=*/nullptr,
                                     /*guess_P=*/nullptr,
                                     /*guess_P_beta=*/nullptr, &cpmd_guess_p_k_);
            } catch (...) {
                solver_->max_iter = saved_max;
                throw;
            }
            solver_->max_iter = saved_max;
        };

        fixed_density_eval_k();
        if (static_cast<int>(last_.S_k.size()) != nk
            || static_cast<int>(last_.H_ks_k.size()) != nk)
            throw std::runtime_error(
                "GPWEvaluator: multi-k CPMD fixed-density missing H_ks_k/S_k");
        for (int ik = 0; ik < nk; ++ik)
            cpmd_k_[static_cast<std::size_t>(ik)].align_to_overlap(
                last_.S_k[static_cast<std::size_t>(ik)]);
        fixed_density_eval_k();

        frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

        if (compute_forces) {
            // Populate C_occ_k + eps_k = diag(C† H C) for k-summed forces.
            const int n_occ = solver_->n_occ;
            last_.C_occ_k.resize(static_cast<std::size_t>(nk));
            last_.eps_k.resize(static_cast<std::size_t>(nk));
            for (int ik = 0; ik < nk; ++ik) {
                const auto& C =
                    cpmd_k_[static_cast<std::size_t>(ik)].C();
                const auto& H = last_.H_ks_k[static_cast<std::size_t>(ik)];
                const int nb = static_cast<int>(C.shape(0));
                std::vector<Cplx> flat(
                    static_cast<std::size_t>(nb)
                    * static_cast<std::size_t>(n_occ));
                for (int mu = 0; mu < nb; ++mu)
                    for (int o = 0; o < n_occ; ++o)
                        flat[static_cast<std::size_t>(mu) * n_occ
                             + static_cast<std::size_t>(o)] = C(mu, o);
                last_.C_occ_k[static_cast<std::size_t>(ik)] = std::move(flat);
                // ε_i ≈ Re (C† H C)_ii  (BO surface → exact eigenvalues).
                std::vector<double> eps(static_cast<std::size_t>(n_occ), 0.0);
                for (int o = 0; o < n_occ; ++o) {
                    Cplx s(0, 0);
                    for (int mu = 0; mu < nb; ++mu) {
                        Cplx hco(0, 0);
                        for (int nu = 0; nu < nb; ++nu)
                            hco += H(mu, nu) * C(nu, o);
                        s += std::conj(C(mu, o)) * hco;
                    }
                    eps[static_cast<std::size_t>(o)] = s.real();
                }
                last_.eps_k[static_cast<std::size_t>(ik)] = std::move(eps);
            }
            last_.converged = true;
            const auto F = gpw_forces_from_scf(atoms, *solver_, last_);
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

        // Electronic steps: T_e = Σ_k w_k T_e^{(k)}.
        double te_ha = 0.0;
        for (int ik = 0; ik < nk; ++ik) {
            cpmd_k_[static_cast<std::size_t>(ik)].advance(
                last_.H_ks_k[static_cast<std::size_t>(ik)],
                last_.S_k[static_cast<std::size_t>(ik)], cpmd_dt_fs_,
                /*occ=*/2.0);
            const double wk =
                (ik < static_cast<int>(last_.k_weights.size()))
                    ? last_.k_weights[static_cast<std::size_t>(ik)]
                    : (1.0 / static_cast<double>(nk));
            te_ha += wk * cpmd_k_[static_cast<std::size_t>(ik)].kinetic_ha();
        }
        last_te_ev_ = te_ha * ATV::units::kHartreeInEv;
        last_etherm_ev_ = 0.0;
        publish_electronic_energies(frame, last_te_ev_, last_etherm_ev_);
        return;
    }
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
            double eth_ha = 0.0;
            cpmd_->advance(last_.H_ks, last_.S, cpmd_dt_fs_, /*occ=*/1.0);
            te_ha += cpmd_->kinetic_ha();
            eth_ha += cpmd_->thermostat_energy_ha();
            if (cpmd_beta_ && cpmd_beta_->ready()) {
                cpmd_beta_->advance(last_.H_ks_beta, last_.S, cpmd_dt_fs_,
                                    /*occ=*/1.0);
                te_ha += cpmd_beta_->kinetic_ha();
                eth_ha += cpmd_beta_->thermostat_energy_ha();
            }
            last_te_ev_ = te_ha * ATV::units::kHartreeInEv;
            last_etherm_ev_ = eth_ha * ATV::units::kHartreeInEv;
            publish_electronic_energies(frame, last_te_ev_, last_etherm_ev_);
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

        cpmd_->advance(last_.H_ks, last_.S, cpmd_dt_fs_, /*occ=*/2.0);
        last_te_ev_ = cpmd_->kinetic_ha() * ATV::units::kHartreeInEv;
        last_etherm_ev_ =
            cpmd_->thermostat_energy_ha() * ATV::units::kHartreeInEv;
        publish_electronic_energies(frame, last_te_ev_, last_etherm_ev_);
        return;
    }

    // ── Second-generation Car–Parrinello (ASPC + partial SCF) ───────────
    // Cold steps (insufficient history): full max_iter SCF, record P.
    // Warm steps: ASPC predict → SCF with corrector budget → forces at that
    // density → Kolafa correct + record (Kühne PRL 98, 066401 (2007)).
    // Γ RKS: total real P. Γ UKS: dual real P_α / P_β.
    // multi-k RKS: complex P_k per k. multi-k UKS: dual P_k^α / P_k^β.
    if (second_gen_cp) {
        const bool multi_k_sgcp =
            (k_mesh[0] != 1 || k_mesh[1] != 1 || k_mesh[2] != 1);
        const bool uks_sgcp = (n_alpha + n_beta > 0);

        // ── multi-k complex ASPC on P_k ────────────────────────────────
        if (multi_k_sgcp) {
            const std::vector<linalg::CplxMatrix>* guess_p_k = nullptr;
            const std::vector<linalg::CplxMatrix>* guess_p_k_beta = nullptr;
            std::vector<linalg::CplxMatrix> P_pred_k;
            std::vector<linalg::CplxMatrix> P_pred_k_beta;
            bool have_pred = false;
            if (!sgcp_aspc_k_.empty()
                && static_cast<int>(sgcp_aspc_k_.size()) > 0
                && sgcp_aspc_k_[0].history_size() >= second_gen_cp_k + 1) {
                // All k share the MD step count; require every channel warm.
                have_pred = true;
                P_pred_k.resize(sgcp_aspc_k_.size());
                for (std::size_t ik = 0; ik < sgcp_aspc_k_.size(); ++ik) {
                    if (!sgcp_aspc_k_[ik].predict(P_pred_k[ik])) {
                        have_pred = false;
                        break;
                    }
                }
                if (have_pred && uks_sgcp) {
                    if (sgcp_aspc_k_beta_.size() != sgcp_aspc_k_.size())
                        have_pred = false;
                    else {
                        P_pred_k_beta.resize(sgcp_aspc_k_beta_.size());
                        for (std::size_t ik = 0; ik < sgcp_aspc_k_beta_.size();
                             ++ik) {
                            if (!sgcp_aspc_k_beta_[ik].predict(
                                    P_pred_k_beta[ik])) {
                                have_pred = false;
                                break;
                            }
                        }
                    }
                }
                if (have_pred) {
                    sgcp_guess_p_k_ = P_pred_k;
                    guess_p_k = &sgcp_guess_p_k_;
                    if (uks_sgcp) {
                        sgcp_guess_p_k_beta_ = P_pred_k_beta;
                        guess_p_k_beta = &sgcp_guess_p_k_beta_;
                    }
                }
            }

            const int saved_max = solver_->max_iter;
            if (have_pred)
                solver_->max_iter = second_gen_cp_max_iter;
            try {
                last_ = solver_->run(atoms, /*guess_C=*/nullptr,
                                     /*guess_P=*/nullptr,
                                     /*guess_P_beta=*/nullptr, guess_p_k,
                                     guess_p_k_beta);
            } catch (...) {
                solver_->max_iter = saved_max;
                throw;
            }
            solver_->max_iter = saved_max;

            auto reseat_k =
                [&](std::vector<electronic::ASPCPredictorCplx>& xs) {
                    if (static_cast<int>(xs.size()) != last_.n_kpts) {
                        xs.clear();
                        xs.reserve(static_cast<std::size_t>(last_.n_kpts));
                        for (int ik = 0; ik < last_.n_kpts; ++ik)
                            xs.emplace_back(second_gen_cp_k);
                    }
                };

            const bool have_pk =
                last_.n_kpts > 0
                && static_cast<int>(last_.P_k.size()) == last_.n_kpts;
            const bool have_pk_b =
                !uks_sgcp
                || (have_pk
                    && static_cast<int>(last_.P_k_beta.size())
                           == last_.n_kpts);
            if (have_pk && have_pk_b) {
                reseat_k(sgcp_aspc_k_);
                for (int ik = 0; ik < last_.n_kpts; ++ik) {
                    const auto& Pk = last_.P_k[static_cast<std::size_t>(ik)];
                    auto& aspc = sgcp_aspc_k_[static_cast<std::size_t>(ik)];
                    if (have_pred
                        && static_cast<int>(P_pred_k.size()) == last_.n_kpts) {
                        linalg::CplxMatrix P_corr;
                        aspc.correct(Pk, P_pred_k[static_cast<std::size_t>(ik)],
                                     P_corr);
                        aspc.record(P_corr);
                    } else {
                        aspc.record(Pk);
                    }
                }
                if (uks_sgcp) {
                    reseat_k(sgcp_aspc_k_beta_);
                    for (int ik = 0; ik < last_.n_kpts; ++ik) {
                        const auto& Pk =
                            last_.P_k_beta[static_cast<std::size_t>(ik)];
                        auto& aspc =
                            sgcp_aspc_k_beta_[static_cast<std::size_t>(ik)];
                        if (have_pred
                            && static_cast<int>(P_pred_k_beta.size())
                                   == last_.n_kpts) {
                            linalg::CplxMatrix P_corr;
                            aspc.correct(
                                Pk, P_pred_k_beta[static_cast<std::size_t>(ik)],
                                P_corr);
                            aspc.record(P_corr);
                        } else {
                            aspc.record(Pk);
                        }
                    }
                } else {
                    sgcp_aspc_k_beta_.clear();
                }
            } else {
                for (auto& a : sgcp_aspc_k_)
                    a.reset();
                for (auto& a : sgcp_aspc_k_beta_)
                    a.reset();
            }

            frame.meta.get(kEnergy) =
                last_.e_total * ATV::units::kHartreeInEv;

            if (compute_forces) {
                // Multi-k forces at the post-corrector density. Partial SCF
                // may leave converged=false; SGCP semantics are forces at
                // the corrector state (same as Γ-UKS SGCP).
                if (!(last_.n_kpts > 0
                      && static_cast<int>(last_.C_occ_k.size())
                             == last_.n_kpts
                      && static_cast<int>(last_.eps_k.size())
                             == last_.n_kpts
                      && static_cast<int>(last_.P_k.size()) == last_.n_kpts))
                    throw std::runtime_error(
                        "GPWEvaluator: multi-k second_gen_cp force path needs "
                        "P_k / C_occ_k / eps_k from the corrector SCF");
                if (uks_sgcp
                    && (static_cast<int>(last_.P_k_beta.size()) != last_.n_kpts
                        || static_cast<int>(last_.C_occ_k_beta.size())
                               != last_.n_kpts
                        || static_cast<int>(last_.eps_k_beta.size())
                               != last_.n_kpts))
                    throw std::runtime_error(
                        "GPWEvaluator: multi-k UKS second_gen_cp force path "
                        "needs dual-spin k payloads from the corrector SCF");
                GpwResult force_scf = last_;
                force_scf.converged = true;
                const auto F = gpw_forces_from_scf(atoms, *solver_, force_scf);
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

        // ── Γ real ASPC ────────────────────────────────────────────────
        if (!sgcp_aspc_)
            throw std::logic_error(
                "GPWEvaluator: second_gen_cp Γ path missing ASPC state");
        xt::xtensor<double, 2> P_pred;
        xt::xtensor<double, 2> P_pred_beta;
        const bool have_pred_a = sgcp_aspc_->predict(P_pred);
        bool have_pred_b = false;
        if (uks_sgcp && sgcp_aspc_beta_)
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

        auto record_channel = [](electronic::ASPCPredictor& aspc,
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
                // orbitals without the fixed-density payload. DIIS path can
                // re-eval with max_iter=0; OT forbids fixed-density eval —
                // OT exit must export H_ks/S itself (second-gen-cpmd-ot).
                if (last_.H_ks.size() == 0 || last_.S.size() == 0) {
                    if (use_ot || solver_->use_ot)
                        throw std::runtime_error(
                            "GPWEvaluator: second_gen_cp OT corrector did not "
                            "export H_ks/S (required for forces at density)");
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
            if (uks_sgcp && sgcp_aspc_beta_) {
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
            } else if (use_ot || solver_->use_ot) {
                // OT corrector: prefer forces-from-scf when C/eps present
                // (partial OT may leave converged=false). Fall back to
                // density path when H_ks was exported.
                if (last_.P.size() > 0 && last_.C.size() > 0
                    && !last_.eps.empty()) {
                    GpwResult force_scf = last_;
                    force_scf.converged = true;
                    F = gpw_forces_from_scf(atoms, *solver_, force_scf);
                } else if (last_.C.size() > 0 && last_.P.size() > 0
                           && last_.H_ks.size() > 0) {
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
                } else {
                    throw std::runtime_error(
                        "GPWEvaluator: second_gen_cp OT force path needs "
                        "P/C(/eps or H_ks) from the corrector SCF");
                }
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

    // ── Zero-SCF linearized XL-BOMD (shadow Born–Oppenheimer forces) ────
    // After the cold full-SCF seed (xl history ready): one H[n] build from
    // the auxiliary density n, one diagonalization → D = q[n], shadow
    // energy E = Tr(H D) − Tr(H n) + E_KS[n], forces at (D, H[n]). No
    // iterative SCF (iters = 0). Cold seed falls through to the SCF-XL
    // path below (first advance seeds history with D_scf).
    if (xl_bomd_linearized && xl_ && xl_->ready()) {
        const linalg::RealMatrix& n_aux = xl_->aux();
        const int saved_max = solver_->max_iter;
        solver_->max_iter = 0;
        try {
            last_ = solver_->run(atoms, /*guess_C=*/nullptr, &n_aux);
        } catch (...) {
            solver_->max_iter = saved_max;
            throw;
        }
        solver_->max_iter = saved_max;

        if (last_.H_ks.size() == 0 || last_.S.size() == 0
            || last_.P.size() == 0)
            throw std::runtime_error(
                "GPWEvaluator: xl_bomd_linearized fixed-density step needs "
                "H_ks, S, and P from max_iter=0");

        const int nb = static_cast<int>(last_.H_ks.shape(0));
        const int n_occ = solver_->n_occ;
        if (n_occ < 1)
            throw std::runtime_error(
                "GPWEvaluator: xl_bomd_linearized needs n_occ >= 1");

        // One diagonalization of H[n] → D = q[n] (zero iterative SCF).
        auto [eig, C] =
            linalg::diagonalize_fock_generalized(last_.H_ks, last_.S);
        const linalg::RealMatrix D =
            linalg::build_density_matrix(C, n_occ, nb, 2.0);

        // Shadow energy: Tr(H D) − Tr(H n) + E_KS[n]. Recovers E_KS at D=n.
        auto tr_HM = [nb](const linalg::RealMatrix& H,
                          const linalg::RealMatrix& M) {
            double t = 0.0;
            for (int i = 0; i < nb; ++i)
                for (int j = 0; j < nb; ++j)
                    t += H(static_cast<std::size_t>(i),
                           static_cast<std::size_t>(j))
                         * M(static_cast<std::size_t>(j),
                             static_cast<std::size_t>(i));
            return t;
        };
        const double e_ks_n = last_.e_total;
        const double tr_Hn = tr_HM(last_.H_ks, last_.P);
        const double tr_HD = tr_HM(last_.H_ks, D);
        last_.e_total = tr_HD - tr_Hn + e_ks_n;

        last_.P = D;
        last_.C = C;
        last_.eps.assign(eig.begin(), eig.end());
        last_.converged = true; // stationary for the shadow functional
        last_.iters = 0;        // zero iterative SCF

        xl_->advance(D);

        frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

        if (compute_forces) {
            linalg::RealMatrix C_occ = xt::zeros<double>(
                {static_cast<std::size_t>(nb),
                 static_cast<std::size_t>(n_occ)});
            for (int mu = 0; mu < nb; ++mu)
                for (int i = 0; i < n_occ; ++i)
                    C_occ(static_cast<std::size_t>(mu),
                          static_cast<std::size_t>(i)) = C(
                        static_cast<std::size_t>(mu),
                        static_cast<std::size_t>(i));
            // Forces of U(R,n) at fixed n: density D, KS matrix H[n].
            const auto F = gpw_forces_at_density(atoms, *solver_, D, C_occ,
                                                 last_.H_ks);
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

    // BOMD / FD-phonon warm-start: reuse previous MOs and spin densities.
    // GpwDft::run accepts optional guess_C; empty last_.C → core guess.
    // Multi-k has no real C warm-start (Bloch C_k live on the k path).
    //
    // UKS critical path: guess_C alone seeds *both* spins from α-like
    // columns (see GpwDft::run), which is a cold open-shell start. Pass
    // the converged P_α / P_β so successive force evaluations (Phonon
    // displacements, UKS AIMD) re-enter near the previous solution.
    const ATV::cpu::linalg::RealMatrix* guess = nullptr;
    if (last_.C.size() > 0 && last_.n_kpts == 0) {
        guess = &last_.C;
    }
    // XL-BOMD: propagated auxiliary density is the SCF initial-density
    // guess (Niklasson dissipative recursion, seeded on the first step).
    // Γ RKS: total P. Γ UKS: P_α / P_β. multi-k RKS: P_k. multi-k UKS:
    // dual P_k^α / P_k^β.
    // (Linearized warm path returned above when aux is ready.)
    const ATV::cpu::linalg::RealMatrix* guess_p =
        (xl_ && xl_->ready()) ? &xl_->aux() : nullptr;
    const ATV::cpu::linalg::RealMatrix* guess_p_beta =
        (xl_beta_ && xl_beta_->ready()) ? &xl_beta_->aux() : nullptr;
    if (guess_p == nullptr && last_.n_kpts == 0 && last_.P.size() > 0
        && last_.converged) {
        if (n_alpha + n_beta > 0) {
            if (last_.P_beta.size() > 0) {
                guess_p = &last_.P;
                guess_p_beta = &last_.P_beta;
            }
        } else {
            guess_p = &last_.P;
        }
    }
    const std::vector<ATV::cpu::linalg::CplxMatrix>* guess_p_k = nullptr;
    const std::vector<ATV::cpu::linalg::CplxMatrix>* guess_p_k_beta = nullptr;
    if (!xl_k_.empty() && xl_k_[0].ready()) {
        xl_guess_p_k_.resize(xl_k_.size());
        for (std::size_t ik = 0; ik < xl_k_.size(); ++ik)
            xl_guess_p_k_[ik] = xl_k_[ik].aux();
        guess_p_k = &xl_guess_p_k_;
        if (!xl_k_beta_.empty() && xl_k_beta_[0].ready()
            && xl_k_beta_.size() == xl_k_.size()) {
            xl_guess_p_k_beta_.resize(xl_k_beta_.size());
            for (std::size_t ik = 0; ik < xl_k_beta_.size(); ++ik)
                xl_guess_p_k_beta_[ik] = xl_k_beta_[ik].aux();
            guess_p_k_beta = &xl_guess_p_k_beta_;
        }
    }
    last_ = solver_->run(atoms, guess, guess_p, guess_p_beta, guess_p_k,
                         guess_p_k_beta);
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
    // Multi-k XL: advance one complex propagator per k from SCF P_k
    // (RKS total or UKS α). Multi-k UKS also advances β from P_k_beta.
    if (xl_bomd && last_.n_kpts > 0) {
        const bool uks_k = (n_alpha + n_beta > 0);
        const bool have_a =
            last_.converged
            && static_cast<int>(last_.P_k.size()) == last_.n_kpts
            && last_.n_kpts > 0;
        const bool have_b =
            !uks_k
            || (have_a
                && static_cast<int>(last_.P_k_beta.size()) == last_.n_kpts);
        if (have_a && have_b) {
            auto reseat = [&](std::vector<electronic::XlBomdCplx>& xs) {
                if (static_cast<int>(xs.size()) != last_.n_kpts) {
                    xs.clear();
                    xs.reserve(static_cast<std::size_t>(last_.n_kpts));
                    for (int ik = 0; ik < last_.n_kpts; ++ik)
                        xs.emplace_back(xl_bomd_k);
                }
            };
            reseat(xl_k_);
            for (int ik = 0; ik < last_.n_kpts; ++ik)
                xl_k_[static_cast<std::size_t>(ik)].advance(
                    last_.P_k[static_cast<std::size_t>(ik)]);
            if (uks_k) {
                reseat(xl_k_beta_);
                for (int ik = 0; ik < last_.n_kpts; ++ik)
                    xl_k_beta_[static_cast<std::size_t>(ik)].advance(
                        last_.P_k_beta[static_cast<std::size_t>(ik)]);
            } else {
                xl_k_beta_.clear();
            }
        } else {
            for (auto& xl : xl_k_)
                xl.reset();
            for (auto& xl : xl_k_beta_)
                xl.reset();
        }
    }
    frame.meta.get(kEnergy) = last_.e_total * ATV::units::kHartreeInEv;

    if (compute_forces) {
        // Prefer forces at the density just converged (no second SCF). UKS has
        // no guess_C warm-start, so a re-SCF would double wall; RKS cold
        // re-SCF is kept only when the energy result lacks orbitals (should
        // not happen after a normal converged run). Multi-k uses C_occ_k.
        std::vector<std::array<double, 3>> F;
        const bool gamma_payload =
            last_.converged && last_.P.size() > 0 && last_.C.size() > 0
            && !last_.eps.empty();
        const bool k_payload =
            last_.converged && last_.n_kpts > 0
            && static_cast<int>(last_.C_occ_k.size()) == last_.n_kpts
            && static_cast<int>(last_.eps_k.size()) == last_.n_kpts;
        if (gamma_payload || k_payload) {
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
    // Multi-k RKS: seed one CpmdCplx per sampled k from Bloch C_k / C_occ_k.
    if (cpmd && cpmd_k_.empty()
        && (k_mesh[0] != 1 || k_mesh[1] != 1 || k_mesh[2] != 1)) {
        if (!last_.converged || last_.n_kpts < 1)
            throw std::runtime_error(
                "GPWEvaluator: multi-k cpmd seed requires a converged multi-k "
                "SCF");
        const int nk = last_.n_kpts;
        const int n_occ = solver_->n_occ;
        if (n_occ < 1)
            throw std::runtime_error(
                "GPWEvaluator: multi-k RKS CPMD seed requires n_occ >= 1");
        // Prefer full C_k; fall back to flat C_occ_k.
        const bool have_C_k =
            static_cast<int>(last_.C_k.size()) == nk
            && last_.C_k[0].size() > 0;
        const bool have_C_occ =
            static_cast<int>(last_.C_occ_k.size()) == nk
            && !last_.C_occ_k[0].empty();
        if (!have_C_k && !have_C_occ)
            throw std::runtime_error(
                "GPWEvaluator: multi-k cpmd seed missing C_k / C_occ_k");
        cpmd_k_.clear();
        cpmd_k_.reserve(static_cast<std::size_t>(nk));
        for (int ik = 0; ik < nk; ++ik) {
            linalg::CplxMatrix C_occ;
            if (have_C_k) {
                const auto& Cfull = last_.C_k[static_cast<std::size_t>(ik)];
                const int nb = static_cast<int>(Cfull.shape(0));
                if (static_cast<int>(Cfull.shape(1)) < n_occ)
                    throw std::runtime_error(
                        "GPWEvaluator: multi-k cpmd seed: C_k columns < n_occ");
                C_occ = xt::zeros<std::complex<double>>(
                    {static_cast<std::size_t>(nb),
                     static_cast<std::size_t>(n_occ)});
                for (int mu = 0; mu < nb; ++mu)
                    for (int o = 0; o < n_occ; ++o)
                        C_occ(mu, o) = Cfull(mu, o);
            } else {
                const auto& flat = last_.C_occ_k[static_cast<std::size_t>(ik)];
                // Infer nb from flat size / n_occ.
                if (static_cast<int>(flat.size()) % n_occ != 0)
                    throw std::runtime_error(
                        "GPWEvaluator: multi-k cpmd seed: C_occ_k size not "
                        "divisible by n_occ");
                const int nb = static_cast<int>(flat.size()) / n_occ;
                C_occ = xt::zeros<std::complex<double>>(
                    {static_cast<std::size_t>(nb),
                     static_cast<std::size_t>(n_occ)});
                for (int mu = 0; mu < nb; ++mu)
                    for (int o = 0; o < n_occ; ++o)
                        C_occ(mu, o) =
                            flat[static_cast<std::size_t>(mu) * n_occ
                                 + static_cast<std::size_t>(o)];
            }
            electronic::CpmdCplx ch(cpmd_mu);
            ch.seed(C_occ);
            cpmd_k_.push_back(std::move(ch));
        }
        // Seed step: the electronic coefficients have not been propagated, so
        // both electronic energies are exactly 0 (a re-`build()` after a warm
        // run must reset the thermostat readback too, not leave it stale).
        last_te_ev_ = 0.0;
        last_etherm_ev_ = 0.0;
        publish_electronic_energies(frame, last_te_ev_, last_etherm_ev_);
    } else if (cpmd_ && !cpmd_->ready()) {
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
        // Seed step: see the multi-k branch above — both terms are 0 here.
        last_te_ev_ = 0.0;
        last_etherm_ev_ = 0.0;
        publish_electronic_energies(frame, last_te_ev_, last_etherm_ev_);
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
