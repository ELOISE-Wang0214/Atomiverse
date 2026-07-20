/// @file cpu/scf/mean_field_method.cpp
/// @brief `MeanFieldMethod<D>` template implementation — Frame API; a
///        CPUDriver<scf::Phase> loop over the sub-spec-02 stages, funneled
///        through the private drive_scf_. Explicit instantiations for the four
///        concrete methods live here; the per-method ctor defaults +
///        resolve_occupation live in rhf/uhf/rks/uks.cpp.

#include "atomiverse/cpu/scf/methods.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/reducers/xnorm.hpp>
#include <xtensor/reducers/xreducer.hpp>
#include <xtensor/views/xview.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/linalg/eigensolve.h" // CanonicalOrth, build_orthogonalization_matrix
#include "atomiverse/cpu/scf/aspc_predictor.h"
#include "atomiverse/cpu/scf/internal/fci_pack.h"
#include "atomiverse/cpu/scf/internal/gradient.h"
#include "atomiverse/cpu/scf/internal/hamiltonian_builder.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/internal/level_shift.h" // kLevelShiftSkipThreshold
#include "atomiverse/cpu/scf/internal/minao.h"
#include "atomiverse/cpu/scf/internal/mp2_pack.h"
#include "atomiverse/cpu/scf/scf_guess.h"     // apply_homo_lumo_mixing
#include "atomiverse/cpu/scf/spin_analysis.h" // compute_spin_contamination

namespace ATV {
namespace cpu {
namespace scf {

// ── PositionGuard: unwrap for the integral build, restore on exit ─────────────

template <typename D>
MeanFieldMethod<D>::PositionGuard::PositionGuard(Frame& frame, bool active)
    : frame_(frame), active_(active) {
    if (!active_)
        return;
    const int n = static_cast<int>(frame_.nlocal());
    const StorageReal* __restrict px = frame_.span(kAtom_x);
    const StorageReal* __restrict py = frame_.span(kAtom_y);
    const StorageReal* __restrict pz = frame_.span(kAtom_z);
    x_.assign(px, px + n);
    y_.assign(py, py + n);
    z_.assign(pz, pz + n);
    if (n <= 1)
        return;
    if (frame_.simbox.length(0) <= 0.0 || frame_.simbox.length(1) <= 0.0
        || frame_.simbox.length(2) <= 0.0)
        return;
    StorageReal* __restrict sp_x = frame_.span(kAtom_x);
    StorageReal* __restrict sp_y = frame_.span(kAtom_y);
    StorageReal* __restrict sp_z = frame_.span(kAtom_z);
    const double x0 = sp_x[0], y0 = sp_y[0], z0 = sp_z[0];
    for (int i = 1; i < n; ++i) {
        double dx = sp_x[i] - x0, dy = sp_y[i] - y0, dz = sp_z[i] - z0;
        frame_.simbox.min_image(dx, dy, dz);
        sp_x[i] = x0 + dx;
        sp_y[i] = y0 + dy;
        sp_z[i] = z0 + dz;
    }
}

template <typename D>
MeanFieldMethod<D>::PositionGuard::~PositionGuard() {
    if (!active_)
        return;
    StorageReal* __restrict px = frame_.span(kAtom_x);
    StorageReal* __restrict py = frame_.span(kAtom_y);
    StorageReal* __restrict pz = frame_.span(kAtom_z);
    const int n = static_cast<int>(x_.size());
    for (int i = 0; i < n; ++i) {
        px[i] = x_[static_cast<std::size_t>(i)];
        py[i] = y_[static_cast<std::size_t>(i)];
        pz[i] = z_[static_cast<std::size_t>(i)];
    }
}

// ── Lifecycle ─────────────────────────────────────────────────────────────────

template <typename D>
MeanFieldMethod<D>::MeanFieldMethod() = default;
template <typename D>
MeanFieldMethod<D>::~MeanFieldMethod() = default;
template <typename D>
MeanFieldMethod<D>::MeanFieldMethod(MeanFieldMethod&&) noexcept = default;
template <typename D>
MeanFieldMethod<D>&
MeanFieldMethod<D>::operator=(MeanFieldMethod&&) noexcept = default;

template <typename D>
void MeanFieldMethod<D>::set_predictor_(int order, int corrector_steps,
                                        double corrector_tol, bool fallback) {
    use_predictor_ = true;
    corrector_steps_ = corrector_steps;
    corrector_tol_ = corrector_tol;
    fallback_full_scf_ = fallback;
    aspc_alpha_ = std::make_unique<ASPCPredictor>(order);
    aspc_beta_ = std::make_unique<ASPCPredictor>(order);
}

// ── Stage wiring + Plan compile (lazy, once) ──────────────────────────────────

template <typename D>
void MeanFieldMethod<D>::reset_built_driver_() {
    // Drop the wired stages and every trace of their Plan / criterion
    // registration, so ensure_driver_built_ can rebind cleanly against the
    // current object. Reaches the loop-node internals directly (they are
    // protected precisely so a subclass can drive its own loop).
    fock_stage_.reset();
    mix_stage_.reset();
    diag_stage_.reset();
    density_stage_.reset();
    converged_stage_.reset();
    this->pending_.clear();
    this->owned_components_.clear();
    this->plan_ = {};
    this->criteria_.clear();
    this->criteria_owned_.clear();
    this->backend_.reset();
    driver_built_ = false;
}

template <typename D>
void MeanFieldMethod<D>::ensure_driver_built_() {
    if (driver_built_) {
        if (wired_self_ == this)
            return;
        // Relocated since it wired (a move after a run) — the stages point at
        // the moved-from object's members. Drop the stale wiring and rebind.
        reset_built_driver_();
    }

    // Wire each stage to this method's PLAINLY-NAMED members by reference (no
    // state bag); the discrete method config is captured by value now that all
    // knobs are finalized. restricted_ / is_hf_ are set by the caller
    // (drive_scf_) before this runs.
    fock_stage_ = std::make_unique<FockStage>(
        restricted_, is_hf_, xc_type_, grid_level, n_radial, n_angular,
        hybrid_alpha_override, fixed_grid, *integrals_sp_, P_alpha_, P_beta_,
        F_alpha_, F_beta_, J_, e_xc_grid_, e_exact_exchange_, int_n_alpha_,
        int_n_beta_);
    mix_stage_ = std::make_unique<MixStage>(
        restricted_, use_diis, diis_start_iter, diis_max_vectors, level_shift,
        S_, X_, P_alpha_, P_beta_, iter_, F_alpha_, F_beta_, err_alpha_,
        err_beta_);
    diag_stage_ = std::make_unique<DiagonalizeStage>(
        restricted_, F_alpha_, F_beta_, S_, X_, n_dropped_, C_alpha_, C_beta_,
        eig_alpha_, eig_beta_);
    density_stage_ = std::make_unique<DensityStage>(
        restricted_, is_hf_, damping_factor, C_alpha_, C_beta_, F_alpha_,
        F_beta_, Hcore_, J_, n_alpha_, n_beta_, n_basis_, diis_active_,
        e_xc_grid_, e_exact_exchange_, v_nuc_, P_alpha_, P_beta_, e_total_,
        e_one_, e_coulomb_);
    converged_stage_ = std::make_unique<ScfConverged>(
        tol, (conv_tol_grad > 0) ? conv_tol_grad : std::sqrt(tol), restricted_,
        e_total_, err_alpha_, err_beta_);

    this->on(ATV::scf::Phase::fock).add_ref(*fock_stage_);
    this->on(ATV::scf::Phase::mix).add_ref(*mix_stage_);
    this->on(ATV::scf::Phase::diagonalize).add_ref(*diag_stage_);
    this->on(ATV::scf::Phase::density).add_ref(*density_stage_);
    this->add_criterion_ref(*converged_stage_);

    // Whole-molecule serial evaluator — no spatial decomposition.
    this->set_serial_on_replica();
    this->build(); // CPUDriver::build — compile the Plan (un-hidden via using).
    driver_built_ = true;
    wired_self_ = this;
}

// ── Nested driver lifecycle ───────────────────────────────────────────────────

template <typename D>
void MeanFieldMethod<D>::build_on(cpu::Backend& rt) {
    // Stage wiring needs geometry (IntegralPack dimensions); deferred to the
    // first drive_scf_. Parent Runtime does not own this method's Plan.
    (void)rt;
}

template <typename D>
void MeanFieldMethod<D>::destroy_on(cpu::Backend& rt) {
    (void)rt;
    reset_built_driver_();
}

template <typename D>
int MeanFieldMethod<D>::run_nested_erased(cpu::Context& ctx, cpu::Frame& f,
                                          int budget) {
    (void)ctx;
    DriveOpts opts;
    opts.unwrap = true;
    opts.forces = compute_force;
    opts.write_energy = true;
    if (budget > 0)
        opts.max_iter = budget;
    drive_scf_(f, nullptr, opts);
    return result_sp_->iterations;
}

// ── THE funnel ────────────────────────────────────────────────────────────────

template <typename D>
SCFResult& MeanFieldMethod<D>::drive_scf_(Frame& frame,
                                          const SCFInitialGuess* guess,
                                          const DriveOpts& opts) {
    // ASPC predict / correct / record — gated only on use_predictor_ (and
    // opts.apply_aspc so the corrector/fallback recursion does not re-enter).
    SCFInitialGuess density_guess;
    const SCFInitialGuess* effective_guess = guess;
    bool predicted = false;

    if (opts.apply_aspc && use_predictor_ && aspc_alpha_ && aspc_beta_) {
        xt::xtensor<double, 2> Pa_guess, Pb_guess;
        const bool have_alpha = aspc_alpha_->predict(Pa_guess);
        const bool have_beta = aspc_beta_->predict(Pb_guess);
        if (have_alpha && have_beta) {
            density_guess.P_alpha = std::move(Pa_guess);
            density_guess.P_beta = std::move(Pb_guess);
            density_guess.has_alpha = true;
            density_guess.has_beta = true;
            effective_guess = &density_guess;
            predicted = true;
        } else if (has_last_density_ && !guess) {
            density_guess.P_alpha = last_P_alpha_;
            density_guess.P_beta = last_P_beta_;
            density_guess.has_alpha = true;
            density_guess.has_beta = true;
            effective_guess = &density_guess;
        }
    } else if (opts.apply_aspc && has_last_density_ && !guess) {
        // Warm-start from the previous outer-step density (even without ASPC).
        density_guess.P_alpha = last_P_alpha_;
        density_guess.P_beta = last_P_beta_;
        density_guess.has_alpha = true;
        density_guess.has_beta = true;
        effective_guess = &density_guess;
    }

    if (opts.apply_aspc && use_predictor_ && predicted) {
        // ASPC corrector: short SCF at loosened tol; fall back to full solve.
        DriveOpts copts = opts;
        copts.max_iter = corrector_steps_;
        copts.tol_override = corrector_tol_;
        copts.apply_aspc = false;
        drive_scf_core_(frame, effective_guess, copts);
        if (!result_sp_->converged && fallback_full_scf_) {
            DriveOpts full = opts;
            full.apply_aspc = false;
            drive_scf_core_(frame, effective_guess, full);
        }
    } else {
        DriveOpts core = opts;
        core.apply_aspc = false;
        drive_scf_core_(frame, effective_guess, core);
    }

    const SCFResult& result = *result_sp_;
    if (opts.apply_aspc && use_predictor_ && aspc_alpha_ && aspc_beta_) {
        if (predicted) {
            xt::xtensor<double, 2> Pa_corr, Pb_corr;
            aspc_alpha_->correct(result.P_alpha, density_guess.P_alpha,
                                 Pa_corr);
            aspc_beta_->correct(result.P_beta, density_guess.P_beta, Pb_corr);
            aspc_alpha_->record(Pa_corr);
            aspc_beta_->record(Pb_corr);
        } else {
            aspc_alpha_->record(result.P_alpha);
            aspc_beta_->record(result.P_beta);
        }
    }
    if (opts.apply_aspc) {
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_ = result.P_beta;
    }
    return *result_sp_;
}

template <typename D>
SCFResult& MeanFieldMethod<D>::drive_scf_core_(Frame& frame,
                                               const SCFInitialGuess* guess,
                                               const DriveOpts& opts) {
    // ── Prologue: per-run assembly ────────────────────────────────────────────
    PositionGuard guard(frame, opts.unwrap);

    *integrals_sp_ = IntegralPack{frame, basis, charge, compute_eri};
    const IntegralPack& ints = *integrals_sp_;
    resolve_occupation(ints.n_electrons);

    restricted_ = (spin_ == SpinType::Restricted);
    is_hf_ = (method_ == MeanFieldType::HF);
    const bool restricted = restricted_;
    const bool is_hf = is_hf_;
    const int n = ints.n_basis;
    const int n_alpha = n_alpha_;
    const int n_beta = n_beta_;

    // Occupation validity guard (no-fallback discipline).
    if (n_alpha < 0 || n_beta < 0 || n_alpha + n_beta != ints.n_electrons
        || (restricted && n_alpha != n_beta)) {
        throw std::invalid_argument(
            "MeanFieldMethod::run: inconsistent occupation configuration "
            "(n_electrons="
            + std::to_string(ints.n_electrons) + ", n_alpha="
            + std::to_string(n_alpha) + ", n_beta=" + std::to_string(n_beta)
            + "): requires n_alpha >= 0, n_beta >= 0, n_alpha + n_beta == "
              "n_electrons, and n_alpha == n_beta when Restricted. Check the "
              "molecular charge / spin multiplicity.");
    }

    ensure_driver_built_();

    S_ = ints.S;
    Hcore_ = ints.Hcore;
    n_basis_ = n;
    v_nuc_ = ints.V_nuc;
    diis_pairs_ = 0;
    e_xc_grid_ = 0;
    e_exact_exchange_ = 0;
    int_n_alpha_ = 0;
    int_n_beta_ = 0;

    // Orthogonalization matrix X + linear-dependence handling —
    // canonical-orth path only when columns drop.
    internal::EigenSolver eigen;
    internal::DensityBuilder density;
    n_dropped_ = 0;
    if (lindep_threshold > 0) {
        linalg::CanonicalOrth ortho(S_, lindep_threshold);
        n_dropped_ = ortho.n_dropped;
        X_ = (ortho.n_dropped > 0) ? std::move(ortho.X)
                                   : std::move(ortho.lowdin);
    } else {
        X_ = linalg::build_orthogonalization_matrix(S_);
    }

    auto diagonalize = [&](const RealMatrix& F) {
        return (n_dropped_ > 0) ? eigen.step_orthogonal(F, X_)
                                : eigen.step(F, S_);
    };

    // Initial guess, reusing diagonalize + density once.
    auto [eig0, C0] = diagonalize(Hcore_);

    RealMatrix P_alpha, P_beta;
    const bool valid_guess_alpha =
        guess && guess->has_alpha && guess->P_alpha.dimension() == 2
        && static_cast<int>(guess->P_alpha.shape()[0]) == n
        && static_cast<int>(guess->P_alpha.shape()[1]) == n;
    const bool valid_guess_beta =
        guess && guess->has_beta && guess->P_beta.dimension() == 2
        && static_cast<int>(guess->P_beta.shape()[0]) == n
        && static_cast<int>(guess->P_beta.shape()[1]) == n;
    const bool want_minao = (init_guess == InitGuess::MINAO)
                            && !valid_guess_alpha && minao_supported(ints);
    RealMatrix P_minao;
    if (want_minao)
        P_minao = build_minao_density(ints);

    if (restricted && valid_guess_alpha) {
        P_alpha = 0.5 * (guess->P_alpha + xt::transpose(guess->P_alpha));
        P_beta = P_alpha;
    } else if (!restricted && valid_guess_alpha && valid_guess_beta) {
        P_alpha = 0.5 * (guess->P_alpha + xt::transpose(guess->P_alpha));
        P_beta = 0.5 * (guess->P_beta + xt::transpose(guess->P_beta));
    } else if (want_minao && restricted) {
        P_alpha = P_minao;
        P_beta = P_alpha;
    } else if (want_minao) {
        P_alpha = 0.5 * P_minao;
        P_beta = 0.5 * P_minao;
        const bool is_ks = (method_ == MeanFieldType::KS);
        const bool skip_mixing = is_ks && (n_alpha == n_beta);
        if (method_ == MeanFieldType::HF && n_alpha > 0 && n_beta > 0
            && !skip_mixing && n_alpha != n_beta) {
            RealMatrix C_alpha_init = C0;
            apply_homo_lumo_mixing(C_alpha_init, n_alpha, n);
            P_alpha = density.step(C_alpha_init, n_alpha, n);
            P_beta = density.step(C0, n_beta, n);
        }
    } else if (restricted) {
        P_alpha = density.step(C0, n_alpha, n, double(2));
        P_beta = P_alpha;
    } else {
        RealMatrix C_alpha_init = C0;
        const bool is_ks = (method_ == MeanFieldType::KS);
        const bool skip_mixing = is_ks && (n_alpha == n_beta);
        if (n_alpha > 0 && n_beta > 0 && !skip_mixing)
            apply_homo_lumo_mixing(C_alpha_init, n_alpha, n);
        P_alpha = density.step(C_alpha_init, n_alpha, n);
        P_beta = density.step(C0, n_beta, n);
    }

    P_alpha_ = std::move(P_alpha);
    P_beta_ = std::move(P_beta);
    C_alpha_ = C0;
    C_beta_ = C0;
    eig_alpha_ = eig0;
    eig_beta_ = eig0;

    *result_sp_ = SCFResult{};
    SCFResult& res = *result_sp_;
    res.n_alpha = n_alpha;
    res.n_beta = n_beta;
    res.nuclear_repulsion = ints.V_nuc;

    mix_stage_->reset();
    converged_stage_->reset();
    // Retarget the criterion for this run: a loosened threshold for the ASPC
    // corrector (opts.tol_override), else the method's own tol — restoring the
    // old corrector_loop's tolerance behavior without a second criterion.
    const double eff_tol = (opts.tol_override > 0) ? opts.tol_override : tol;
    const double eff_grad =
        (conv_tol_grad > 0) ? conv_tol_grad : std::sqrt(eff_tol);
    converged_stage_->set_tolerances(eff_tol, eff_grad);

    // ── Loop: fock → mix → diagonalize → density → check ──────────────────────
    const int budget = (opts.max_iter > 0) ? opts.max_iter : max_iter;
    const bool shift_active =
        std::abs(level_shift) > internal::kLevelShiftSkipThreshold;

    const auto saved_window = frame.resolution_window();
    setup_pass(frame);

    int executed = 0;
    for (iter_ = 0; iter_ < budget; ++iter_) {
        // diis_active is (diis_->size() >= 2) AFTER the MixStage's step this
        // iteration. DIISMixer::step pushes one
        // pair then trims to max_vectors, so size == min(pushes, max_vectors).
        if (use_diis && iter_ >= diis_start_iter) {
            diis_pairs_ = std::min(diis_pairs_ + 1, diis_max_vectors);
            diis_active_ = (diis_pairs_ >= 2);
        } else {
            diis_active_ = false;
        }

        mutable_outer_status().reset();
        step_once(iter_);
        ++executed;
        if (poll_criteria_(iter_) != ATV::Verdict::Continue)
            break;
    }

    teardown_pass();
    frame.restore_resolution_window(saved_window);

    const bool did_converge = this->converged();

    // ── Epilogue: de-shift, spin, result fill, energy / forces ────────────────
    // De-shifted orbital energies: rebuild one UNSHIFTED Fock from the
    // converged density and diagonalize it.
    if (did_converge && shift_active) {
        internal::HamiltonianBuilder deshift;
        deshift.prepare(ints, is_hf, grid_level, n_radial, n_angular,
                        fixed_grid);
        internal::FockResult fock =
            deshift.step(ints, restricted, is_hf, xc_type_,
                         hybrid_alpha_override, P_alpha_, P_beta_);
        std::tie(eig_alpha_, C_alpha_) = diagonalize(fock.F_alpha);
        if (restricted) {
            eig_beta_ = eig_alpha_;
            C_beta_ = C_alpha_;
        } else {
            std::tie(eig_beta_, C_beta_) = diagonalize(fock.F_beta);
        }
    }

    res.converged = did_converge;
    res.iterations = executed;
    res.spin = spin_;
    res.energy = e_total_;
    res.one_electron_energy = e_one_;
    res.coulomb_energy = e_coulomb_;
    res.xc_energy = e_xc_grid_;
    // HF carries exchange inside F (recovered by subtraction); KS sets it from
    // the hybrid term so the decomposition closes (zero for pure KS).
    res.exchange_energy =
        is_hf ? ((e_total_ - v_nuc_) - e_one_ - e_coulomb_ - e_xc_grid_)
              : e_exact_exchange_;
    res.residual = converged_stage_->diis_error();
    res.eigenvalues_alpha = eig_alpha_;
    res.eigenvalues_beta = eig_beta_;
    res.P_alpha = P_alpha_;
    res.P_beta = restricted ? P_alpha_ : P_beta_;
    res.C_alpha = C_alpha_;
    res.C_beta = restricted ? C_alpha_ : C_beta_;
    res.integrated_n_alpha = int_n_alpha_;
    res.integrated_n_beta = int_n_beta_;

    if (!restricted)
        res.spin_contamination = compute_spin_contamination(
            res.C_alpha, res.C_beta, S_, n_alpha, n_beta);

    *has_result_sp_ = true;

    if (opts.forces) {
        internal::ScfGradient gradient{method_,
                                       spin_,
                                       xc_type_,
                                       grid_level,
                                       n_radial,
                                       n_angular,
                                       force_native_two_electron_deriv,
                                       xc_grad_project_translation};
        last_grad_ = gradient.compute(frame, ints, res);
        const int natm = static_cast<int>(frame.nlocal());
        AccReal* __restrict fx = frame.span(kAtom_fx);
        AccReal* __restrict fy = frame.span(kAtom_fy);
        AccReal* __restrict fz = frame.span(kAtom_fz);
        for (int A = 0; A < natm; ++A) {
            fx[A] = last_grad_.forces(A, 0);
            fy[A] = last_grad_.forces(A, 1);
            fz[A] = last_grad_.forces(A, 2);
        }
        frame.meta.get(kEnergy) = last_grad_.energy;
    } else if (opts.write_energy) {
        frame.meta.get(kEnergy) = res.energy * units::kHartreeInEv;
    }

    if (on_evaluate_)
        on_evaluate_(res);

    return res;
}

// ── Entry points ──────────────────────────────────────────────────────────────

template <typename D>
SCFResult& MeanFieldMethod<D>::run(Frame& frame, const SCFInitialGuess* guess,
                                   bool write_energy) {
    DriveOpts opts;
    opts.unwrap = true;
    opts.write_energy = write_energy;
    opts.forces = false;
    return drive_scf_(frame, guess, opts);
}

template <typename D>
SCFResult& MeanFieldMethod<D>::run(const Frame& frame,
                                   const SCFInitialGuess* guess) {
    // SCF is read-only w.r.t. the frame: the stages compute on private member
    // matrices, never frame fields; this entry unwraps nothing and writes
    // nothing back. The const_cast is therefore never a mutation of the
    // caller's object (the backend only binds the pointer and touches the
    // const-qualified resolution window).
    DriveOpts opts;
    opts.unwrap = false;
    opts.write_energy = false;
    opts.forces = false;
    return drive_scf_(const_cast<Frame&>(frame), guess, opts);
}

template <typename D>
GradientResult MeanFieldMethod<D>::forces(Frame& frame) {
    DriveOpts opts;
    opts.unwrap = false;
    opts.write_energy = false;
    opts.forces = true;
    drive_scf_(frame, nullptr, opts);
    return last_grad_;
}

template <typename D>
MP2Result MeanFieldMethod<D>::mp2(Frame& frame) {
    const bool prev = compute_eri;
    compute_eri = true;
    run(frame);
    compute_eri = prev;
    return mp2_energy(*integrals_sp_, *result_sp_);
}

template <typename D>
FCIResult MeanFieldMethod<D>::fci(Frame& frame) {
    const bool prev = compute_eri;
    compute_eri = true;
    run(frame);
    compute_eri = prev;
    return fci_energy(*integrals_sp_, *result_sp_);
}

// explicit instantiations
template class MeanFieldMethod<RHF>;
template class MeanFieldMethod<UHF>;
template class MeanFieldMethod<RKS>;
template class MeanFieldMethod<UKS>;

} // namespace scf
} // namespace cpu
} // namespace ATV
