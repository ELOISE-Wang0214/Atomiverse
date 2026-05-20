/// @file cpu/scf/scf_solver.cpp
/// @brief SCFSolver implementation — orchestrates the 5 internal helpers
///        through the SCF self-consistency cycle.
///
/// `scf_solve()` is a thin functional wrapper that constructs a transient
/// `SCFSolver` and runs it; both produce bit-identical results.

#include "atomiverse/cpu/scf/scf_solver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include <xtensor/core/xmath.hpp>
#include <xtensor/reducers/xnorm.hpp>
#include <xtensor/reducers/xreducer.hpp>
#include <xtensor/views/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/scf/minao.h"
#include "atomiverse/cpu/scf/scf_utils.h"
#include "atomiverse/core/log.h"

namespace ATV {
namespace cpu {
namespace scf {

SCFSolver::SCFSolver(SCFConfig config)
    : config_(std::move(config)),
      diis_(config_.diis_max_vectors),
      convergence_(config_.tol,
                   (config_.conv_tol_grad > 0)
                       ? config_.conv_tol_grad
                       : std::sqrt(config_.tol)) {}

void SCFSolver::build(const IntegralPack& reference) {
    // Geometry-invariant pre-allocation only. Everything else (X = S^{-1/2},
    // shell-pair lists, KS grid) is geometry-dependent — atoms move during
    // geometry optimization / AIMD, so those structures are rebuilt per
    // run() to stay consistent with the scf_solve() functional contract.
    eigen_.prepare(reference.n_basis);
    density_.prepare(reference.n_basis);
    built_ = true;
}

SCFResult& SCFSolver::run(const IntegralPack& integrals,
                          const SCFInitialGuess* guess) {
    const auto t_total = LogClock::now();
    const int n         = integrals.n_basis;
    const RealMatrix S_cm     = integrals.S;
    const RealMatrix Hcore_cm = integrals.Hcore;

    const bool restricted = (config_.spin == SpinType::Restricted);
    const bool is_hf      = (config_.method == MeanFieldType::HF);
    const int  n_alpha    = config_.n_alpha;
    const int  n_beta     = config_.n_beta;
    const char* method_label = is_hf ? "HF" : "KS";
    const char* spin_label   = restricted ? "restricted" : "unrestricted";

    ATV_INFO(
        "[scf] start method={} spin={} basis={} n_alpha={} n_beta={} tol={:.1e} max_iter={}",
        method_label, spin_label, n, n_alpha, n_beta, config_.tol, config_.max_iter);

    // Lazy build — convenience for callers that skip explicit build().
    if (!built_) build(integrals);

    // Geometry-dependent setup: orthogonalization matrix, shell-pair
    // lists, KS grid. All depend on current atom positions and so are
    // rebuilt each run(). (See build() — only n_basis-sized workspaces
    // are pre-allocated there.)
    X_ = build_orthogonalization_matrix(S_cm);
    hamiltonian_.prepare(integrals, config_);

    // Reset stateful helpers for this run.
    diis_.reset();
    convergence_.reset();

    // ── Initial guess ───────────────────────────────────────────────────
    auto [eig0, C0] = diagonalize_fock_generalized(Hcore_cm, S_cm);

    RealMatrix P_alpha, P_beta;
    const bool valid_guess_alpha = guess && guess->has_alpha
        && guess->P_alpha.dimension() == 2
        && static_cast<int>(guess->P_alpha.shape()[0]) == n
        && static_cast<int>(guess->P_alpha.shape()[1]) == n;
    const bool valid_guess_beta = guess && guess->has_beta
        && guess->P_beta.dimension() == 2
        && static_cast<int>(guess->P_beta.shape()[0]) == n
        && static_cast<int>(guess->P_beta.shape()[1]) == n;

    const bool want_minao = (config_.init_guess == InitGuess::MINAO)
                          && !valid_guess_alpha
                          && minao_supported(integrals);
    RealMatrix P_minao;
    if (want_minao) {
        const auto t_guess = LogClock::now();
        P_minao = build_minao_density(integrals);
        ATV_INFO("[scf] initial guess source=minao trace_ne={:.4f} ms={:.2f}",
                 xt::sum(P_minao * S_cm)(), elapsed_ms(t_guess));
    }

    if (restricted && valid_guess_alpha) {
        P_alpha = 0.5 * (guess->P_alpha + xt::transpose(guess->P_alpha));
        P_beta = P_alpha;
        ATV_INFO("[scf] initial guess source=density_history restricted=true");
    } else if (!restricted && valid_guess_alpha && valid_guess_beta) {
        P_alpha = 0.5 * (guess->P_alpha + xt::transpose(guess->P_alpha));
        P_beta  = 0.5 * (guess->P_beta + xt::transpose(guess->P_beta));
        ATV_INFO("[scf] initial guess source=density_history restricted=false");
    } else if (want_minao && restricted) {
        P_alpha = P_minao;
        P_beta  = P_alpha;
    } else if (want_minao) {
        P_alpha = 0.5 * P_minao;
        P_beta  = 0.5 * P_minao;
        const bool is_ks = (config_.method == MeanFieldType::KS);
        const bool skip_mixing = is_ks && (n_alpha == n_beta);
        if (config_.method == MeanFieldType::HF && n_alpha > 0 && n_beta > 0
            && !skip_mixing && n_alpha != n_beta) {
            RealMatrix C_alpha_init = C0;
            apply_homo_lumo_mixing(C_alpha_init, n_alpha, n);
            P_alpha = density_.step(C_alpha_init, n_alpha, n);
            P_beta  = density_.step(C0, n_beta, n);
        }
    } else if (restricted) {
        P_alpha = density_.step(C0, n_alpha, n, double(2));
        P_beta  = P_alpha;
        ATV_INFO("[scf] initial guess source=hcore");
    } else {
        RealMatrix C_alpha_init = C0;
        const bool is_ks = (config_.method == MeanFieldType::KS);
        const bool skip_mixing = is_ks && (n_alpha == n_beta);
        if (n_alpha > 0 && n_beta > 0 && !skip_mixing)
            apply_homo_lumo_mixing(C_alpha_init, n_alpha, n);
        P_alpha = density_.step(C_alpha_init, n_alpha, n);
        P_beta  = density_.step(C0, n_beta, n);
        ATV_INFO("[scf] initial guess source=hcore");
    }

    double damping = config_.damping_factor;

    RealMatrix C_alpha = C0, C_beta = C0;
    RealVector eig_alpha = eig0, eig_beta = eig0;

    result_ = SCFResult{};
    result_.n_alpha = n_alpha;
    result_.n_beta  = n_beta;
    result_.nuclear_repulsion = integrals.V_nuc;

    // ── SCF loop ────────────────────────────────────────────────────────
    for (int iter = 0; iter < config_.max_iter; ++iter) {
        const auto t_iter = LogClock::now();

        // 1. Build Fock (method-dispatched).
        auto fock = hamiltonian_.step(integrals, config_, P_alpha, P_beta);
        RealMatrix& F_alpha = fock.F_alpha;
        RealMatrix& F_beta  = fock.F_beta;

        // 2. DIIS error + extrapolation.
        RealMatrix err_alpha = internal::diis_error(F_alpha, P_alpha, S_cm, X_);
        RealMatrix err_beta;
        if (!restricted)
            err_beta = internal::diis_error(F_beta, P_beta, S_cm, X_);

        bool diis_active = false;
        if (config_.use_diis && iter >= config_.diis_start_iter) {
            if (restricted) {
                F_alpha = diis_.step(F_alpha, err_alpha);
                F_beta  = F_alpha;
            } else {
                auto [Fa_diis, Fb_diis] = diis_.step(F_alpha, err_alpha,
                                                     F_beta,  err_beta);
                F_alpha = Fa_diis;
                F_beta  = Fb_diis;
            }
            diis_active = (diis_.size() >= 2);
        }

        // 3. Diagonalize.
        std::tie(eig_alpha, C_alpha) = eigen_.step(F_alpha, S_cm);
        if (!restricted) {
            std::tie(eig_beta, C_beta) = eigen_.step(F_beta, S_cm);
        } else {
            C_beta = C_alpha;
            eig_beta = eig_alpha;
        }

        // 4. Build new density.
        RealMatrix P_alpha_new, P_beta_new;
        if (restricted) {
            P_alpha_new = density_.step(C_alpha, n_alpha, n, double(2));
            P_beta_new  = P_alpha_new;
        } else {
            P_alpha_new = density_.step(C_alpha, n_alpha, n);
            P_beta_new  = density_.step(C_beta,  n_beta,  n);
        }

        // 5. Damping (PySCF convention: damping = fraction of OLD retained).
        RealMatrix P_alpha_next, P_beta_next;
        if (!diis_active) {
            P_alpha_next = (double(1) - damping) * P_alpha_new + damping * P_alpha;
            if (!restricted)
                P_beta_next = (double(1) - damping) * P_beta_new + damping * P_beta;
            else
                P_beta_next = P_alpha_next;
        } else {
            P_alpha_next = P_alpha_new;
            P_beta_next  = restricted ? P_alpha_new : P_beta_new;
        }

        // 6. Energy evaluation.
        double E_total;
        double E_one = 0, E_J = 0;

        if (restricted && is_hf) {
            E_total = double(0.5) * xt::sum((Hcore_cm + F_alpha) * P_alpha)()
                    + integrals.V_nuc;
            E_one = xt::sum(Hcore_cm * P_alpha)();
        } else if (restricted) {
            E_one = xt::sum(Hcore_cm * P_alpha_new)();
            E_J   = double(0.5) * xt::sum(fock.J * P_alpha_new)();
            E_total = E_one + E_J + fock.E_xc_grid + integrals.V_nuc;
        } else if (is_hf) {
            E_total = double(0.5) * (
                xt::sum((Hcore_cm + F_alpha) * P_alpha)() +
                xt::sum((Hcore_cm + F_beta)  * P_beta)()
            ) + integrals.V_nuc;
            E_one = xt::sum(Hcore_cm * P_alpha)() + xt::sum(Hcore_cm * P_beta)();
        } else {
            RealMatrix P_total_new = P_alpha_new + P_beta_new;
            E_one   = xt::sum(Hcore_cm * P_alpha_new)() + xt::sum(Hcore_cm * P_beta_new)();
            E_J     = double(0.5) * xt::sum(fock.J * P_total_new)();
            E_total = E_one + E_J + fock.E_xc_grid + integrals.V_nuc;
        }

        // 7. Diagnostics + convergence check.
        double res_a = xt::norm_l2(P_alpha_next - P_alpha)();
        double res_b = restricted ? double(0)
                     : (n_beta > 0 ? xt::norm_l2(P_beta_next - P_beta)() : double(0));
        double residual = std::max(res_a, res_b);

        const RealMatrix* err_beta_ptr = restricted ? nullptr : &err_beta;
        auto verdict = convergence_.step(E_total, err_alpha, err_beta_ptr);

        P_alpha = P_alpha_next;
        P_beta  = P_beta_next;

        // 8. Store this iteration's result.
        result_.energy = E_total;
        result_.one_electron_energy = E_one;
        result_.coulomb_energy = E_J;
        result_.xc_energy = fock.E_xc_grid;
        result_.exchange_energy = (E_total - integrals.V_nuc) - E_one - E_J - fock.E_xc_grid;
        result_.residual = residual;
        result_.iterations = iter + 1;
        result_.eigenvalues_alpha = eig_alpha;
        result_.eigenvalues_beta  = eig_beta;
        result_.P_alpha = P_alpha;
        result_.P_beta  = restricted ? P_alpha : P_beta;
        result_.C_alpha = C_alpha;
        result_.C_beta  = restricted ? C_alpha : C_beta;
        result_.integrated_n_alpha = fock.int_n_alpha;
        result_.integrated_n_beta  = fock.int_n_beta;
        result_.trace.push_back({iter + 1, E_total, residual, diis_active});

        ATV_DEBUG(
            "[scf] iter={} energy_ha={:.10f} dE={:.3e} err={:.3e} residual={:.3e} diis={} iter_ms={:.1f}",
            iter + 1, E_total, verdict.dE, verdict.err_norm, residual,
            diis_active, elapsed_ms(t_iter));

        if (verdict.converged) {
            result_.converged = true;
            ATV_INFO(
                "[scf] converged iter={} energy_ha={:.10f} dE={:.3e} err={:.3e}",
                iter + 1, E_total, verdict.dE, verdict.err_norm);
            break;
        }
    }

    // ── Post-SCF: spin contamination ───────────────────────────────────
    if (!restricted) {
        result_.spin_contamination = compute_spin_contamination(
            result_.C_alpha, result_.C_beta, S_cm, n_alpha, n_beta);
    }

    if (!result_.converged) {
        ATV_WARN(
            "[scf] not converged method={} spin={} iterations={} energy_ha={:.10f} residual={:.3e}",
            method_label, spin_label, result_.iterations, result_.energy, result_.residual);
    }

    ATV_INFO(
        "[scf] done method={} spin={} converged={} iterations={} energy_ha={:.10f} residual={:.3e} total_ms={:.1f}",
        method_label, spin_label, result_.converged, result_.iterations,
        result_.energy, result_.residual, elapsed_ms(t_total));

    return result_;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
