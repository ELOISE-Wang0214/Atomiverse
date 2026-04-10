/// @file cpu/scf/scf_solve.cpp
/// @brief Unified mean-field SCF solver: RHF, UHF, UKS through one engine.
///
/// All three methods share the same SCF loop structure:
///   1. Build orthogonalization matrix X
///   2. Initial guess from Hcore (+ HOMO-LUMO mixing for unrestricted)
///   3. SCF loop:
///      a. Build Fock: method-dependent (G for HF, J+Vxc for KS)
///      b. DIIS extrapolation (combined alpha/beta)
///      c. Diagonalize: per spin channel
///      d. Build density: per spin channel
///      e. Energy evaluation + convergence check
///   4. Post-SCF: spin contamination
///
/// The only method-dependent part is the Fock build (step 3a).
/// Everything else is shared code from scf_utils.h, operators.h, diis.h.

#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/scf_utils.h"
#include "atomiverse/cpu/scf/operators.h"
#include "atomiverse/cpu/scf/shell_pair.h"
#include "atomiverse/cpu/scf/diis.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/xc_grid_engine.h"
#include "atomiverse/core/log.h"

#include <cmath>
#include <algorithm>
#include <vector>

#include <xtensor/reducers/xnorm.hpp>
#include <xtensor/views/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

namespace ATV {
namespace cpu {
namespace scf {

// ── Unified SCF solver ──────────────────────────────────────────────���───

SCFResult scf_solve(const IntegralPack& integrals, const SCFConfig& config,
                    const SCFInitialGuess* guess) {
    const auto t_total = LogClock::now();
    int n = integrals.n_basis;
    RealMatrix S_cm     = integrals.S;
    RealMatrix Hcore_cm = integrals.Hcore;

    bool restricted  = (config.spin == SpinType::Restricted);
    bool is_hf       = (config.method == MeanFieldType::HF);
    int n_alpha      = config.n_alpha;
    int n_beta       = config.n_beta;
    const char* method_label = is_hf ? "HF" : "KS";
    const char* spin_label = restricted ? "restricted" : "unrestricted";

    ATV_INFO(
        "[scf] start method={} spin={} basis={} n_alpha={} n_beta={} tol={:.1e} max_iter={}",
        method_label, spin_label, n, n_alpha, n_beta, config.tol, config.max_iter);

    // ── 0. Build shell pair lists for direct SCF (once per call, reused across iters) ──
    // Shell-block list (libcint/PySCF approach): amortizes Boys+R-table over
    // all angular combinations per primitive quartet.  Used for RHF G-build.
    ShellBlockPairList sp_block;
    ShellPairList sp;  // AO-level list kept for UHF/UKS paths
    if (!integrals.shells.empty()) {
        sp_block = build_shell_block_pair_list(integrals);
        sp = build_shell_pair_list(integrals);
    }

    // ── 1. Orthogonalization ────────────────────────────────────────────
    RealMatrix X = build_orthogonalization_matrix(S_cm);

    // ── 2. Grid for KS methods ──────────────────────────────────────────
    MolecularGrid grid;
    if (!is_hf && !integrals.atoms.empty()) {
        const auto t_grid = LogClock::now();
        ATV_INFO(
            "[scf] grid build start natm={} radial={} angular={}",
            integrals.atoms.size(), config.n_radial, config.n_angular);
        grid = build_molecular_grid(integrals.atoms.data(),
                                    static_cast<int>(integrals.atoms.size()),
                                    config.n_radial, config.n_angular);
        ATV_INFO(
            "[scf] grid build done points={} ms={:.1f}",
            grid.points.size(), elapsed_ms(t_grid));
    }

    // ── 3. Initial guess ────────────────────────────────────────────────
    auto [eig0, C0] = diagonalize_fock(Hcore_cm, X);

    RealMatrix P_alpha, P_beta;
    const bool valid_guess_alpha = guess && guess->has_alpha
        && guess->P_alpha.dimension() == 2
        && static_cast<int>(guess->P_alpha.shape()[0]) == n
        && static_cast<int>(guess->P_alpha.shape()[1]) == n;
    const bool valid_guess_beta = guess && guess->has_beta
        && guess->P_beta.dimension() == 2
        && static_cast<int>(guess->P_beta.shape()[0]) == n
        && static_cast<int>(guess->P_beta.shape()[1]) == n;

    if (restricted && valid_guess_alpha) {
        P_alpha = 0.5 * (guess->P_alpha + xt::transpose(guess->P_alpha));
        P_beta = P_alpha;
        ATV_INFO("[scf] initial guess source=density_history restricted=true");
    } else if (!restricted && valid_guess_alpha && valid_guess_beta) {
        P_alpha = 0.5 * (guess->P_alpha + xt::transpose(guess->P_alpha));
        P_beta = 0.5 * (guess->P_beta + xt::transpose(guess->P_beta));
        ATV_INFO("[scf] initial guess source=density_history restricted=false");
    } else if (restricted) {
        // RHF: P = 2 * C_occ * C_occ^T
        P_alpha = build_density_matrix(C0, n_alpha, n, double(2));
        P_beta  = P_alpha;  // alias — not used independently
    } else {
        // UHF: always apply HOMO-LUMO mixing to break spin symmetry
        // (needed to find broken-symmetry solutions at stretched geometries).
        // UKS: skip mixing for closed-shell singlets — LSDA has no
        // Coulson-Fischer point for H2/STO-3G, so the symmetric solution
        // is correct. Mixing creates a pathological initial guess.
        RealMatrix C_alpha_init = C0;
        bool is_ks = (config.method == MeanFieldType::KS);
        bool skip_mixing = is_ks && (n_alpha == n_beta);
        if (n_alpha > 0 && n_beta > 0 && !skip_mixing)
            apply_homo_lumo_mixing(C_alpha_init, n_alpha, n);
        P_alpha = build_density_matrix(C_alpha_init, n_alpha, n);
        P_beta  = build_density_matrix(C0, n_beta, n);
    }

    // ── 4. Convergence infrastructure ───────────────────────────────────
    DIIS diis(config.diis_max_vectors);
    double damping = config.damping_factor;
    double prev_residual = double(1);  // track residual from previous iteration

    RealMatrix C_alpha = C0, C_beta = C0;
    RealVector eig_alpha = eig0, eig_beta = eig0;

    SCFResult result{};
    result.n_alpha = n_alpha;
    result.n_beta  = n_beta;
    result.nuclear_repulsion = integrals.V_nuc;

    // ── 5. SCF loop ─────────────────────────────────────────────────────
    for (int iter = 0; iter < config.max_iter; ++iter) {
        const auto t_iter = LogClock::now();

        // ── 5a. Build Fock matrices (method-dependent) ──────────────────
        RealMatrix F_alpha, F_beta;
        RealMatrix J;
        double E_xc_grid = 0;
        double int_n_alpha = 0, int_n_beta = 0;

        if (restricted && is_hf) {
            // RHF: PySCF approach — store ERI once, then BLAS contraction.
            // For aspirin STO-3G (73 AOs) the ERI tensor is ~28 MB, fits in
            // L3 cache. Each iteration does one DGEMV instead of 459K shell
            // quartets.  Falls back to direct SCF when ERI is not stored
            // (large systems where n^4 doesn't fit in memory).
            RealMatrix G;
            if (integrals.ERI.size() > 0) {
                build_coulomb_exchange_fused(P_alpha, integrals.ERI, n, G);
            } else {
                build_g_rhf_direct_block(sp_block, P_alpha, G);
            }
            F_alpha = Hcore_cm + G;
            F_beta  = F_alpha;
        } else if (!restricted && is_hf) {
            // UHF: J, K_alpha, K_beta from one direct-SCF ERI pass
            RealMatrix P_total = P_alpha + P_beta;
            RealMatrix K_alpha, K_beta;
            build_jk_uhf_direct(sp, P_total, P_alpha, P_beta, J, K_alpha, K_beta);
            F_alpha = Hcore_cm + J - K_alpha;
            F_beta  = Hcore_cm + J - K_beta;
        } else {
            // UKS: Coulomb J via direct SCF; exchange handled by XC functional
            RealMatrix P_total = P_alpha + P_beta;
            build_j_direct(sp, P_total, J);
            auto xc_result = evaluate_xc_batched(grid, integrals.shells, n,
                                                  P_alpha, P_beta, config.xc_type);
            F_alpha = Hcore_cm + J + xc_result.Vxc_alpha;
            F_beta  = Hcore_cm + J + xc_result.Vxc_beta;
            E_xc_grid = xc_result.E_xc;
            int_n_alpha = xc_result.n_alpha_integrated;
            int_n_beta  = xc_result.n_beta_integrated;
        }

        // ── 5b. DIIS extrapolation ──────────────────────────────────────
        bool diis_active = false;
        if (config.use_diis && iter >= config.diis_start_iter) {
            if (restricted) {
                RealMatrix e = diis_error(F_alpha, P_alpha, S_cm);
                F_alpha = diis.extrapolate(F_alpha, e);
                F_beta = F_alpha;
            } else {
                RealMatrix ea = diis_error(F_alpha, P_alpha, S_cm);
                RealMatrix eb = diis_error(F_beta,  P_beta,  S_cm);
                auto [Fa_diis, Fb_diis] = diis.extrapolate(F_alpha, ea, F_beta, eb);
                F_alpha = Fa_diis;
                F_beta  = Fb_diis;
            }
            diis_active = (diis.size() >= 2);
        }

        // ── 5c. Diagonalize ─────────────────────────────────────────────
        auto [ea_new, Ca_new] = diagonalize_fock(F_alpha, X);
        C_alpha = Ca_new; eig_alpha = ea_new;
        if (!restricted) {
            auto [eb_new, Cb_new] = diagonalize_fock(F_beta, X);
            C_beta = Cb_new; eig_beta = eb_new;
        } else {
            C_beta = C_alpha; eig_beta = eig_alpha;
        }

        // ── 5d. Build density ───────────────────────────────────────────
        RealMatrix P_alpha_new, P_beta_new;
        if (restricted) {
            P_alpha_new = build_density_matrix(C_alpha, n_alpha, n, double(2));
            P_beta_new  = P_alpha_new;
        } else {
            P_alpha_new = build_density_matrix(C_alpha, n_alpha, n);
            P_beta_new  = build_density_matrix(C_beta,  n_beta,  n);
        }

        // ── 5e. Apply damping ────────────────────────────────────────────
        // Pre-DIIS: linear mixing to stabilize the first few iterations.
        // Post-DIIS: no density damping. DIIS on the Fock matrix handles
        // convergence. Density damping with DIIS creates spurious fixed
        // points (Pulay, J. Comp. Chem. 3, 556, 1982).
        RealMatrix P_alpha_next, P_beta_next;
        if (!diis_active) {
            P_alpha_next = damping * P_alpha_new + (double(1) - damping) * P_alpha;
            if (!restricted)
                P_beta_next = damping * P_beta_new + (double(1) - damping) * P_beta;
            else
                P_beta_next = P_alpha_next;
        } else {
            P_alpha_next = P_alpha_new;
            P_beta_next  = restricted ? P_alpha_new : P_beta_new;
        }

        // ── 5f. Energy evaluation ───────────────────────────────────────
        // Use the trace form: E = 0.5 * Σ_σ Tr[(Hcore + F^σ) · P^σ] + V_nuc
        // This is self-consistent with the Fock matrices that were actually built.
        double E_total;
        double E_one = 0, E_J = 0, E_K = 0;

        if (restricted) {
            // RHF: E = 0.5 * Tr[(Hcore + F) * P] + V_nuc
            E_total = double(0.5) * xt::sum((Hcore_cm + F_alpha) * P_alpha_new)()
                    + integrals.V_nuc;
            // Decompose: E_one = Tr[Hcore * P], E_J = 0.5*Tr[J*P], E_K via difference
            E_one = xt::sum(Hcore_cm * P_alpha_new)();
            // For restricted HF, G = J - 0.5K and E_elec = Tr[H*P] + 0.5*Tr[G*P]
            // = E_one + 0.5*Tr[G*P]. E_J and E_K need separate J,K builds (expensive).
            // Approximate via: E_elec = E_total - V_nuc, E_J+E_K = E_elec - E_one
        } else if (is_hf) {
            // UHF: E = 0.5 * [Tr((H+F^α)*P^α) + Tr((H+F^β)*P^β)] + V_nuc
            E_total = double(0.5) * (
                xt::sum((Hcore_cm + F_alpha) * P_alpha_new)() +
                xt::sum((Hcore_cm + F_beta)  * P_beta_new)()
            ) + integrals.V_nuc;
            E_one = xt::sum(Hcore_cm * P_alpha_new)() + xt::sum(Hcore_cm * P_beta_new)();
        } else {
            // UKS: E = Σ_σ Tr[Hcore·P^σ] + 0.5*Tr[J·P_total] + E_xc + V_nuc
            RealMatrix P_total_new = P_alpha_new + P_beta_new;
            E_one = xt::sum(Hcore_cm * P_alpha_new)() + xt::sum(Hcore_cm * P_beta_new)();
            E_J = double(0.5) * xt::sum(J * P_total_new)();
            E_total = E_one + E_J + E_xc_grid + integrals.V_nuc;
        }

        // ── 5g. Convergence check ───────────────────────────────────────
        double res_a = xt::norm_l2(P_alpha_next - P_alpha)();
        double res_b = restricted ? double(0) :
                     (n_beta > 0 ? xt::norm_l2(P_beta_next - P_beta)() : double(0));
        double residual = std::max(res_a, res_b);
        prev_residual = residual;

        P_alpha = P_alpha_next;
        P_beta  = P_beta_next;

        // ── 5h. Store results ───────────────────────────────────────────
        result.energy = E_total;
        result.one_electron_energy = E_one;
        result.coulomb_energy = E_J;
        result.xc_energy = E_xc_grid;
        result.exchange_energy = (E_total - integrals.V_nuc) - E_one - E_J - E_xc_grid;
        result.residual = residual;
        result.iterations = iter + 1;
        result.eigenvalues_alpha = eig_alpha;
        result.eigenvalues_beta  = eig_beta;
        result.P_alpha = P_alpha;
        result.P_beta  = restricted ? P_alpha : P_beta;
        result.C_alpha = C_alpha;
        result.C_beta  = restricted ? C_alpha : C_beta;
        result.integrated_n_alpha = int_n_alpha;
        result.integrated_n_beta  = int_n_beta;
        result.trace.push_back({iter + 1, E_total, residual, diis_active});

        ATV_DEBUG(
            "[scf] iter={} energy_ha={:.10f} residual={:.3e} diis={} iter_ms={:.1f}",
            iter + 1, E_total, residual, diis_active, elapsed_ms(t_iter));

        if (iter > 0 && residual < config.tol) {
            result.converged = true;
            ATV_INFO(
                "[scf] converged iter={} energy_ha={:.10f} residual={:.3e}",
                iter + 1, E_total, residual);
            break;
        }
    }

    // ── 6. Post-SCF: spin contamination ─────────────────────────────────
    if (!restricted) {
        result.spin_contamination = compute_spin_contamination(
            result.C_alpha, result.C_beta, S_cm, n_alpha, n_beta);
    }

    if (!result.converged) {
        ATV_WARN(
            "[scf] not converged method={} spin={} iterations={} energy_ha={:.10f} residual={:.3e}",
            method_label, spin_label, result.iterations, result.energy, result.residual);
    }

    ATV_INFO(
        "[scf] done method={} spin={} converged={} iterations={} energy_ha={:.10f} residual={:.3e} total_ms={:.1f}",
        method_label, spin_label, result.converged, result.iterations,
        result.energy, result.residual, elapsed_ms(t_total));

    return result;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
