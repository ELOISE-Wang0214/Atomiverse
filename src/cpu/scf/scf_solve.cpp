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
#include "atomiverse/cpu/scf/diis.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/xc_grid_engine.h"

#include <cmath>
#include <algorithm>
#include <vector>

#include <xtensor/xnorm.hpp>
#include <xtensor/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

namespace ATV {
namespace cpu {
namespace scf {

// ── Unified SCF solver ──────────────────────────────────────────────���───

SCFResult scf_solve(const IntegralPack& integrals, const SCFConfig& config) {
    Int n = integrals.n_basis;
    RealMatrix S_cm     = integrals.S;
    RealMatrix Hcore_cm = integrals.Hcore;

    bool restricted  = (config.spin == SpinType::Restricted);
    bool is_hf       = (config.method == MeanFieldType::HF);
    int n_alpha      = config.n_alpha;
    int n_beta       = config.n_beta;

    // ── 1. Orthogonalization ────────────────────────────────────────────
    RealMatrix X = build_orthogonalization_matrix(S_cm);

    // ── 2. Grid for KS methods ──────────────────────────────────────────
    MolecularGrid grid;
    if (!is_hf && !integrals.atoms.empty()) {
        grid = build_molecular_grid(integrals.atoms.data(),
                                    static_cast<int>(integrals.atoms.size()),
                                    config.n_radial, config.n_angular);
    }

    // ── 3. Initial guess ────────────────────────────────────────────────
    auto [eig0, C0] = diagonalize_fock(Hcore_cm, X);

    RealMatrix P_alpha, P_beta;
    if (restricted) {
        // RHF: P = 2 * C_occ * C_occ^T
        P_alpha = build_density_matrix(C0, n_alpha, n, Real(2));
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
    Real damping = config.damping_factor;
    Real prev_residual = Real(1);  // track residual from previous iteration

    RealMatrix C_alpha = C0, C_beta = C0;
    RealVector eig_alpha = eig0, eig_beta = eig0;

    SCFResult result{};
    result.n_alpha = n_alpha;
    result.n_beta  = n_beta;
    result.nuclear_repulsion = integrals.V_nuc;

    // ── 5. SCF loop ─────────────────────────────────────────────────────
    for (int iter = 0; iter < config.max_iter; ++iter) {

        // ── 5a. Build Fock matrices (method-dependent) ──────────────────
        RealMatrix F_alpha, F_beta;
        RealMatrix J;
        Real E_xc_grid = 0;
        Real int_n_alpha = 0, int_n_beta = 0;

        if (restricted && is_hf) {
            // RHF: fused G = J - 0.5K
            RealMatrix G;
            build_coulomb_exchange_fused(P_alpha, integrals.ERI, n, G);
            F_alpha = Hcore_cm + G;
            F_beta  = F_alpha;  // same for restricted
        } else if (!restricted && is_hf) {
            // UHF: separate J and K per spin
            RealMatrix P_total = P_alpha + P_beta;
            RealMatrix K_alpha, K_beta;
            build_coulomb(P_total, integrals.ERI, n, J);
            build_exchange(P_alpha, integrals.ERI, n, K_alpha);
            build_exchange(P_beta,  integrals.ERI, n, K_beta);
            F_alpha = Hcore_cm + J - K_alpha;
            F_beta  = Hcore_cm + J - K_beta;
        } else {
            // UKS: J + V_xc from grid
            RealMatrix P_total = P_alpha + P_beta;
            build_coulomb(P_total, integrals.ERI, n, J);
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
            P_alpha_new = build_density_matrix(C_alpha, n_alpha, n, Real(2));
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
            P_alpha_next = damping * P_alpha_new + (Real(1) - damping) * P_alpha;
            if (!restricted)
                P_beta_next = damping * P_beta_new + (Real(1) - damping) * P_beta;
            else
                P_beta_next = P_alpha_next;
        } else {
            P_alpha_next = P_alpha_new;
            P_beta_next  = restricted ? P_alpha_new : P_beta_new;
        }

        // ── 5f. Energy evaluation ───────────────────────────────────────
        // Use the trace form: E = 0.5 * Σ_σ Tr[(Hcore + F^σ) · P^σ] + V_nuc
        // This is self-consistent with the Fock matrices that were actually built.
        Real E_total;
        Real E_one = 0, E_J = 0, E_K = 0;

        if (restricted) {
            // RHF: E = 0.5 * Tr[(Hcore + F) * P] + V_nuc
            E_total = Real(0.5) * xt::sum((Hcore_cm + F_alpha) * P_alpha_new)()
                    + integrals.V_nuc;
            // Decompose: E_one = Tr[Hcore * P], E_J = 0.5*Tr[J*P], E_K via difference
            E_one = xt::sum(Hcore_cm * P_alpha_new)();
            // For restricted HF, G = J - 0.5K and E_elec = Tr[H*P] + 0.5*Tr[G*P]
            // = E_one + 0.5*Tr[G*P]. E_J and E_K need separate J,K builds (expensive).
            // Approximate via: E_elec = E_total - V_nuc, E_J+E_K = E_elec - E_one
        } else if (is_hf) {
            // UHF: E = 0.5 * [Tr((H+F^α)*P^α) + Tr((H+F^β)*P^β)] + V_nuc
            E_total = Real(0.5) * (
                xt::sum((Hcore_cm + F_alpha) * P_alpha_new)() +
                xt::sum((Hcore_cm + F_beta)  * P_beta_new)()
            ) + integrals.V_nuc;
            E_one = xt::sum(Hcore_cm * P_alpha_new)() + xt::sum(Hcore_cm * P_beta_new)();
        } else {
            // UKS: E = Σ_σ Tr[Hcore·P^σ] + 0.5*Tr[J·P_total] + E_xc + V_nuc
            RealMatrix P_total_new = P_alpha_new + P_beta_new;
            E_one = xt::sum(Hcore_cm * P_alpha_new)() + xt::sum(Hcore_cm * P_beta_new)();
            E_J = Real(0.5) * xt::sum(J * P_total_new)();
            E_total = E_one + E_J + E_xc_grid + integrals.V_nuc;
        }

        // ── 5g. Convergence check ───────────────────────────────────────
        Real res_a = xt::norm_l2(P_alpha_next - P_alpha)();
        Real res_b = restricted ? Real(0) :
                     (n_beta > 0 ? xt::norm_l2(P_beta_next - P_beta)() : Real(0));
        Real residual = std::max(res_a, res_b);
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

        if (iter > 0 && residual < config.tol) {
            result.converged = true;
            break;
        }
    }

    // ── 6. Post-SCF: spin contamination ─────────────────────────────────
    if (!restricted) {
        result.spin_contamination = compute_spin_contamination(
            result.C_alpha, result.C_beta, S_cm, n_alpha, n_beta);
    }

    return result;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
