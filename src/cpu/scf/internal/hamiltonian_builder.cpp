/// @file cpu/scf/internal/hamiltonian_builder.cpp
/// @brief Method-dispatched Fock build (RHF / UHF / RKS / UKS).

#include "atomiverse/cpu/scf/internal/hamiltonian_builder.h"

#include <xtensor/reducers/xreducer.hpp>

#include "atomiverse/core/log.h"
#include "atomiverse/cpu/scf/operators.h"
#include "atomiverse/cpu/scf/xc_grid_engine.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

// ── Shared discrete cores ────────────────────────────────────────────────────
// The one site the setup / Fock-build logic lives. Both the config-taking
// overloads (which unpack their config) and the config-free discrete overloads
// (which forward straight through) call these, so the two entries are
// bit-identical by construction — the bodies below moved here verbatim from the
// former config-taking overloads, with each `config.<field>` replaced by the
// matching discrete parameter.

void HamiltonianBuilder::prepare_core(const IntegralPack& integrals, bool is_hf,
                                      GridLevel grid_level, int n_radial,
                                      int n_angular,
                                      const MolecularGrid* fixed_grid) {
    if (!integrals.shells.empty()) {
        sp_block_ = build_shell_block_pair_list(integrals);
        sp_ = build_shell_pair_list(integrals);
    }

    if (!is_hf && !integrals.atoms.empty()) {
        if (fixed_grid != nullptr) {
            // Externally supplied grid: used by the fixed-grid FD KS
            // gradient oracle so the +/- displacements share one grid.
            grid_ = *fixed_grid;
            ATV_INFO(
                "[scf] grid build skipped: using supplied fixed grid points={}",
                grid_.points.size());
        } else {
            const auto t_grid = LogClock::now();
            ATV_INFO(
                "[scf] grid build start natm={} level={} radial={} angular={}",
                integrals.atoms.size(), static_cast<int>(grid_level), n_radial,
                n_angular);
            grid_ =
                build_molecular_grid(integrals.atoms.data(),
                                     static_cast<int>(integrals.atoms.size()),
                                     grid_level, n_radial, n_angular);
            ATV_INFO("[scf] grid build done points={} ms={:.1f}",
                     grid_.points.size(), elapsed_ms(t_grid));
        }
    }
}

FockResult HamiltonianBuilder::step_core(const IntegralPack& integrals,
                                         bool restricted, bool is_hf,
                                         XCType xc_type,
                                         double hybrid_alpha_override,
                                         const RealMatrix& P_alpha,
                                         const RealMatrix& P_beta) const {
    FockResult out;
    const RealMatrix& Hcore = integrals.Hcore;
    const int n = integrals.n_basis;

    if (restricted && is_hf) {
        // RHF: PySCF approach — store ERI once, then BLAS contraction.
        // Falls back to direct SCF when ERI is not stored.
        RealMatrix G;
        if (integrals.ERI.size() > 0) {
            build_coulomb_exchange_fused(P_alpha, integrals.ERI, n, G);
        } else {
            build_g_rhf_direct_block(sp_block_, P_alpha, G);
        }
        out.F_alpha = Hcore + G;
        out.F_beta = out.F_alpha;
    } else if (!restricted && is_hf) {
        RealMatrix P_total = P_alpha + P_beta;
        RealMatrix K_alpha, K_beta;
        build_jk_uhf_direct(sp_, P_total, P_alpha, P_beta, out.J, K_alpha,
                            K_beta);
        out.F_alpha = Hcore + out.J - K_alpha;
        out.F_beta = Hcore + out.J - K_beta;
    } else {
        // KS branch (restricted or unrestricted).
        // Restricted convention: P_alpha is the FULL double-occupied density
        //                        (P_beta aliased). XC takes spin-resolved halves.
        // Unrestricted convention: P_alpha and P_beta each carry one spin.
        RealMatrix P_total_J, P_alpha_xc, P_beta_xc;
        if (restricted) {
            P_total_J = P_alpha;
            P_alpha_xc = double(0.5) * P_alpha;
            P_beta_xc = P_alpha_xc;
        } else {
            P_total_J = P_alpha + P_beta;
            P_alpha_xc = P_alpha;
            P_beta_xc = P_beta;
        }

        // Hybrid exact-exchange admixture. alpha_hf is taken from the
        // functional unless the regression override is set (>= 0). It is
        // applied at exactly ONE site below: the bare K matrix is scaled
        // by alpha_hf once, and both the Fock term (-alpha_hf*K) and the
        // energy term (alpha_hf*E_K) consume that single scaled matrix.
        const double alpha_hf = (hybrid_alpha_override >= double(0))
                                    ? hybrid_alpha_override
                                    : hybrid_exchange_fraction(xc_type);
        const bool is_hybrid =
            (alpha_hf > double(0)) || (hybrid_alpha_override >= double(0));

        auto xc = evaluate_xc_batched(grid_, integrals.shells, n, P_alpha_xc,
                                      P_beta_xc, xc_type);

        if (is_hybrid && restricted) {
            // Restricted KS global hybrid.
            //
            // Convention pin: build_k_direct_block returns the full
            // exchange matrix K with the RHF normalization
            // F = Hcore + J - 0.5*K. The "bare K" (the matrix that, as
            // F = Hcore + J - K, reproduces pure HF) is K_bare = 0.5*K.
            // alpha_hf is applied at exactly ONE site — the scalar
            // 0.5*alpha_hf multiplying K — shared by the Fock term and the
            // energy term (E_K = -0.25*Tr[P*K]).
            const double half_alpha = double(0.5) * alpha_hf;

            if (half_alpha == double(0)) {
                // alpha_hf == 0 regression path (ac-004). J is built with
                // the SAME build_j_direct call as the pure-KS branch and
                // the Fock / energy arithmetic is the byte-for-byte
                // identical expression — so the hybrid path collapses
                // exactly onto the pre-existing pure-DFT result. K is still
                // built (the hybrid branch is genuinely exercised) but the
                // 0.5*0*K term is omitted to avoid any -0.0/+0.0 artifact.
                build_j_direct(sp_, P_total_J, out.J);
                RealMatrix K;
                build_k_direct_block(sp_block_, P_total_J, K);
                out.F_alpha = Hcore + out.J + xc.Vxc_alpha;
                out.F_beta = out.F_alpha;
                out.E_exact_exchange = double(0);
            } else {
                // Production B3LYP / PBE0. Both J and K come from the fast
                // shell-block builders: build_g_rhf_direct_block gives
                // G = J - 0.5*K, build_k_direct_block gives K, hence
                // J = G + 0.5*K. The Fock matrix is then
                //   F = Hcore + J - 0.5*alpha_hf*K + Vxc
                //     = Hcore + G + 0.5*(1 - alpha_hf)*K + Vxc.
                RealMatrix G, K;
                build_g_rhf_direct_block(sp_block_, P_total_J, G);
                build_k_direct_block(sp_block_, P_total_J, K);
                out.J = G + double(0.5) * K; // Coulomb-only, for E_J
                out.F_alpha = Hcore + G
                              + double(0.5) * (double(1) - alpha_hf) * K
                              + xc.Vxc_alpha;
                out.F_beta = out.F_alpha;
                // E_K = -0.5*Tr[P*K_bare] = -0.25*Tr[P*K]; the hybrid
                // contribution to E_total is alpha_hf*E_K.
                const double E_K = double(-0.25) * xt::sum(P_total_J * K)();
                out.E_exact_exchange = alpha_hf * E_K;
            }
        } else if (is_hybrid) {
            // Unrestricted KS hybrid (open-shell). Out of the closed-shell
            // deliverable's scope, but kept correct: per-spin exchange via
            // the AO-level UHF builder. F^sigma = Hcore + J - alpha_hf*K^sigma.
            RealMatrix K_alpha, K_beta;
            build_jk_uhf_direct(sp_, P_total_J, P_alpha_xc, P_beta_xc, out.J,
                                K_alpha, K_beta);
            out.F_alpha = Hcore + out.J - alpha_hf * K_alpha + xc.Vxc_alpha;
            out.F_beta = Hcore + out.J - alpha_hf * K_beta + xc.Vxc_beta;
            const double E_K = double(-0.5)
                               * (xt::sum(P_alpha_xc * K_alpha)()
                                  + xt::sum(P_beta_xc * K_beta)());
            out.E_exact_exchange = alpha_hf * E_K;
        } else {
            // Pure KS: Coulomb only, no exchange matrix. Byte-for-byte
            // identical to the pre-hybrid path.
            build_j_direct(sp_, P_total_J, out.J);
            out.F_alpha = Hcore + out.J + xc.Vxc_alpha;
            out.F_beta = Hcore + out.J + xc.Vxc_beta;
        }

        out.E_xc_grid = xc.E_xc;
        out.int_n_alpha = xc.n_alpha_integrated;
        out.int_n_beta = xc.n_beta_integrated;
    }

    return out;
}

// ── Config-taking overloads: thin unpackers over the shared core ─────────────

void HamiltonianBuilder::prepare(const IntegralPack& integrals,
                                 const SCFConfig& config) {
    prepare_core(integrals, config.method == MeanFieldType::HF,
                 config.grid_level, config.n_radial, config.n_angular,
                 config.fixed_grid);
}

FockResult HamiltonianBuilder::step(const IntegralPack& integrals,
                                    const SCFConfig& config,
                                    const RealMatrix& P_alpha,
                                    const RealMatrix& P_beta) const {
    return step_core(integrals, config.spin == SpinType::Restricted,
                     config.method == MeanFieldType::HF, config.xc_type,
                     config.hybrid_alpha_override, P_alpha, P_beta);
}

// ── Config-free discrete overloads: pass straight through to the shared core ─

void HamiltonianBuilder::prepare(const IntegralPack& integrals, bool is_hf,
                                 GridLevel grid_level, int n_radial,
                                 int n_angular,
                                 const MolecularGrid* fixed_grid) {
    prepare_core(integrals, is_hf, grid_level, n_radial, n_angular, fixed_grid);
}

FockResult HamiltonianBuilder::step(const IntegralPack& integrals,
                                    bool restricted, bool is_hf, XCType xc_type,
                                    double hybrid_alpha_override,
                                    const RealMatrix& P_alpha,
                                    const RealMatrix& P_beta) const {
    return step_core(integrals, restricted, is_hf, xc_type,
                     hybrid_alpha_override, P_alpha, P_beta);
}

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
