/// @file cpu/scf/internal/hamiltonian_builder.cpp
/// @brief Method-dispatched Fock build (RHF / UHF / RKS / UKS).

#include "atomiverse/cpu/scf/internal/hamiltonian_builder.h"

#include "atomiverse/cpu/scf/operators.h"
#include "atomiverse/cpu/scf/xc_grid_engine.h"
#include "atomiverse/core/log.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

void HamiltonianBuilder::prepare(const IntegralPack& integrals,
                                  const SCFConfig& config) {
    if (!integrals.shells.empty()) {
        sp_block_ = build_shell_block_pair_list(integrals);
        sp_       = build_shell_pair_list(integrals);
    }

    const bool is_hf = (config.method == MeanFieldType::HF);
    if (!is_hf && !integrals.atoms.empty()) {
        const auto t_grid = LogClock::now();
        ATV_INFO(
            "[scf] grid build start natm={} radial={} angular={}",
            integrals.atoms.size(), config.n_radial, config.n_angular);
        grid_ = build_molecular_grid(integrals.atoms.data(),
                                     static_cast<int>(integrals.atoms.size()),
                                     config.n_radial, config.n_angular);
        ATV_INFO(
            "[scf] grid build done points={} ms={:.1f}",
            grid_.points.size(), elapsed_ms(t_grid));
    }

    prepared_ = true;
}

FockResult HamiltonianBuilder::step(const IntegralPack& integrals,
                                     const SCFConfig& config,
                                     const RealMatrix& P_alpha,
                                     const RealMatrix& P_beta) const {
    FockResult out;
    const bool restricted = (config.spin == SpinType::Restricted);
    const bool is_hf      = (config.method == MeanFieldType::HF);
    const RealMatrix& Hcore = integrals.Hcore;
    const int n             = integrals.n_basis;

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
        out.F_beta  = out.F_alpha;
    } else if (!restricted && is_hf) {
        RealMatrix P_total = P_alpha + P_beta;
        RealMatrix K_alpha, K_beta;
        build_jk_uhf_direct(sp_, P_total, P_alpha, P_beta,
                            out.J, K_alpha, K_beta);
        out.F_alpha = Hcore + out.J - K_alpha;
        out.F_beta  = Hcore + out.J - K_beta;
    } else {
        // KS branch (restricted or unrestricted).
        // Restricted convention: P_alpha is the FULL double-occupied density
        //                        (P_beta aliased). XC takes spin-resolved halves.
        // Unrestricted convention: P_alpha and P_beta each carry one spin.
        RealMatrix P_total_J, P_alpha_xc, P_beta_xc;
        if (restricted) {
            P_total_J  = P_alpha;
            P_alpha_xc = double(0.5) * P_alpha;
            P_beta_xc  = P_alpha_xc;
        } else {
            P_total_J  = P_alpha + P_beta;
            P_alpha_xc = P_alpha;
            P_beta_xc  = P_beta;
        }
        build_j_direct(sp_, P_total_J, out.J);
        auto xc = evaluate_xc_batched(grid_, integrals.shells, n,
                                      P_alpha_xc, P_beta_xc, config.xc_type);
        out.F_alpha       = Hcore + out.J + xc.Vxc_alpha;
        out.F_beta        = Hcore + out.J + xc.Vxc_beta;
        out.E_xc_grid     = xc.E_xc;
        out.int_n_alpha   = xc.n_alpha_integrated;
        out.int_n_beta    = xc.n_beta_integrated;
    }

    return out;
}

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
