/// @file cpu/scf/operators.h
/// @brief Two-electron operator builders: stored-ERI (legacy) and direct-SCF variants.
///
/// **Direct-SCF functions** (preferred for CPU):
///   build_g_rhf_direct   — RHF G = J - 0.5*K, one ERI pass, 8-fold symmetry + Schwarz screening
///   build_jk_uhf_direct  — UHF J, K_alpha, K_beta in one ERI pass
///   build_j_direct       — UKS Coulomb J only
///
/// **Stored-ERI functions** (legacy, kept for CUDA/test parity):
///   build_coulomb, build_exchange, build_coulomb_exchange_fused
///
/// @see cpu::scf::scf_solve which uses the direct-SCF path internally
/// @see cpu::scf::ShellPairList for the precomputed pair list input
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/scf_utils.h"
#include "atomiverse/cpu/scf/shell_pair.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Build Coulomb matrix J from ERI and total density.
/// @details J_ij = sum_kl P_kl (ij|kl). Exploits (i,j) symmetry.
/// @param P_total Total density matrix P_alpha + P_beta [n x n].
/// @param ERI Two-electron integrals in chemist's notation [n^4].
/// @param n Basis dimension.
/// @param[out] J Coulomb matrix [n x n].
void build_coulomb(const RealMatrix& P_total,
                   const xt::xtensor<double, 4>& ERI, int n,
                   RealMatrix& J);

/// @brief Build exchange matrix K from ERI and per-spin density.
/// @details K_ij = sum_kl P_kl (ik|jl). Exploits (i,j) symmetry.
/// @param P_spin Per-spin density matrix [n x n].
/// @param ERI Two-electron integrals in chemist's notation [n^4].
/// @param n Basis dimension.
/// @param[out] K Exchange matrix [n x n].
void build_exchange(const RealMatrix& P_spin,
                    const xt::xtensor<double, 4>& ERI, int n,
                    RealMatrix& K);

/// @brief Build fused Coulomb + exchange G = J - 0.5*K (RHF shortcut).
/// @details G_ij = sum_kl P_kl [(ij|kl) - 0.5*(ik|jl)].
///          P is the total RHF density (includes factor of 2).
/// @param P_total Total RHF density matrix [n x n].
/// @param ERI Two-electron integrals in chemist's notation [n^4].
/// @param n Basis dimension.
/// @param[out] G Fused two-electron matrix [n x n].
void build_coulomb_exchange_fused(const RealMatrix& P_total,
                                  const xt::xtensor<double, 4>& ERI, int n,
                                  RealMatrix& G);

// ── Direct-SCF operators (no stored ERI tensor) ──────────────────────────

/// @brief Direct RHF G = J - 0.5*K via on-the-fly ERI with 8-fold symmetry and
///        Cauchy-Schwarz screening.
///
/// Evaluates canonical shell quartets (pair_ij >= pair_kl, i>=j, k>=l) once
/// and scatters contributions to G using full 8-fold permutation symmetry.
/// Thread-local reduction ensures race-free OpenMP parallelism.
///
/// @param sp      Shell pair list with Schwarz bounds (from build_shell_pair_list).
/// @param P_total Total RHF density matrix P = 2 * P_occ [n x n].
/// @param[out] G  Fused G = J - 0.5*K [n x n].
void build_g_rhf_direct(const ShellPairList& sp,
                         const RealMatrix& P_total,
                         RealMatrix& G);

/// @brief Direct UHF J, K_alpha, K_beta via one on-the-fly ERI pass.
///
/// Computes three matrices in a single quartet loop:
///   J       from P_total = P_alpha + P_beta  (Coulomb)
///   K_alpha from P_alpha                     (alpha exchange)
///   K_beta  from P_beta                      (beta exchange)
///
/// @param sp       Shell pair list with Schwarz bounds.
/// @param P_total  Total density P_alpha + P_beta [n x n].
/// @param P_alpha  Alpha spin density [n x n].
/// @param P_beta   Beta spin density [n x n].
/// @param[out] J       Coulomb matrix [n x n].
/// @param[out] K_alpha Alpha exchange matrix [n x n].
/// @param[out] K_beta  Beta exchange matrix [n x n].
void build_jk_uhf_direct(const ShellPairList& sp,
                          const RealMatrix& P_total,
                          const RealMatrix& P_alpha,
                          const RealMatrix& P_beta,
                          RealMatrix& J,
                          RealMatrix& K_alpha,
                          RealMatrix& K_beta);

/// @brief Direct Coulomb-only J matrix for UKS (exchange via XC functional).
///
/// @param sp      Shell pair list with Schwarz bounds.
/// @param P_total Total density P_alpha + P_beta [n x n].
/// @param[out] J  Coulomb matrix [n x n].
void build_j_direct(const ShellPairList& sp,
                     const RealMatrix& P_total,
                     RealMatrix& J);

// ── Shell-block direct-SCF operators (libcint/PySCF approach) ────────────

/// @brief Shell-block direct RHF G = J - 0.5*K (libcint/PySCF approach).
///
/// Iterates canonical shell-block pair quartets (block_left >= block_right)
/// and calls native::eri_block() to evaluate all AO ERIs within a shell-block
/// at once, amortizing Boys function and R-table computation.
///
/// Key performance gain over the AO-level build_g_rhf_direct:
///   For a (pp|pp) shell-block quartet with 3 primitives (STO-3G), this
///   evaluates 81 primitive quartets × (1 Boys+R-table per quartet + 81
///   cheap inner angular loops) vs 81×81 = 6561 full eri_prim_md() calls.
///
/// OpenMP parallelism via thread-local G reduction (same as AO-level variant).
/// Scatter uses the same 8-fold ERI symmetry with the exchange mirror-permutation
/// fix from the AO-level build_g_rhf_direct.
///
/// @param sp      Shell block pair list (from build_shell_block_pair_list).
/// @param P_total Total RHF density matrix P = 2*P_occ [n x n].
/// @param[out] G  Fused G = J - 0.5*K [n x n].
void build_g_rhf_direct_block(const ShellBlockPairList& sp,
                               const RealMatrix& P_total,
                               RealMatrix& G);

} // namespace scf
} // namespace cpu
} // namespace ATV
