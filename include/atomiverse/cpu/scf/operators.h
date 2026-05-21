/// @file cpu/scf/operators.h
/// @brief Two-electron operator builders: stored-ERI and direct-SCF variants.
///
/// **Direct-SCF functions** (default for CPU — no stored ERI tensor):
///   build_g_rhf_direct   — RHF G = J - 0.5*K, one ERI pass, 8-fold symmetry + Schwarz screening
///   build_jk_uhf_direct  — UHF J, K_alpha, K_beta in one ERI pass
///   build_j_direct       — UKS Coulomb J only
///   build_g_rhf_direct_block — shell-block RHF G (production RHF/RKS path)
///
/// **Stored-ERI function** (used when `IntegralPack::ERI` is materialized):
///   build_coulomb_exchange_fused — RHF Fock build, selected at runtime by
///   `HamiltonianBuilder` whenever a dense ERI tensor is available.
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

/// @brief Build fused Coulomb + exchange G = J - 0.5*K (RHF stored-ERI path).
/// @details G_ij = sum_kl P_kl [(ij|kl) - 0.5*(ik|jl)].
///          P is the total RHF density (includes factor of 2).
///          Selected by `HamiltonianBuilder` when `IntegralPack::ERI` holds
///          a materialized dense tensor; otherwise the direct-SCF
///          `build_g_rhf_direct_block` path runs.
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
void build_g_rhf_direct(const ShellPairList& sp, const RealMatrix& P_total,
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
void build_jk_uhf_direct(const ShellPairList& sp, const RealMatrix& P_total,
                         const RealMatrix& P_alpha, const RealMatrix& P_beta,
                         RealMatrix& J, RealMatrix& K_alpha,
                         RealMatrix& K_beta);

/// @brief Direct Coulomb-only J matrix for UKS (exchange via XC functional).
///
/// @param sp      Shell pair list with Schwarz bounds.
/// @param P_total Total density P_alpha + P_beta [n x n].
/// @param[out] J  Coulomb matrix [n x n].
void build_j_direct(const ShellPairList& sp, const RealMatrix& P_total,
                    RealMatrix& J);

// ── Shell-block direct-SCF operators (libcint/PySCF approach) ────────────

/// @brief Shell-block direct exchange-only matrix K (KS-hybrid path).
///
/// Builds the full exchange matrix `K_ij = sum_kl P_kl (ik|jl)` (NO 0.5
/// factor, no Coulomb) using the same shell-block ERI amortization as
/// `build_g_rhf_direct_block`. `K` is normalized so that the pure RHF Fock
/// is `F = Hcore + J - 0.5*K`; a KS global hybrid uses
/// `F = Hcore + J - 0.5*alpha_hf*K + Vxc`.
///
/// Used by `HamiltonianBuilder` for the B3LYP / PBE0 restricted hybrid
/// Fock build. Far faster than the AO-level `build_jk_uhf_direct` exchange
/// path on benzene-scale systems (shell-block Boys/R-table amortization).
///
/// @param sp Shell block pair list (from build_shell_block_pair_list).
/// @param P  Total RHF density matrix P = 2*P_occ [n x n].
/// @param[out] K Exchange matrix [n x n], symmetric, no 0.5 factor.
void build_k_direct_block(const ShellBlockPairList& sp, const RealMatrix& P,
                          RealMatrix& K);

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
                              const RealMatrix& P_total, RealMatrix& G);

} // namespace scf
} // namespace cpu
} // namespace ATV
