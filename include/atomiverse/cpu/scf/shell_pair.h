/// @file cpu/scf/shell_pair.h
/// @brief ShellPairList: canonical shell pairs with Cauchy-Schwarz bounds for direct SCF.
///
/// Pre-computes Schwarz bounds Q_ij = sqrt((ij|ij)) for all canonical pairs (i>=j).
/// Used by direct-SCF Fock builders to skip negligible shell-quartet contributions.
///
/// @see build_shell_pair_list
/// @see build_g_rhf_direct, build_jk_uhf_direct, build_j_direct
#pragma once

#include <vector>
#include <array>

#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

struct IntegralPack;  // forward declare
struct ShellList;     // forward declare

/// @brief Single canonical shell pair with Cauchy-Schwarz bound.
struct ShellPair {
    int i;          ///< Basis function index (i >= j).
    int j;          ///< Basis function index.
    double schwarz; ///< Cauchy-Schwarz bound sqrt((ij|ij)).
};

/// @brief Precomputed canonical shell pair list for direct SCF.
///
/// Owns copies of shells, atoms, and shell_data so it can outlive the IntegralPack.
/// Only pairs with schwarz > threshold are retained.
struct ShellPairList {
    std::vector<ShellPair>            pairs;      ///< Screened canonical pairs.
    std::vector<native::Shell>        shells;     ///< Owned shell copies.
    std::vector<native::Atom>         atoms;      ///< Owned atom copies.
    std::vector<std::vector<double>>  shell_data; ///< Backing [α..., c...] per shell (length 2*n_prim).
    int    n_basis    = 0;                        ///< Total basis functions.
    double threshold  = 1e-10;                    ///< Schwarz screening threshold.
};

/// @brief Build the canonical shell pair list from an IntegralPack.
///
/// Computes Q_ij = sqrt((ij|ij)) for all i >= j, retaining pairs where Q_ij > thr.
/// The returned list owns all shell data (safe to use after the IntegralPack is destroyed).
///
/// @param pack Source integral pack (uses pack.shells, pack.shell_data, pack.atoms).
/// @param thr  Cauchy-Schwarz screening threshold [default 1e-10].
/// @return ShellPairList with screened pairs and owned data.
ShellPairList build_shell_pair_list(const IntegralPack& pack, double thr = 1e-10);

// ── Shell-block pair list (libcint/PySCF approach) ────────────────────

/// @brief Canonical shell-block pair with block-level Cauchy-Schwarz bound.
///
/// The Schwarz bound is Q_{AB} = sqrt(max_{ia,ib} |(ia ib|ia ib)|), covering
/// all AO combinations within the block pair.  Used to screen entire block
/// quartets before computing eri_block().
struct ShellBlockPair {
    int si;          ///< Block index i (si >= sj, canonical ordering).
    int sj;          ///< Block index j.
    double schwarz;  ///< Block-level Schwarz bound sqrt(max_diag (AB|AB)).
};

/// @brief Precomputed canonical shell-block pair list for direct SCF.
///
/// Groups AO-level shells into shell blocks (s-shells: n_cart=1;
/// p-shells: n_cart=3 grouping px/py/pz together).  For STO-3G aspirin
/// (73 AOs), this gives 47 shell blocks vs 73 AO-level shells.
///
/// The list owns all shell-block data so it can outlive the IntegralPack.
struct ShellBlockPairList {
    std::vector<ShellBlockPair>           pairs;       ///< Screened canonical block pairs.
    std::vector<native::ShellBlock>       blocks;      ///< Owned shell-block descriptors.
    std::vector<std::vector<double>>      block_data;  ///< Backing [α..., c...] per block (length 2*n_prim).
    int    n_basis   = 0;                              ///< Total AO basis functions.
    double threshold = 1e-10;                          ///< Schwarz screening threshold.
};

/// @brief Build the canonical shell-block pair list from an IntegralPack.
///
/// Groups consecutive AO-level shells into shell blocks, then computes
/// block-level Schwarz bounds Q_{AB} = sqrt(max_{ia,ib} (ia ib|ia ib))
/// using eri_block().  Pairs with Q_{AB} ≤ thr are dropped.
///
/// @param pack Source integral pack.
/// @param thr  Cauchy-Schwarz screening threshold [default 1e-10].
/// @return ShellBlockPairList with owned data.
ShellBlockPairList build_shell_block_pair_list(const IntegralPack& pack,
                                               double thr = 1e-10);

/// @brief Group a flat vector of AO-level `native::Shell` entries into
///        libcint-style l-blocks without materializing Schwarz bounds.
///
/// Consecutive shells sharing an angular momentum (same `l`) collapse
/// into one `native::ShellBlock` with `n_cart = n_cart_of(l)`. Callers
/// feed this with `ShellList::shells` (one-electron integral path) or
/// `IntegralPack::shells` (Fock-build Schwarz path) — both use the same
/// flat layout.
///
/// Returned blocks alias the input's `alpha` / `coeff` pointers, so the
/// source storage must outlive the returned vector.
[[nodiscard]] std::vector<native::ShellBlock>
blocks_from_shells(const std::vector<native::Shell>& shells);

} // namespace scf
} // namespace cpu
} // namespace ATV
