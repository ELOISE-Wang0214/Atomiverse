/// @file cpu/scf/internal/hamiltonian_builder.h
/// @brief internal::HamiltonianBuilder — F = Hcore + J − αK + V_xc.
///
/// **Internal SCF helper.** Does NOT inherit `cpu::ComponentBase`.
///
/// Encapsulates the method-dispatched Fock build (RHF / UHF / RKS / UKS)
/// + ownership of the precomputed shell-pair lists and integration grid.
/// Held by `cpu::scf::SCFSolver`; populated by `SCFSolver::build()` once
/// per molecular-shape (n_basis fixed) and reused across all
/// `SCFSolver::run()` calls and SCF iterations.
#pragma once

#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"   // SCFConfig, MeanFieldType, SpinType
#include "atomiverse/cpu/scf/scf_utils.h"
#include "atomiverse/cpu/scf/shell_pair.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

struct FockResult {
    RealMatrix F_alpha;          ///< α-channel Fock matrix [n × n].
    RealMatrix F_beta;           ///< β-channel Fock matrix [n × n].
    RealMatrix J;                ///< Coulomb J built from the total density.
    double     E_xc_grid     = 0;///< Grid-integrated E_xc (KS only).
    double     int_n_alpha   = 0;///< Grid-integrated α electron count (KS only).
    double     int_n_beta    = 0;///< Grid-integrated β electron count (KS only).
};

class HamiltonianBuilder {
public:
    HamiltonianBuilder() = default;

    /// @brief One-time setup. Builds shell-pair lists + (KS only) the
    /// integration grid. Called from `SCFSolver::build()`.
    /// @param integrals Reference IntegralPack (ownership stays with caller).
    /// @param config    SCF configuration (drives KS-vs-HF dispatch + grid size).
    void prepare(const IntegralPack& integrals, const SCFConfig& config);

    /// @brief Build F from the current α/β densities.
    /// @param integrals Same IntegralPack passed to prepare() — caller-owned.
    /// @param config    SCF config (read-only here; ownership unchanged).
    /// @param P_alpha   α density (or full closed-shell density when restricted).
    /// @param P_beta    β density (alias of P_alpha when restricted; total = P_alpha
    ///                  for the restricted KS energy formula).
    /// @return Filled FockResult.
    FockResult step(const IntegralPack& integrals,
                    const SCFConfig& config,
                    const RealMatrix& P_alpha,
                    const RealMatrix& P_beta) const;

private:
    ShellBlockPairList sp_block_;  ///< RHF/RKS shell-block list (libcint-style).
    ShellPairList      sp_;        ///< AO-level list for UHF/UKS J/K paths.
    MolecularGrid      grid_;      ///< KS XC grid (empty for HF).
    bool               prepared_ = false;
};

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
