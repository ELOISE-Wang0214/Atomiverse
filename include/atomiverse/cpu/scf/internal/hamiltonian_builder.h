/// @file cpu/scf/internal/hamiltonian_builder.h
/// @brief internal::HamiltonianBuilder — F = Hcore + J − αK + V_xc.
///
/// **Internal SCF helper.** Does NOT inherit `cpu::Component`.
///
/// Encapsulates the method-dispatched Fock build (RHF / UHF / RKS / UKS)
/// + ownership of the precomputed shell-pair lists and integration grid.
/// Held by `cpu::scf::SCFSolver`; populated by `SCFSolver::build()` once
/// per molecular-shape (n_basis fixed) and reused across all
/// `SCFSolver::run()` calls and SCF iterations.
#pragma once

#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/internal/scf_types.h" // SCFConfig, MeanFieldType, SpinType
#include "atomiverse/cpu/scf/internal/shell_pair.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

struct FockResult {
    RealMatrix F_alpha;     ///< α-channel Fock matrix [n × n].
    RealMatrix F_beta;      ///< β-channel Fock matrix [n × n].
    RealMatrix J;           ///< Coulomb J built from the total density.
    double E_xc_grid = 0;   ///< Grid-integrated E_xc (KS only) [Hartree].
    double int_n_alpha = 0; ///< Grid-integrated α electron count (KS only).
    double int_n_beta = 0;  ///< Grid-integrated β electron count (KS only).
    /// @brief Exact (Hartree-Fock) exchange energy admixture for a KS
    ///        global hybrid: `alpha_hf * E_K` [Hartree]. Zero for pure
    ///        KS and for HF (HF carries its exchange inside F directly).
    ///        The Fock builder applies `alpha_hf` at a single site so this
    ///        term and the `-alpha_hf*K` Fock term use the same scaled K.
    double E_exact_exchange = 0;
};

class HamiltonianBuilder {
  public:
    HamiltonianBuilder() = default;

    /// @brief One-time setup. Builds shell-pair lists + (KS only) the
    /// integration grid. Called from `SCFSolver::build()`.
    /// @details Thin unpacker: forwards `config`'s method/grid fields to the
    ///          shared private `prepare_core`, so this and the config-free
    ///          `prepare` overload build byte-for-byte the same lists + grid.
    /// @param integrals Reference IntegralPack (ownership stays with caller).
    /// @param config    SCF configuration (drives KS-vs-HF dispatch + grid size).
    void prepare(const IntegralPack& integrals, const SCFConfig& config);

    /// @brief Build F from the current α/β densities.
    /// @details Thin unpacker: forwards `config`'s spin/method/xc fields to the
    ///          shared private `step_core`, so this and the config-free `step`
    ///          overload produce byte-for-byte the same F and energies.
    /// @param integrals Same IntegralPack passed to prepare() — caller-owned.
    /// @param config    SCF config (read-only here; ownership unchanged).
    /// @param P_alpha   α density (or full closed-shell density when restricted).
    /// @param P_beta    β density (alias of P_alpha when restricted; total = P_alpha
    ///                  for the restricted KS energy formula).
    /// @return Filled FockResult.
    FockResult step(const IntegralPack& integrals, const SCFConfig& config,
                    const RealMatrix& P_alpha, const RealMatrix& P_beta) const;

    /// @brief One-time setup from discrete method parameters — the config-free
    ///        entry used by the driver-side `FockStage`.
    /// @details Passes the discrete parameters straight through to the shared
    ///          private core `prepare_core` (the same core the config-taking
    ///          overload unpacks into), so the shell-pair lists and integration
    ///          grid are byte-for-byte the same as the shipping solver path. No
    ///          config object is constructed on this entry.
    /// @param integrals Reference IntegralPack (ownership stays with caller).
    /// @param is_hf     true → Hartree-Fock (no grid); false → Kohn-Sham
    ///                  (builds the XC integration grid).
    /// @param grid_level Named DFT integration-grid quality level (KS only).
    /// @param n_radial   Radial grid points per atom (SG0 only).
    /// @param n_angular  Angular grid points per atom (SG0 only).
    /// @param fixed_grid Optional external grid (KS only; not owned).
    void prepare(const IntegralPack& integrals, bool is_hf,
                 GridLevel grid_level, int n_radial, int n_angular,
                 const MolecularGrid* fixed_grid);

    /// @brief Build F from the current α/β densities via the config-free entry.
    /// @details Passes the discrete parameters straight through to the shared
    ///          private core `step_core` (the same core the config-taking
    ///          overload unpacks into), so the Fock matrices and energies are
    ///          byte-for-byte the same as the shipping solver path. No config
    ///          object is constructed on this entry.
    /// @param integrals Same IntegralPack passed to prepare() — caller-owned.
    /// @param restricted true → restricted spin (α mirrors β); false → open.
    /// @param is_hf      true → Hartree-Fock; false → Kohn-Sham.
    /// @param xc_type    Exchange-correlation functional (KS only).
    /// @param hybrid_alpha_override Regression override for the exact-exchange
    ///                   fraction; `< 0` uses the functional's own α_hf.
    /// @param P_alpha    α density (full density when restricted).
    /// @param P_beta     β density (alias of P_alpha when restricted).
    /// @return Filled FockResult.
    FockResult step(const IntegralPack& integrals, bool restricted, bool is_hf,
                    XCType xc_type, double hybrid_alpha_override,
                    const RealMatrix& P_alpha, const RealMatrix& P_beta) const;

  private:
    /// @brief Shared one-time-setup core over discrete parameters.
    /// @details Holds the shell-pair-list + (KS) grid build. Both the
    ///          config-taking `prepare` (which unpacks its config) and the
    ///          discrete `prepare` (which forwards directly) call this — the
    ///          single site the setup logic lives, so the two entries are
    ///          bit-identical by construction.
    /// @param integrals  Reference IntegralPack (ownership stays with caller).
    /// @param is_hf      true → Hartree-Fock (no grid); false → Kohn-Sham.
    /// @param grid_level Named DFT integration-grid quality level (KS only).
    /// @param n_radial   Radial grid points per atom (SG0 only).
    /// @param n_angular  Angular grid points per atom (SG0 only).
    /// @param fixed_grid Optional external grid (KS only; not owned).
    void prepare_core(const IntegralPack& integrals, bool is_hf,
                      GridLevel grid_level, int n_radial, int n_angular,
                      const MolecularGrid* fixed_grid);

    /// @brief Shared Fock-build core over discrete parameters.
    /// @details Holds the method-dispatched F = Hcore + J − αK + V_xc build.
    ///          Both the config-taking `step` (which unpacks its config) and the
    ///          discrete `step` (which forwards directly) call this — the
    ///          single site the Fock arithmetic lives, so the two entries are
    ///          bit-identical by construction.
    /// @param integrals  Same IntegralPack passed to prepare() — caller-owned.
    /// @param restricted true → restricted spin (α mirrors β); false → open.
    /// @param is_hf      true → Hartree-Fock; false → Kohn-Sham.
    /// @param xc_type    Exchange-correlation functional (KS only).
    /// @param hybrid_alpha_override Regression override for the exact-exchange
    ///                   fraction; `< 0` uses the functional's own α_hf.
    /// @param P_alpha    α density (full density when restricted).
    /// @param P_beta     β density (alias of P_alpha when restricted).
    /// @return Filled FockResult.
    FockResult step_core(const IntegralPack& integrals, bool restricted,
                         bool is_hf, XCType xc_type,
                         double hybrid_alpha_override,
                         const RealMatrix& P_alpha,
                         const RealMatrix& P_beta) const;

    ShellBlockPairList sp_block_; ///< RHF/RKS shell-block list (libcint-style).
    ShellPairList sp_;            ///< AO-level list for UHF/UKS J/K paths.
    MolecularGrid grid_;          ///< KS XC grid (empty for HF).
};

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
