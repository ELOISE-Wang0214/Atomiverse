/// @file cpu/scf/uhf.h
/// @brief Unrestricted Hartree–Fock (UHF) concrete mean-field method.
#pragma once

#include <string>

#include "atomiverse/cpu/scf/delta_scf.h"
#include "atomiverse/cpu/scf/eom_sf_ccsd.h"
#include "atomiverse/cpu/scf/mean_field_method.h"
#include "atomiverse/cpu/scf/scf_guess.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Unrestricted Hartree–Fock — the same mean-field, single-determinant
///        method as RHF but with independent spatial orbitals for the spin-up
///        (α) and spin-down (β) electrons, so open-shell / spin-polarized
///        systems (unequal α and β counts) are representable.
class UHF : public MeanFieldMethod<UHF> {
  public:
    int n_alpha = 0; ///< Alpha (spin-up) electrons (must be set before run).
    int n_beta = 0;  ///< Beta (spin-down) electrons (must be set before run).

    /// @brief Construct with the UHF defaults (unrestricted HF; Hcore initial
    ///        guess; 0.7 density damping; DIIS capped at 6 vectors and started
    ///        at iteration 5; energy tolerance 1e-8 Ha). Set `n_alpha` /
    ///        `n_beta` before `run`.
    UHF();
    /// @brief Construct with an explicit α / β electron count.
    /// @param n_alpha Spin-up electrons.
    /// @param n_beta  Spin-down electrons.
    explicit UHF(int n_alpha, int n_beta);
    std::string name() const { return "uhf"; }

    /// @brief EOM-SF-CCSD spin-flip excitation energies on a UHF → UCCSD
    ///        reference (drives SCF + UCCSD + SF-EOM).
    ///
    /// High-spin UHF (set `n_alpha` / `n_beta`), unrestricted spin-orbital
    /// CCSD, then equation-of-motion spin-flip singles+doubles. Matches
    /// PySCF `cc.UCCSD(mf).eomsf_ccsd(nroots=n_states)` on the same
    /// geometry / basis (cart). Product surface: ω in eV.
    ///
    /// @param frame    Geometry [Å] + Z; energy meta written from UHF.
    /// @param n_states Number of lowest SF roots (dimensionless, ≥ 1).
    /// @param guess    Optional density guess for the reference UHF (ROHF
    ///                 warm-start for multi-minima cases such as triplet O₂).
    /// @return SF excitation energies [eV] plus UCCSD energetics [Ha].
    /// @throws std::invalid_argument when SCF does not converge, n_states
    ///         is invalid, or the SF manifold is empty (no fallback).
    EOMSFCCSDResult eom_sf_ccsd(Frame& frame, int n_states = 1,
                                const SCFInitialGuess* guess = nullptr);

    /// @brief ΔSCF vertical excitation with MOM/IMOM occupation locking.
    ///
    /// Runs a ground-state UHF on `frame`, promotes one α electron from MO
    /// index `from` to virtual `to` via
    /// `SlaterDeterminant::with_occupation`, seeds the non-Aufbau density,
    /// and re-SCFs with `mom = style` (default `IMOM` = PySCF
    /// `scf.addons.mom_occ_`). Frame energy meta is left at the **excited**
    /// state; ground/excited Hartree energies and ΔE [eV] are on the return.
    ///
    /// @param frame Geometry [Å]; mutated only for energy meta on the
    ///              excited-state run.
    /// @param from  Occupied α MO index to depopulate (0-based; f=1).
    /// @param to    Virtual α MO index to populate (0-based; f=0).
    /// @param style `MomStyle::IMOM` (default) or `MomStyle::MOM`.
    /// @param guess Optional density guess for the ground-state UHF only.
    /// @return `DeltaScfResult` with e_ground / e_excited [Ha], delta_e_ev.
    /// @throws std::invalid_argument on bad indices or style==Off.
    /// @throws std::runtime_error when either SCF does not converge.
    DeltaScfResult delta_scf(Frame& frame, int from, int to,
                             MomStyle style = MomStyle::IMOM,
                             const SCFInitialGuess* guess = nullptr);

  protected:
    /// @brief Occupation is taken directly from the `n_alpha` / `n_beta`
    ///        members; the electron-count argument is unused for the open-shell
    ///        case.
    void resolve_occupation(int n_electrons) override;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
