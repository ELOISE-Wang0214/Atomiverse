/// @file cpu/scf/uks.h
/// @brief Unrestricted Kohn–Sham (UKS) concrete mean-field method.
#pragma once

#include <string>

#include "atomiverse/cpu/scf/mean_field_method.h"
#include "atomiverse/cpu/scf/tdhf.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Unrestricted Kohn–Sham — Kohn–Sham density-functional theory with
///        independent spatial orbitals for the spin-up (α) and spin-down (β)
///        electrons, for open-shell / spin-polarized systems.
class UKS : public MeanFieldMethod<UKS> {
  public:
    int n_alpha = 0; ///< Alpha (spin-up) electrons (set before run).
    int n_beta = 0;  ///< Beta (spin-down) electrons (set before run).
    XCType xc =
        XCType::LSDA; ///< Exchange–correlation functional (default LSDA).

    /// @brief Construct with the UKS defaults (unrestricted KS; LSDA functional;
    ///        Hcore initial guess; 0.5 density damping; DIIS capped at 6
    ///        vectors; energy tolerance 1e-8 Ha). Set `n_alpha` / `n_beta`
    ///        before `run`.
    UKS();
    /// @brief Construct UKS with explicit α / β counts and functional.
    /// @param n_alpha Spin-up electrons.
    /// @param n_beta  Spin-down electrons.
    /// @param xc      Exchange–correlation functional (default LSDA).
    UKS(int n_alpha, int n_beta, XCType xc = XCType::LSDA);
    std::string name() const { return "uks"; }

    /// @brief Spin-conserving Tamm–Dancoff TDDFT on the most recent `run()`.
    ///
    /// Open-shell Casida **A** on the unrestricted KS reference (Hirata &
    /// Head-Gordon TDA; Dreuw & Head-Gordon, Chem. Rev. 105, 4009 (2005)):
    ///
    ///   \(A_{ia\sigma,jb\sigma'}=\delta_{\sigma\sigma'}\delta_{ij}\delta_{ab}
    ///   (\varepsilon_a^\sigma-\varepsilon_i^\sigma)
    ///   +(ia\sigma|jb\sigma')
    ///   -\alpha_{\mathrm{hf}}\delta_{\sigma\sigma'}(ij\sigma|ab\sigma')
    ///   +(ia\sigma|f^{\mathrm{xc}}_{\sigma\sigma'}|jb\sigma')\)
    ///
    /// For pure LSDA (`α_hf = 0`) the exchange subtraction vanishes and the
    /// adiabatic \(f^{\mathrm{xc}}\) is the spin-resolved second derivative of
    /// Slater+VWN5 on the KS grid (same path as PySCF `tdscf.TDA` on
    /// `dft.UKS`). There is no closed-shell singlet/triplet channel split —
    /// roots are Ms-conserving excitations of the open-shell reference.
    ///
    /// Dense path only (`n_ov_α + n_ov_β ≤ 512`). Dense `eigh` of **A**.
    /// @param n_states Number of lowest roots (1 ≤ n_states ≤ n_ov).
    /// @return Excitation energies [eV], length-gauge transition dipoles and
    ///         oscillator strengths; `tda = true`; `alda = true` for LSDA.
    /// @throws std::runtime_error when called before `run()`, on a restricted
    ///         determinant, or when n_ov exceeds 512 (no fallback).
    /// @throws std::invalid_argument when n_states is outside [1, n_ov].
    ExcitedStatesResult tda(int n_states) const;

  protected:
    /// @brief Occupation is taken directly from the `n_alpha` / `n_beta`
    ///        members; also latches the selected `xc` functional. The
    ///        electron-count argument is unused for the open-shell case.
    void resolve_occupation(int n_electrons) override;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
