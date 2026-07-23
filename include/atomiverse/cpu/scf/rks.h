/// @file cpu/scf/rks.h
/// @brief Restricted Kohn–Sham (RKS) concrete mean-field method.
#pragma once

#include <string>

#include "atomiverse/cpu/scf/mean_field_method.h"
#include "atomiverse/cpu/scf/tdhf.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Restricted Kohn–Sham — Kohn–Sham density-functional theory (a
///        non-interacting orbital reference plus an exchange–correlation
///        functional that supplies the many-body effects) in its closed-shell
///        form, with spin-up (α) and spin-down (β) electrons sharing one set of
///        spatial orbitals. Defaults to the PBE functional on the MatchPySCF
///        integration grid.
class RKS : public MeanFieldMethod<RKS> {
  public:
    XCType xc = XCType::PBE; ///< Exchange–correlation functional (default PBE).

    /// @brief Construct with the RKS defaults (restricted KS; PBE functional;
    ///        MINAO initial guess; no damping; MatchPySCF grid; energy
    ///        tolerance 1e-9 Ha).
    RKS();
    /// @brief Construct RKS with an explicit exchange–correlation functional.
    /// @param xc Functional to use (e.g. LSDA, PBE, B3LYP, PBE0).
    explicit RKS(XCType xc);
    std::string name() const { return "rks"; }

    /// @brief Vertical TDA excitation energies of the most recent `run()`.
    ///
    /// Casida TDA on the converged KS orbitals:
    /// - **LSDA:** full MO-product ALDA kernel
    ///   \((ia|f^{\mathrm{xc}}|jb)=\int\phi_i\phi_a\,f_{\mathrm{ALDA}}(n)\,
    ///   \phi_j\phi_b\,d\mathbf{r}\) (Dirac-exchange second derivative;
    ///   Casida 1995) — `ExcitedStatesResult::alda = true`.
    /// - **Global hybrid (B3LYP / PBE0):** exact-exchange fraction
    ///   \(\alpha_{\mathrm{hf}}\) on the Casida **A** matrix
    ///   (\(^{1}A=\delta\delta\Delta\varepsilon+2(ia|jb)-\alpha_{\mathrm{hf}}(ij|ab)\))
    ///   plus the semi-local \(f^{\mathrm{xc}}\) kernel from a finite
    ///   difference of \(V_{\mathrm{xc}}\) on the KS grid (Dreuw &
    ///   Head-Gordon, Chem. Rev. 105, 4009 (2005)) —
    ///   `alda = false`.
    /// - **Pure GGA (PBE):** HF-like CIS@RKS fallback (`alda = false`).
    /// @param n_states Number of lowest excited states.
    /// @param channel  Singlet or triplet channel.
    /// @return Excitation energies [eV], transition properties;
    ///         `alda == true` for the LSDA full-grid path only.
    ExcitedStatesResult
    tda(int n_states,
        ExcitationChannel channel = ExcitationChannel::singlet) const;

  protected:
    /// @brief Closed-shell occupation (n_alpha = n_beta = n_electrons / 2);
    ///        also latches the selected `xc` functional into the method state.
    void resolve_occupation(int n_electrons) override;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
