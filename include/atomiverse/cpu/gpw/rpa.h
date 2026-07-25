/// @file cpu/gpw/rpa.h
/// @brief Γ-point GPW direct RPA (dRPA) correlation energy result (public
///        value type).
///
/// User path: `GpwDft::direct_rpa(atoms)` after a converged closed-shell
/// `run()` — the GPW analogue of molecular `RHF::direct_rpa()`
/// (`cpu/scf/rpa.h`). Product default reuses the RI-GPW plane-wave
/// assembly of \((ia|jb)\) from `GpwDft::mp2` (`mp2_ri = true`); the
/// Poisson foothold remains available with `mp2_ri = false`.
///
/// Correlation energy via the Furche plasmon formula on the full
/// Coulomb-only (\(\alpha_{\mathrm{hf}}=0\)) Casida spectrum:
/// \f[
///   E_{\mathrm{c}}^{\mathrm{dRPA}}
///     = \tfrac12\Bigl(\sum_s \Omega_s - \mathrm{Tr}\,A\Bigr),
/// \f]
/// with singlet spin-adapted \(A,B\) and \(\Omega_s\) the positive
/// eigenvalues of \((A-B)^{1/2}(A+B)(A-B)^{1/2}\). Total energy is the
/// dRPA@KS sum \(E_{\mathrm{SCF}}+E_c\) (no EXX replacement — periodic
/// exact exchange remains blocked).
///
/// References: Furche, J. Chem. Phys. **114**, 5982 (2001); Ren et al.,
/// New J. Phys. **14**, 053020 (2012); Del Ben et al., J. Chem. Theory
/// Comput. **8**, 4177 (2012) for the GPW \((ia|jb)\) route.
#pragma once

#include <vector>

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Γ-point GPW direct RPA correlation energy (closed-shell RKS;
///        RI-GPW or Poisson foothold for the MO repulsion tensor).
struct RpaResult {
    /// @brief dRPA correlation energy \(E_c\) [Hartree] (typically < 0).
    double e_corr = 0.0;
    /// @brief Converged GPW SCF total energy [Hartree].
    double e_scf = 0.0;
    /// @brief e_scf + e_corr [Hartree].
    double e_total = 0.0;
    /// @brief \(\mathrm{Tr}\,A\) of the singlet dRPA \(A\) matrix [Hartree].
    double tr_A = 0.0;
    /// @brief Correlated occupied spatial orbitals.
    int n_occ = 0;
    /// @brief Virtual spatial orbitals in the spectrum (after any cap).
    int n_virt = 0;
    /// @brief true when (ia|jb) used the RI-GPW G-space product path.
    bool ri = false;
    /// @brief Positive dRPA excitation energies \(\Omega_s\) [Hartree],
    ///        ascending (full spectrum, length \(n_{\mathrm{occ}}\cdot
    ///        n_{\mathrm{virt}}\)).
    std::vector<double> omega_ha;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
