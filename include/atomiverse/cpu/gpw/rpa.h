/// @file cpu/gpw/rpa.h
/// @brief GPW direct RPA (dRPA) correlation energy result (public value type).
///
/// User path: `GpwDft::direct_rpa(atoms)` after a converged closed-shell
/// `run()` — the GPW analogue of molecular `RHF::direct_rpa()`
/// (`cpu/scf/rpa.h`). Supports **Γ** (real MOs; RI-GPW or Poisson) and
/// **multi-k** (complex Bloch + finite-\(q\) inter-\(k\) Casida; RI-GPW
/// product only).
///
/// Correlation energy via the Furche plasmon formula on the
/// Coulomb-only (\(\alpha_{\mathrm{hf}}=0\)) Casida spectrum:
/// \f[
///   E_{\mathrm{c}}^{\mathrm{dRPA}}
///     = \tfrac12\Bigl(\sum_s \Omega_s - \mathrm{Tr}\,A\Bigr),
/// \f]
/// multi-k product: \(E_c=\sum_q w_q\,\tfrac12\bigl(\sum_s
/// \Omega_{q,s}-\mathrm{Tr}\,A_q\bigr)\) with transitions \(i_k\to a_{k+q}\)
/// and Coulomb at transfer \(q\). Total energy is the dRPA@KS sum
/// \(E_{\mathrm{SCF}}+E_c\) (no EXX replacement — periodic exact exchange
/// remains blocked).
///
/// References: Furche, J. Chem. Phys. **114**, 5982 (2001); Ren et al.,
/// New J. Phys. **14**, 053020 (2012); Del Ben et al., J. Chem. Theory
/// Comput. **8**, 4177 (2012) for the GPW \((ia|jb)\) route.
#pragma once

#include <vector>

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief GPW direct RPA correlation energy (closed-shell RKS; Γ or multi-k).
struct RpaResult {
    /// @brief dRPA correlation energy \(E_c\) [Hartree] (typically < 0).
    double e_corr = 0.0;
    /// @brief Converged GPW SCF total energy [Hartree].
    double e_scf = 0.0;
    /// @brief e_scf + e_corr [Hartree].
    double e_total = 0.0;
    /// @brief \(\mathrm{Tr}\,A\) (Γ) or \(\sum_q w_q\,\mathrm{Tr}\,A_q\)
    ///        (multi-k finite-\(q\)) [Hartree].
    double tr_A = 0.0;
    /// @brief Correlated occupied spatial orbitals.
    int n_occ = 0;
    /// @brief Virtual spatial orbitals in the spectrum (after any cap).
    int n_virt = 0;
    /// @brief true when (ia|jb) used the RI-GPW G-space product path.
    bool ri = false;
    /// @brief Sampled k-points used in the RPA sum (0 = real Γ path).
    int n_kpts = 0;
    /// @brief Positive dRPA excitation energies \(\Omega_s\) [Hartree],
    ///        ascending. Γ: length \(n_{\mathrm{occ}}\cdot n_{\mathrm{virt}}\).
    ///        Multi-k: all finite-\(q\) block roots concatenated then sorted
    ///        (length \(N_k^2\cdot n_{\mathrm{occ}}\cdot n_{\mathrm{virt}}\)
    ///        with \(q\)-mesh = \(k\)-mesh).
    std::vector<double> omega_ha;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
