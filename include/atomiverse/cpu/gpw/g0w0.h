/// @file cpu/gpw/g0w0.h
/// @brief Γ-point GPW spectral G0W0 quasiparticle result (public value type).
///
/// User path: `GpwDft::g0w0(atoms)` after a converged closed-shell
/// `run()` — the GPW analogue of molecular `RHF::g0w0()`
/// (`cpu/scf/g0w0.h`). Screening is direct RPA (dRPA) on the same
/// RI-GPW / Poisson `(ia|jb)` path as `GpwDft::direct_rpa`; the
/// correlation self-energy is the spectral (frequency-exact) sum over
/// RPA poles. Exact-exchange diagonals use truncated (MT) Coulomb so
/// the Γ-only foothold is G=0-safe without a periodic-HF product path.
///
/// Quasiparticle equation (Newton, orbit-by-orbit):
/// \f[
///   \omega = \varepsilon_p + \mathrm{Re}\,\Sigma_{c,p}(\omega)
///            + \Sigma_{x,p} - v_{\mathrm{xc},p}.
/// \f]
///
/// References: Hybertsen & Louie, Phys. Rev. B **34**, 5390 (1986);
/// Furche, J. Chem. Phys. **114**, 5982 (2001); molecular twin
/// `RHF::g0w0` / PySCF `gw.gw_exact.GWExact` form.
#pragma once

#include <vector>

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Diagonal spectral G0W0 quasiparticle energies of a closed-shell
///        Γ-point GPW RKS reference.
struct G0W0Result {
    /// @brief Quasiparticle orbital energies \(\varepsilon_p^{\mathrm{QP}}\)
    ///        [Hartree], same order as the KS MO spectrum (not re-sorted).
    std::vector<double> qp_energy_ha;

    /// @brief Mean-field (KS) orbital energies \(\varepsilon_p\) [Hartree].
    std::vector<double> mf_energy_ha;

    /// @brief Exact-exchange diagonals \(\Sigma_{x,p}\) [Hartree].
    std::vector<double> sigma_x_ha;

    /// @brief Local XC diagonals \(v_{\mathrm{xc},p}\) [Hartree].
    std::vector<double> vxc_ha;

    /// @brief HOMO→LUMO quasiparticle gap
    ///        \((\varepsilon_{\mathrm{LUMO}}^{\mathrm{QP}}
    ///        -\varepsilon_{\mathrm{HOMO}}^{\mathrm{QP}})\) [eV].
    double gap_ev = 0.0;

    /// @brief KS HOMO→LUMO gap [eV] (diagnostic).
    double mf_gap_ev = 0.0;

    /// @brief Vertical ionization potential \(-\varepsilon_{\mathrm{HOMO}}^{\mathrm{QP}}\)
    ///        [eV].
    double ip_ev = 0.0;

    /// @brief Vertical electron affinity \(\varepsilon_{\mathrm{LUMO}}^{\mathrm{QP}}\)
    ///        [eV].
    double ea_ev = 0.0;

    /// @brief Number of occupied spatial MOs (closed-shell).
    int n_occ = 0;

    /// @brief Virtual spatial MOs in the dRPA screening space (after cap).
    int n_virt = 0;

    /// @brief true when `(ia|jb)` used the RI-GPW product path.
    bool ri = false;

    /// @brief true when every orbital's Newton QP solve converged.
    bool converged = false;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
