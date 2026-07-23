/// @file cpu/scf/adc.h
/// @brief Strict ADC(2) excitation / attachment energies (public API).
///
/// User paths on a converged closed-shell RHF reference:
/// - `RHF::adc2(n_states)` — electronic excitation (EE) manifold
/// - `RHF::ea_adc2(n_states)` — electron-attachment (EA) manifold
///
/// Intermediate-state representation of the PP / EA propagator through second
/// order (Schirmer, Phys. Rev. A **26**, 2395 (1982); Trofimov & Schirmer,
/// J. Phys. B **28**, 2299 (1995)). Strict ADC(2) (not ADC(2)-x): second-order
/// singles block, first-order singles–doubles coupling, zeroth-order doubles
/// diagonal.
#pragma once

#include <vector>

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Strict ADC(2) electronic excitation energies on a closed-shell RHF
///        reference (singlet manifold, PySCF `adc.ADC` / `method='adc(2)'` /
///        `method_type='ee'` convention).
struct ADC2Result {
    /// @brief ADC(2) excitation energies ω [eV], ascending (Davidson Ritz).
    std::vector<double> energy_ev;

    /// @brief MP2 ground-state correlation energy of the same reference
    ///        [Hartree] (PySCF `e_corr` / doubles-only MP2).
    double e_corr_mp2_ha = 0.0;
};

/// @brief Strict EA-ADC(2) electron-attachment energies on a closed-shell RHF
///        reference (PySCF `adc.ADC` / `method='adc(2)'` / `method_type='ea'`).
///
/// Singles space is the virtual (1p) manifold; doubles are 2p–1h. Reported
/// energies are the EA eigenvalues ω [eV] of the EA-ADC matrix (positive for
/// endothermic attachment, e.g. H₂/STO-3G LUMO).
struct EAADC2Result {
    /// @brief EA-ADC(2) attachment energies ω [eV], ascending (Davidson Ritz).
    std::vector<double> energy_ev;

    /// @brief MP2 ground-state correlation energy of the same N-electron
    ///        reference [Hartree] (PySCF `e_corr` / doubles-only MP2).
    double e_corr_mp2_ha = 0.0;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
