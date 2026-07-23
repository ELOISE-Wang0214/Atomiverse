/// @file cpu/scf/adc.h
/// @brief Strict ADC(2) EE and IP results (public API, result-only header).
///
/// User paths on a converged closed-shell RHF reference:
/// - `RHF::adc2(n_states)` — electronic excitation (EE) manifold
/// - `RHF::ip_adc2(n_states)` — ionization potential (IP) manifold
///
/// Intermediate-state representation of the PP / 1p propagators through
/// second order (Schirmer, Phys. Rev. A **26**, 2395 (1982); Trofimov &
/// Schirmer, J. Phys. B **28**, 2299 (1995)). Strict ADC(2) (not ADC(2)-x):
/// second-order 1h/1p–1h block, first-order singles–doubles coupling,
/// zeroth-order 2h1p doubles diagonal.
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

/// @brief Strict ADC(2) ionization potentials on a closed-shell RHF
///        reference (PySCF `adc.ADC` / `method='adc(2)'` /
///        `method_type='ip'` convention).
struct IPADC2Result {
    /// @brief ADC(2) ionization potentials IP [eV], ascending.
    std::vector<double> energy_ev;

    /// @brief MP2 ground-state correlation energy of the same reference
    ///        [Hartree] (PySCF `e_corr` / doubles-only MP2).
    double e_corr_mp2_ha = 0.0;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
