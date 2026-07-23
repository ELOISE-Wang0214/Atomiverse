/// @file cpu/scf/adc.h
/// @brief Strict ADC(2) excitation energies (public API, result-only header).
///
/// User path: `RHF::adc2(n_states)` on a converged closed-shell RHF reference.
/// Intermediate-state representation of the PP propagator through second order
/// (Schirmer, Phys. Rev. A **26**, 2395 (1982); Trofimov & Schirmer, J. Phys. B
/// **28**, 2299 (1995)). Strict ADC(2) (not ADC(2)-x): second-order singles
/// block, first-order singles–doubles coupling, zeroth-order doubles diagonal.
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

} // namespace scf
} // namespace cpu
} // namespace ATV
