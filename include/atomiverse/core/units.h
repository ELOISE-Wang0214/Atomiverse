/// @file core/units.h
/// @brief Central unit / physical-constant registry — the single source
///        of truth for every conversion factor used anywhere in ATV.
///
/// Organization mirrors `pint`'s mental model: each physical dimension
/// owns a small namespace, inside which live (a) named conversion
/// factors and (b) `X_from_Y(value)` helpers. New call sites prefer
/// the helper to an ad-hoc multiplication so intent is visible at the
/// call site:
///
/// @code
///   using namespace ATV::units;
///   double r_bohr   = length::bohr_from_angstrom(1.0);
///   double f_metal  = force::ev_per_angstrom_from_hartree_per_bohr(f_au);
///   double kT       = energy::ev_from_kelvin(298.15);
/// @endcode
///
/// ATV's two internal unit systems:
///   * MD / thermo / minimizer — Å, eV, fs, amu, eV/Å, K  (LAMMPS "metal")
///   * SCF core                — Bohr, Hartree, atomic-unit forces
///
/// Conversion between the two systems happens only at the Frame ↔
/// IntegralPack boundary and at the SCF-force → MD-force handoff.
///
/// Every round-trip to PySCF is pinned to PySCF's `pyscf/data/nist.py`
/// CODATA-2006 values (BOHR = 0.52917721092 Å,
/// HARTREE2EV = 27.211386245988) so that integrals agree bit-for-bit.
///
/// Back-compat: historical `ATV::k*` flat aliases remain at the bottom
/// of this file for call sites that haven't been migrated.
#pragma once

namespace ATV {
namespace units {

// ── Math ─────────────────────────────────────────────────────────────
inline constexpr double pi = 3.14159265358979323846;

// ── Base reference values (CODATA 2006, pinned to PySCF) ────────────
inline constexpr double bohr_per_angstrom_factor    = 1.8897261245650618; // Bohr = Å × this
inline constexpr double angstrom_per_bohr_factor    = 0.52917721092;      // Å    = Bohr × this
inline constexpr double ev_per_hartree_factor       = 27.211386245988;
inline constexpr double hartree_per_ev_factor       = 1.0 / 27.211386245988;
inline constexpr double boltzmann_ev_per_kelvin     = 8.617333262e-5;     // CODATA 2018

// ── Length ───────────────────────────────────────────────────────────
namespace length {
    /// Bohr ← Angstrom.
    [[nodiscard]] inline constexpr double bohr_from_angstrom(double x) {
        return x * bohr_per_angstrom_factor;
    }
    /// Angstrom ← Bohr.
    [[nodiscard]] inline constexpr double angstrom_from_bohr(double x) {
        return x * angstrom_per_bohr_factor;
    }
}

// ── Energy ───────────────────────────────────────────────────────────
namespace energy {
    /// eV ← Hartree.
    [[nodiscard]] inline constexpr double ev_from_hartree(double x) {
        return x * ev_per_hartree_factor;
    }
    /// Hartree ← eV.
    [[nodiscard]] inline constexpr double hartree_from_ev(double x) {
        return x * hartree_per_ev_factor;
    }
    /// eV ← kelvin (thermal energy k_B·T).
    [[nodiscard]] inline constexpr double ev_from_kelvin(double T_kelvin) {
        return boltzmann_ev_per_kelvin * T_kelvin;
    }
}

// ── Force ────────────────────────────────────────────────────────────
// The SCF gradient returns forces in Hartree/Bohr; the MD driver
// expects eV/Å. This sub-namespace owns that boundary conversion.
namespace force {
    /// 1 Hartree/Bohr in eV/Å.
    inline constexpr double ev_per_angstrom_per_hartree_per_bohr =
        ev_per_hartree_factor / angstrom_per_bohr_factor;

    /// eV/Å ← Hartree/Bohr.
    [[nodiscard]] inline constexpr double ev_per_angstrom_from_hartree_per_bohr(double f_au) {
        return f_au * ev_per_angstrom_per_hartree_per_bohr;
    }
    /// Hartree/Bohr ← eV/Å.
    [[nodiscard]] inline constexpr double hartree_per_bohr_from_ev_per_angstrom(double f_metal) {
        return f_metal / ev_per_angstrom_per_hartree_per_bohr;
    }
}

// ── Angle ────────────────────────────────────────────────────────────
namespace angle {
    inline constexpr double rad_per_degree = pi / 180.0;
    inline constexpr double degree_per_rad = 180.0 / pi;

    [[nodiscard]] inline constexpr double rad_from_degree(double x) {
        return x * rad_per_degree;
    }
    [[nodiscard]] inline constexpr double degree_from_rad(double x) {
        return x * degree_per_rad;
    }
}

// ── MD "metal" internal unit coupling ────────────────────────────────
//
// Every atomistic step lives in Å / eV / fs / amu / K. Two conversion
// factors make `F/m → a` and `½ m v² → KE` work in that system without
// any reference to SI:
//
//   force_to_accel             [F/m × this → a, units Å/fs²]
//   mass_velocity2_to_energy   [m × v² × this → E, units eV]
//
// These are the same numbers LAMMPS uses for `metal` units, derived by
// requiring Newton's second law and kinetic-energy formula to be
// dimensionally consistent under the chosen base units.
inline constexpr double force_to_accel             = 1.0 / 103.6427;
inline constexpr double mass_velocity2_to_energy   = 103.6427;

} // namespace units

// ─────────────────────────────────────────────────────────────────────
// Back-compat: flat `ATV::k*` names for call sites that predate the
// `units::` namespace. Prefer the namespaced helpers in new code.
// ─────────────────────────────────────────────────────────────────────

inline constexpr double kPI              = units::pi;
inline constexpr double kDegToRad        = units::angle::rad_per_degree;
inline constexpr double kAngstromToBohr  = units::bohr_per_angstrom_factor;
inline constexpr double kBohr2Ang        = units::angstrom_per_bohr_factor;
inline constexpr double kHartree2eV      = units::ev_per_hartree_factor;
inline constexpr double kForceAU2Metal   = units::force::ev_per_angstrom_per_hartree_per_bohr;
inline constexpr double kBoltzmann       = units::boltzmann_ev_per_kelvin;
inline constexpr double kFtm2v           = units::force_to_accel;
inline constexpr double kMvv2e           = units::mass_velocity2_to_energy;

} // namespace ATV
