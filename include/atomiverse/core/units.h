/// @file core/units.h
/// @brief Central unit / physical-constant registry — the single source
///        of truth for every conversion factor used in ATV.
///
/// Each physical dimension is an enum with one value per supported unit.
/// One `convert(x, from, to)` overload per dimension performs the
/// actual conversion. Internally each unit stores a single factor "how
/// much of this unit equals one of the dimension's base unit"; any
/// `A → B` is then `x * (A_in_base / B_in_base)` — a single multiply
/// and divide. Adding a new unit is one line in the factor table; the
/// number of round-trip constants stays linear, not quadratic.
///
/// @code
///   using namespace ATV::units;
///   double r_bohr = convert(r_ang, Length::Angstrom, Length::Bohr);
///   double e_ev   = convert(e_ha,  Energy::Hartree,  Energy::eV);
///   double f_ev_A = convert(f_au,  Force::HartreePerBohr, Force::eVPerAngstrom);
/// @endcode
///
/// ATV's two internal unit systems:
///   * MD / thermo / minimizer — Å, eV, fs, amu, eV/Å, K  (LAMMPS "metal")
///   * SCF core                — Bohr, Hartree, atomic-unit forces
///
/// Conversion between the two happens only at the Frame ↔ IntegralPack
/// boundary and at the SCF-force → MD-force handoff.
///
/// CODATA 2006 round-trip values (BOHR = 0.52917721092 Å,
/// HARTREE2EV = 27.211386245988) are pinned to PySCF's
/// `pyscf/data/nist.py` so integrals match the ground-truth path
/// bit-for-bit.
#pragma once

namespace ATV {
namespace units {

// ── Dimension enums ──────────────────────────────────────────────────

enum class Length { Bohr, Angstrom };
enum class Energy { Hartree, eV };
enum class Angle  { Rad, Degree };
enum class Force  { HartreePerBohr, eVPerAngstrom };

// ── Factor-to-base tables (internal) ────────────────────────────────
//
// Base unit per dimension: the first enumerator (Bohr / Hartree / Rad
// / HartreePerBohr). `in_base(X)` returns how many base-units equal
// one X, so `convert` reduces to `x * in_base(from) / in_base(to)`.

namespace detail {

[[nodiscard]] inline constexpr double in_base(Length u) {
    switch (u) {
        case Length::Bohr:     return 1.0;
        case Length::Angstrom: return 1.8897261245650618;  // PySCF CODATA 2006
    }
    return 0.0;  // unreachable
}

[[nodiscard]] inline constexpr double in_base(Energy u) {
    switch (u) {
        case Energy::Hartree: return 1.0;
        case Energy::eV:      return 1.0 / 27.211386245988;  // Ha per eV
    }
    return 0.0;
}

[[nodiscard]] inline constexpr double in_base(Angle u) {
    switch (u) {
        case Angle::Rad:    return 1.0;
        case Angle::Degree: return 3.14159265358979323846 / 180.0;
    }
    return 0.0;
}

[[nodiscard]] inline constexpr double in_base(Force u) {
    // Base = Hartree/Bohr. 1 eV/Å in Hartree/Bohr:
    //    (1/27.211…) Hartree / (1/0.52917…) Bohr
    //  = 0.52917… / 27.211…
    switch (u) {
        case Force::HartreePerBohr:  return 1.0;
        case Force::eVPerAngstrom:   return 0.52917721092 / 27.211386245988;
    }
    return 0.0;
}

} // namespace detail

// ── Single convert per dimension ─────────────────────────────────────

template <typename Dim>
[[nodiscard]] inline constexpr double convert(double x, Dim from, Dim to) {
    return x * detail::in_base(from) / detail::in_base(to);
}

// ── Named physical constants (not unit conversions) ─────────────────

/// π.
inline constexpr double pi = 3.14159265358979323846;

/// Boltzmann constant in MD-internal units: eV per kelvin (CODATA 2018).
/// Thermal energy: `kT [eV] = kB_ev_per_k * T [K]`.
inline constexpr double kB_ev_per_k = 8.617333262e-5;

// ── MD "metal" internal coupling (LAMMPS convention) ────────────────
//
// Å / eV / fs / amu / K is not SI-consistent; these two factors make
// Newton's second law and kinetic-energy formula work without any SI
// intermediate:
//
//   acceleration [Å/fs²] = force_to_accel  * (force[eV/Å] / mass[amu])
//   KE [eV]              = 0.5 * mass_v2_to_energy * mass[amu] * v² [Å²/fs²]
inline constexpr double force_to_accel    = 1.0 / 103.6427;
inline constexpr double mass_v2_to_energy = 103.6427;

} // namespace units

// ─────────────────────────────────────────────────────────────────────
// Back-compat: flat `ATV::k*` names for call sites that predate the
// enum-based `units::convert()` API. New code uses the namespaced
// form above; these aliases keep ~95 legacy call sites compiling.
// ─────────────────────────────────────────────────────────────────────

inline constexpr double kPI             = units::pi;
inline constexpr double kDegToRad       = units::detail::in_base(units::Angle::Degree);
inline constexpr double kAngstromToBohr = units::detail::in_base(units::Length::Angstrom);
inline constexpr double kBohr2Ang       = 1.0 / kAngstromToBohr;
inline constexpr double kHartree2eV     = 1.0 / units::detail::in_base(units::Energy::eV);
inline constexpr double kForceAU2Metal  = 1.0 / units::detail::in_base(units::Force::eVPerAngstrom);
inline constexpr double kBoltzmann      = units::kB_ev_per_k;
inline constexpr double kFtm2v          = units::force_to_accel;
inline constexpr double kMvv2e          = units::mass_v2_to_energy;

} // namespace ATV
