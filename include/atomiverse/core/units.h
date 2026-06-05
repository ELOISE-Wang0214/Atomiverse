/// @file core/units.h
/// @brief Central unit / physical-constant registry — the single source
///        of truth for every conversion factor used in ATV.
///
/// Each physical dimension is an enum with one value per supported unit.
/// A single `convert(x, from, to)` overload per dimension performs any
/// conversion in that dimension. Internally one factor per unit is
/// stored (value of `1 [unit]` measured in the dimension's base unit);
/// conversion reduces to `x * factor(from) / factor(to)` — O(N)
/// constants, O(1) call, O(N²) pair support "for free".
///
/// @code
///   using namespace ATV::units;
///   double r_bohr = convert(r_ang, Length::Angstrom, Length::Bohr);
///   double e_ev   = convert(e_ha,  Energy::Hartree,  Energy::eV);
///   double f_ev_A = convert(f_au,  Force::HartreePerBohr, Force::eVPerAngstrom);
/// @endcode
///
/// Base units (first enumerator in each dimension):
///   Length: Bohr      Energy: Hartree      Angle: Rad      Force: Hartree/Bohr
///
/// CODATA 2006 round-trip values (BOHR = 0.52917721092 Å,
/// HARTREE2EV = 27.211386245988) are pinned to PySCF's
/// `pyscf/data/nist.py` so integrals match the ground-truth path
/// bit-for-bit.
#pragma once

namespace ATV {
namespace units {

// ── Named factors (value of "1 X" in the dimension's base unit) ─────
//
// These are the single source of truth for every numeric conversion.
// Both the enum-dispatched convert() below and the legacy ATV::k*
// aliases at the bottom of this header read from these constants.

inline constexpr double kAngstromInBohr = 1.8897261245650618;
inline constexpr double kEvInHartree = 1.0 / 27.211386245988;
inline constexpr double kDegreeInRad = 3.14159265358979323846 / 180.0;
inline constexpr double kEvPerAngstromInHartreePerBohr =
    (1.0 / 27.211386245988)
    / 1.8897261245650618; // = (1 eV / 1 Å) in (Hartree / Bohr)

// ── Dimension enums ──────────────────────────────────────────────────

enum class Length { Bohr, Angstrom };
enum class Energy { Hartree, eV };
enum class Angle { Rad, Degree };
enum class Force { HartreePerBohr, eVPerAngstrom };

// ── Enum-dispatched factor lookup (internal) ────────────────────────

namespace detail {

[[nodiscard]] inline constexpr double factor(Length u) {
    switch (u) {
    case Length::Bohr:
        return 1.0;
    case Length::Angstrom:
        return kAngstromInBohr;
    }
    return 0.0; // unreachable
}

[[nodiscard]] inline constexpr double factor(Energy u) {
    switch (u) {
    case Energy::Hartree:
        return 1.0;
    case Energy::eV:
        return kEvInHartree;
    }
    return 0.0;
}

[[nodiscard]] inline constexpr double factor(Angle u) {
    switch (u) {
    case Angle::Rad:
        return 1.0;
    case Angle::Degree:
        return kDegreeInRad;
    }
    return 0.0;
}

[[nodiscard]] inline constexpr double factor(Force u) {
    switch (u) {
    case Force::HartreePerBohr:
        return 1.0;
    case Force::eVPerAngstrom:
        return kEvPerAngstromInHartreePerBohr;
    }
    return 0.0;
}

} // namespace detail

// ── Single public convert() per dimension ────────────────────────────

template <typename Dim>
[[nodiscard]] inline constexpr double convert(double x, Dim from, Dim to) {
    return x * detail::factor(from) / detail::factor(to);
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
//   acceleration [Å/fs²] = force_to_accel    * (force[eV/Å] / mass[amu])
//   KE [eV]              = 0.5 * mass_v2_to_energy * mass[amu] * v² [Å²/fs²]
inline constexpr double force_to_accel = 1.0 / 103.6427;
inline constexpr double mass_v2_to_energy = 103.6427;

// Electrostatic coupling in MD "metal" units (LAMMPS `qqr2e`):
//
//   E [eV] = coulomb * q_i q_j / r   with q in elementary charge, r in Å
//
// Equals 1/(4πε₀) · e²/Å expressed in eV (CODATA): 14.399645 eV·Å/e².
inline constexpr double coulomb = 14.399645;

} // namespace units

// ─────────────────────────────────────────────────────────────────────
// Back-compat: flat `ATV::k*` names for legacy call sites. Read
// directly from the named factors above — no `units::detail::`
// involvement.
// ─────────────────────────────────────────────────────────────────────

inline constexpr double kPI = units::pi;
inline constexpr double kDegToRad = units::kDegreeInRad;
inline constexpr double kAngstromToBohr = units::kAngstromInBohr;
inline constexpr double kBohr2Ang = 1.0 / units::kAngstromInBohr;
inline constexpr double kHartree2eV = 1.0 / units::kEvInHartree;
inline constexpr double kForceAU2Metal =
    1.0 / units::kEvPerAngstromInHartreePerBohr;
inline constexpr double kBoltzmann = units::kB_ev_per_k;
inline constexpr double kFtm2v = units::force_to_accel;
inline constexpr double kMvv2e = units::mass_v2_to_energy;

/// @brief Electrostatic constant 1/(4πε₀)·e²/Å in eV (LAMMPS `qqr2e`,
/// metal units): `E[eV] = kCoulomb · q_i q_j / r` for q[e], r[Å].
inline constexpr double kCoulomb = units::coulomb;

/// @brief Angstrom to Bohr (QE convention: 1 Bohr = 0.5291772083 A).
constexpr double kAngstrom2Bohr = 1.8897259886;

/// @brief Rydberg to Hartree (exact: 1 Ry = 1/2 Hartree).
constexpr double kRy2Hartree = 0.5;

/// @brief Rydberg to eV (QE convention).
constexpr double kRy2eV = 13.605693;

} // namespace ATV
