/// @file cpu/scf/scf_enums.h
/// @brief The light SCF selector enums — zero-dependency public header so
///        internal headers can take the enums without dragging in result
///        types (rehomed from the since-deleted scf_result.h, 02-rehome).
#pragma once

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Mean-field flavor: Hartree-Fock or Kohn-Sham (dimensionless tag).
enum class MeanFieldType { HF, KS };

/// @brief Spin treatment: restricted (shared spatial orbitals),
///        unrestricted (independent alpha/beta channels),
///        restricted-open (single MO set with split alpha/beta occupation —
///        ROHF; Roothaan effective-Fock coupling), or generalized (GHF:
///        two-component real spinor orbitals φ_i = (φ_i^α; φ_i^β), one
///        2n×2n density with a spin-off-diagonal P^αβ block — the
///        noncollinear-spin mean field; Fukutome, Int. J. Quantum Chem. 20,
///        955 (1981), DOI 10.1002/qua.560200502).
///
/// Paths that dispatch on this enum and do not implement the
/// restricted-open or generalized case (analytic gradients, MP2/FCI packs,
/// SOSCF, level shift) throw explicitly for `RestrictedOpen` /
/// `Generalized` — no silent Restricted/Unrestricted branch is ever taken.
enum class SpinType { Restricted, Unrestricted, RestrictedOpen, Generalized };

/// @brief SCF initial-guess flavor.
///
/// - `Hcore`:   diagonalize the bare core Hamiltonian (worst quality).
/// - `MINAO`:   project MINAO atomic fragment densities (PySCF default).
/// - `Hueckel`: diagonalize the extended-Hückel Fock built from the static
///              VSIP table (`internal::HueckelGuess`); Lehtola,
///              JCTC 15, 1593 (2019), DOI 10.1021/acs.jctc.8b01089.
/// - `SAP`:     diagonalize F = T + V_SAP, the superposition of tabulated
///              atomic potentials on the DFT grid (`internal::SapGuess`);
///              same reference.
enum class InitGuess { Hcore, MINAO, Hueckel, SAP };

/// @brief SCF Fock-mixer family for the mix stage.
///
/// - `CDIIS`:       Pulay commutator-DIIS only (the shipping default).
/// - `EDIIS_CDIIS`: energy-DIIS far from convergence, permanently switching
///                  to CDIIS once max|DIIS error| < `mixer_switch_tol`
///                  (Kudin-Scuseria-Cancès, DOI 10.1063/1.1470195).
/// - `ADIIS_CDIIS`: augmented-Roothaan-Hall DIIS far from convergence with
///                  the same permanent CDIIS switch (Hu-Yang,
///                  DOI 10.1063/1.3304922).
///
/// EDIIS/ADIIS are far-field mixers only — run to the end they stall well
/// above tight convergence, hence the escalation-to-CDIIS design.
enum class MixerKind { CDIIS, EDIIS_CDIIS, ADIIS_CDIIS };

} // namespace scf
} // namespace cpu
} // namespace ATV
