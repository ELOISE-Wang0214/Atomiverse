/// @file cpu/scf/scf_enums.h
/// @brief The light SCF selector enums — zero-dependency public header so
///        internal headers can take the enums without dragging in result
///        types (rehomed verbatim from scf_result.h).
#pragma once

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Mean-field flavor: Hartree-Fock or Kohn-Sham (dimensionless tag).
enum class MeanFieldType { HF, KS };

/// @brief Spin treatment: restricted (shared spatial orbitals) or
///        unrestricted (independent alpha/beta channels).
enum class SpinType { Restricted, Unrestricted };

/// @brief SCF initial-guess flavor: bare core Hamiltonian or MINAO
///        superposition of atomic densities.
enum class InitGuess { Hcore, MINAO };

} // namespace scf
} // namespace cpu
} // namespace ATV
