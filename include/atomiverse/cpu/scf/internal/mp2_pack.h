/// @file cpu/scf/internal/mp2_pack.h
/// @brief Private pack-based MP2 entry (unit tests / method implementation).
#pragma once

#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/slater_determinant.h" // SlaterDeterminant

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Closed-shell canonical MP2 energy from a converged reference.
/// @param ints     Integral pack with the ERI tensor materialized.
/// @param det      Converged closed-shell reference determinant (MO
///                 coefficients, orbital energies, electron counts).
/// @param e_scf_ha Converged SCF total energy [Hartree] (reported as the
///                 reference part of `MP2Result::e_tot`).
/// @param n_frozen Number of lowest occupied spatial MOs excluded from the
///                 correlation sum (frozen core). 0 = all occupied correlated.
///                 Matches PySCF `mp.MP2(mf, frozen=n_frozen)`.
/// @return MP2 correlation and total energies [Hartree].
MP2Result mp2_energy(const IntegralPack& ints, const SlaterDeterminant& det,
                     double e_scf_ha, int n_frozen = 0);

/// @brief Unrestricted canonical MP2 (UMP2) on a UHF/UKS-style reference.
///
/// Spin-orbital second-order energy decomposed into same-spin αα/ββ
/// (antisymmetrized integrals, factor ¼) and opposite-spin αβ channels:
///
/// \f[
///   E^{(2)}_{\sigma\sigma}
///     = \tfrac14\sum_{ijab}
///       \frac{|(ia|jb)-(ib|ja)|^2}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b},
/// \quad
///   E^{(2)}_{\alpha\beta}
///     = \sum_{i_\alpha j_\beta a_\alpha b_\beta}
///       \frac{(i_\alpha a_\alpha|j_\beta b_\beta)^2}
///            {\varepsilon_{i\alpha}+\varepsilon_{j\beta}
///             -\varepsilon_{a\alpha}-\varepsilon_{b\beta}}.
/// \f]
///
/// Matches raw (unprojected) PySCF `mp.UMP2`. No frozen core.
/// @param ints     Integral pack with the ERI tensor materialized.
/// @param det      Converged unrestricted determinant (independent α/β MOs).
/// @param e_scf_ha Converged UHF total energy [Hartree].
/// @return UMP2 correlation and total energies [Hartree].
MP2Result ump2_energy(const IntegralPack& ints, const SlaterDeterminant& det,
                      double e_scf_ha);

/// @brief Closed-shell RI-MP2 correlation energy via Coulomb-metric density
///        fitting (same `DensityFitting` workspace as SCF RI-J/JK). Does not
///        require a materialized 4-center ERI tensor.
/// @param ints     Integral pack (shells + atoms in Bohr).
/// @param det      Converged closed-shell RHF determinant.
/// @param e_scf_ha SCF energy [Hartree].
/// @param aux      Auxiliary basis (def2 universal JKFIT).
/// @param n_frozen Lowest occupied spatial MOs frozen (same convention as
///                 `mp2_energy`).
/// @return MP2Result with RI (ia|jb)≈Σ_P B_ia^P B_jb^P [Hartree].
MP2Result ri_mp2_energy(const IntegralPack& ints, const SlaterDeterminant& det,
                        double e_scf_ha, AuxBasis aux, int n_frozen = 0);

} // namespace scf
} // namespace cpu
} // namespace ATV
