/// @file cpu/semiempirical/nddo_system.h
/// @brief NDDO reference-shaped NDDO system indexing and population helpers.
#pragma once

#include <cstddef>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

struct NDDOAtom {
    int atomic_number = 0;
    int first_ao = 0; ///< 0-based analogue of NDDO reference nfirst(i)-1.
    int last_ao = -1; ///< 0-based analogue of NDDO reference nlast(i)-1.
    int n_orbitals = 0;
    double core_charge =
        0.0; ///< NDDO reference `tore`: valence electron/core charge.
};

struct NDDOSystem {
    std::vector<NDDOAtom> atoms;
    int n_orbitals = 0;
    int n_electrons = 0;
    int n_alpha = 0;
    int n_beta = 0;
    int charge = 0;
};

/// @brief Number of AM1 valence orbitals for atomic number `z`.
int am1_orbital_count(int z) noexcept;

/// @brief NDDO reference AM1 valence/core charge (`tore`) for atomic number `z`.
double am1_core_charge(int z) noexcept;

/// @brief Build a NDDO reference-shaped AO indexing snapshot from an Atomiverse frame.
NDDOSystem build_am1_system(const cpu::Frame& frame, int charge,
                            int multiplicity);

/// @brief Packed lower-triangle index for 0-based AO indices.
std::size_t packed_index(int i, int j) noexcept;

/// @brief NDDO reference `chrge`: sum density diagonal entries over each atom's AO range.
std::vector<double> atomic_populations_from_packed_density(
    const NDDOSystem& system, const std::vector<double>& packed_density);

/// @brief Net Mulliken-like atom charges: core charge minus electron population.
std::vector<double>
atomic_charges_from_packed_density(const NDDOSystem& system,
                                   const std::vector<double>& packed_density);

/// @brief NDDO reference initial diagonal density guess (`pdiag`) for AM1 RHF setup.
std::vector<double> initial_diagonal_density_guess(const NDDOSystem& system);

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
