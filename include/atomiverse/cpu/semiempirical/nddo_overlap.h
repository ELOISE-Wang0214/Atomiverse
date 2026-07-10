/// @file cpu/semiempirical/nddo_overlap.h
/// @brief NDDO reference-shaped diatomic overlap and h1elec helpers.
#pragma once

#include <array>

namespace ATV {
namespace cpu {
namespace semiempirical {

struct NDDOLocalMatrix9 {
    std::array<double, 81> values{};

    double at(int row, int col) const;
    double& at(int row, int col);
};

/// @brief Port of AM1 reference `diat`: local AO overlap for a two-atom vector.
NDDOLocalMatrix9 diatomic_overlap(int z_i, int z_j,
                                  const std::array<double, 3>& x_j_minus_i);

/// @brief Port of AM1 reference `h1elec`: two-center one-electron matrix block.
NDDOLocalMatrix9 one_electron_two_center(int z_i, int z_j,
                                         const std::array<double, 3>& x_i,
                                         const std::array<double, 3>& x_j);

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
