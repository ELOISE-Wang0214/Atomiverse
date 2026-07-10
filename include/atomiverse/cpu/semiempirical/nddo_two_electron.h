/// @file cpu/semiempirical/nddo_two_electron.h
/// @brief NDDO reference-shaped two-center two-electron integral blocks.
#pragma once

#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/nddo_system.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

struct NDDOTwoCenterIntegralBlock {
    int atom_i = -1;
    int atom_j = -1;
    int n_pairs_i = 0;
    int n_pairs_j = 0;
    std::vector<double> values;

    double at(int pair_i, int pair_j) const;
};

/// @brief Port of AM1 reference `reppd/rotatd/w2mat` for AM1 s/p atoms.
NDDOTwoCenterIntegralBlock build_two_center_integral_block(
    const NDDOSystem& system, const cpu::Frame& frame, int atom_i, int atom_j);

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
