/// @file cpu/pw/internal/kpoint.h
/// @brief Monkhorst-Pack k-point mesh generator.
#pragma once

#include "atomiverse/core/real.h"
#include "atomiverse/core/vec3.h"
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

struct KPointGrid {
    std::vector<Vec3> kpoints; // fractional coordinates
    std::vector<Real> weights;
};

KPointGrid generate_mp_mesh(int nk1, int nk2, int nk3, Real shift1 = 0,
                            Real shift2 = 0, Real shift3 = 0);

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
