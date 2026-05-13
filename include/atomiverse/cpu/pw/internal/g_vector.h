/// @file cpu/pw/internal/g_vector.h
/// @brief G-vector enumeration for plane-wave basis.
#pragma once

#include <vector>
#include "atomiverse/core/real.h"
#include "atomiverse/core/vec3.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/box.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

struct GVectorList {
    std::vector<Vec3> g;        // G-vectors in Cartesian coords [1/Bohr]
    std::vector<Real> ekin;     // kinetic energy 0.5*|k+G|^2 [Hartree]
    int nx, ny, nz;             // FFT grid dimensions

    int count() const { return static_cast<int>(g.size()); }
};

GVectorList enumerate_g_vectors(const Box& cell_bohr, Real ecut, const Vec3& k);

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
