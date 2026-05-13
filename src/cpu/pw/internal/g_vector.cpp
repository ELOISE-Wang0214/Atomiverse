/// @file cpu/pw/internal/g_vector.cpp
/// @brief G-vector enumeration implementation.

#include "atomiverse/cpu/pw/internal/g_vector.h"
#include <cmath>
#include <algorithm>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

GVectorList enumerate_g_vectors(const Box& cell_bohr, Real ecut, const Vec3& k) {
    auto recip = cell_bohr.recip_cell();
    Vec3 b0{recip(0,0), recip(0,1), recip(0,2)};
    Vec3 b1{recip(1,0), recip(1,1), recip(1,2)};
    Vec3 b2{recip(2,0), recip(2,1), recip(2,2)};

    Real b0n = b0.norm(), b1n = b1.norm(), b2n = b2.norm();
    Real gmax = std::sqrt(Real(2) * ecut);
    int hmax = static_cast<int>(gmax / b0n) + 1;
    int kmax = static_cast<int>(gmax / b1n) + 1;
    int lmax = static_cast<int>(gmax / b2n) + 1;

    struct GV { Vec3 v; Real e; };
    std::vector<GV> tmp;

    for (int h = -hmax; h <= hmax; ++h) {
        for (int kk = -kmax; kk <= kmax; ++kk) {
            for (int l = -lmax; l <= lmax; ++l) {
                Vec3 g = b0 * Real(h) + b1 * Real(kk) + b2 * Real(l);
                Vec3 kpg = {g[0] + k[0], g[1] + k[1], g[2] + k[2]};
                Real e = Real(0.5) * kpg.norm_squared();
                if (e <= ecut + Real(1e-12)) {
                    tmp.push_back({g, e});
                }
            }
        }
    }

    std::sort(tmp.begin(), tmp.end(),
        [](const GV& a, const GV& b) { return a.v.norm_squared() < b.v.norm_squared(); });

    GVectorList result;
    result.g.reserve(tmp.size());
    result.ekin.reserve(tmp.size());
    for (auto& t : tmp) { result.g.push_back(t.v); result.ekin.push_back(t.e); }

    // FFT grid: find max |h|,|k|,|l|
    int max_h = 0, max_k = 0, max_l = 0;
    Real ib0 = Real(1) / b0n, ib1 = Real(1) / b1n, ib2 = Real(1) / b2n;
    for (auto& gv : result.g) {
        Real hh = gv.dot(b0) * ib0 / (b0.norm() * ib0); // project onto b0
        int hi = static_cast<int>(std::abs(gv[0]) / b0n) + 1;
        int ki = static_cast<int>(std::abs(gv[1]) / b1n) + 1;
        int li = static_cast<int>(std::abs(gv[2]) / b2n) + 1;
        if (hi > max_h) max_h = hi;
        if (ki > max_k) max_k = ki;
        if (li > max_l) max_l = li;
    }
    result.nx = 2 * (max_h + 1);
    result.ny = 2 * (max_k + 1);
    result.nz = 2 * (max_l + 1);

    return result;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
