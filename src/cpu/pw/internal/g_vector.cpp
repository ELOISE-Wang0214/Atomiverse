/// @file cpu/pw/internal/g_vector.cpp
/// @brief G-vector enumeration implementation.

#include "atomiverse/cpu/pw/internal/g_vector.h"
#include <algorithm>
#include <cmath>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

GVectorList enumerate_g_vectors(const ::ATV::Box& cell_bohr, Real ecut,
                                const Vec3& k) {
    auto recip = cell_bohr.recip_cell();
    Vec3 b0{recip(0, 0), recip(0, 1), recip(0, 2)};
    Vec3 b1{recip(1, 0), recip(1, 1), recip(1, 2)};
    Vec3 b2{recip(2, 0), recip(2, 1), recip(2, 2)};

    Real b0n = b0.norm(), b1n = b1.norm(), b2n = b2.norm();
    Real gmax = std::sqrt(Real(2) * ecut);
    int hmax = static_cast<int>(gmax / b0n) + 1;
    int kmax = static_cast<int>(gmax / b1n) + 1;
    int lmax = static_cast<int>(gmax / b2n) + 1;

    struct GV {
        Vec3 v;
        Real e;
        int hi, ki, li;
    };
    std::vector<GV> tmp;

    for (int h = -hmax; h <= hmax; ++h) {
        for (int kk = -kmax; kk <= kmax; ++kk) {
            for (int l = -lmax; l <= lmax; ++l) {
                Vec3 g = b0 * Real(h) + b1 * Real(kk) + b2 * Real(l);
                Vec3 kpg = {g[0] + k[0], g[1] + k[1], g[2] + k[2]};
                Real e = Real(0.5) * kpg.norm_squared();
                if (e <= ecut + Real(1e-12)) {
                    tmp.push_back({g, e, h, kk, l});
                }
            }
        }
    }

    std::sort(tmp.begin(), tmp.end(), [](const GV& a, const GV& b) {
        return a.v.norm_squared() < b.v.norm_squared();
    });

    GVectorList result;
    result.g.reserve(tmp.size());
    result.ekin.reserve(tmp.size());
    result.h.reserve(tmp.size());
    result.k.reserve(tmp.size());
    result.l.reserve(tmp.size());
    for (auto& t : tmp) {
        result.g.push_back(t.v);
        result.ekin.push_back(t.e);
        result.h.push_back(t.hi);
        result.k.push_back(t.ki);
        result.l.push_back(t.li);
    }

    // FFT grid dimensions: must be even (FFTW requirement)
    int nx = 2 * (hmax + 1);
    int ny = 2 * (kmax + 1);
    int nz = 2 * (lmax + 1);
    if (nx % 2)
        ++nx;
    if (ny % 2)
        ++ny;
    if (nz % 2)
        ++nz;
    result.nx = nx;
    result.ny = ny;
    result.nz = nz;

    return result;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
