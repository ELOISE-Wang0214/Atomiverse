/// @file cpu/pw/internal/kpoint.cpp
/// @brief Monkhorst-Pack k-point mesh implementation.

#include "atomiverse/cpu/pw/internal/kpoint.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

KPointGrid generate_mp_mesh(int nk1, int nk2, int nk3,
                            Real shift1, Real shift2, Real shift3) {
    KPointGrid kg;
    int nk = nk1 * nk2 * nk3;
    kg.kpoints.reserve(nk);
    kg.weights.resize(nk, Real(1) / Real(nk));

    Real dn1 = Real(1) / Real(nk1), dn2 = Real(1) / Real(nk2), dn3 = Real(1) / Real(nk3);
    for (int i = 0; i < nk1; ++i)
        for (int j = 0; j < nk2; ++j)
            for (int k = 0; k < nk3; ++k) {
                kg.kpoints.push_back({
                    (Real(i) + shift1) * dn1,
                    (Real(j) + shift2) * dn2,
                    (Real(k) + shift3) * dn3
                });
            }
    return kg;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
