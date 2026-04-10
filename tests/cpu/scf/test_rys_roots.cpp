/// @file tests/cpu/scf/test_rys_roots.cpp
/// @brief Byte-level validation of native::rys_roots against libcint's
///        CINTrys_roots across the domain used by CPU SCF.
///
/// The native ERI engine delegates Boys→roots conversion to a ported copy
/// of libcint's closed-form rys_root1..5. This test asserts that our copy
/// matches libcint element-wise at floating-point precision across:
///   - nroots ∈ {1..5}
///   - x  ∈ {0, small, near 1, medium, near branch boundaries, large}
/// Any drift here implies an incomplete port, so tolerance is tight.

#if defined(ATV_HAS_RHF) && defined(ATV_HAS_LIBCINT)

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "atomiverse/cpu/scf/rys_roots.h"

extern "C" {
void CINTrys_roots(int nroots, double x, double *u, double *w);
}

namespace {

// Sample x at enough branch boundaries to cover every piecewise polynomial
// segment inside rys_root1..5. Values straddle libcint's branch points
// (SMALLX_LIMIT=3e-7, X<1, X<3, X<5, X<10, X<15, X<33, X<40, X>=40).
const std::vector<double> kSampleX = {
    0.0,            // tails into the small-x polynomial.
    1e-10,          // clearly below SMALLX_LIMIT.
    1e-8,           // still small-x.
    5e-7,           // just above SMALLX_LIMIT (first piecewise branch).
    1e-4,           // small polynomial branch.
    0.3,  0.7,      // X<1 branch.
    1.0,  1.5,  2.0, 2.5,  // X<3 branch.
    3.0,  4.0,      // X<5 branch.
    5.5,  7.5,  9.5,       // X<10 branch.
    10.5, 12.5, 14.5,      // X<15 branch.
    16.0, 20.0, 25.0, 30.0, 32.0,   // X<33 branch.
    33.5, 37.0, 39.0,              // X<40 / X<47 branches.
    41.0, 50.0, 80.0,              // large-x branches.
    120.0,                         // deep large-x polynomial cutoff.
};

constexpr double kTol = 2e-14;  // element-wise absolute tolerance vs libcint.

} // namespace

TEST(RysRoots, MatchesLibcint) {
    for (int nroots = 1; nroots <= 5; ++nroots) {
        for (double x : kSampleX) {
            double u_ref[8] = {};
            double w_ref[8] = {};
            double u_atv[8] = {};
            double w_atv[8] = {};

            CINTrys_roots(nroots, x, u_ref, w_ref);
            ATV::cpu::scf::native::rys_roots(nroots, x, u_atv, w_atv);

            for (int i = 0; i < nroots; ++i) {
                EXPECT_NEAR(u_atv[i], u_ref[i], kTol)
                    << "root drift at nroots=" << nroots << " x=" << x
                    << " i=" << i;
                EXPECT_NEAR(w_atv[i], w_ref[i], kTol)
                    << "weight drift at nroots=" << nroots << " x=" << x
                    << " i=" << i;
            }
        }
    }
}

#endif // ATV_HAS_RHF && ATV_HAS_LIBCINT
