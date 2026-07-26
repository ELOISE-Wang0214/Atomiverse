/// @file tests/cpu/electronic/test_xl_bomd.cpp
/// @brief Unit tests: XlBomd / XlBomdCplx dissipative auxiliary-density
///        recursion (Niklasson JCP 130, 214109 (2009) Table I).
///
/// Mirrors include/atomiverse/cpu/electronic/xl_bomd.h — pure dense algebra on
/// the propagator class, no SCF and no GPW types. The GpwDft guess_P plumbing
/// and the GPWEvaluator build() guards live in
/// tests/cpu/gpw/test_gpw_xl_bomd.cpp; the e2e XL-BOMD NVE envelope is
/// regressions/dynamics/xl-bomd-water.

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/generators/xbuilder.hpp>

#include "atomiverse/cpu/electronic/xl_bomd.h"
#include "atomiverse/cpu/linalg/matrix.h"

using ATV::cpu::electronic::XlBomd;
using ATV::cpu::linalg::RealMatrix;

namespace {

RealMatrix constant_matrix(int n, double v) {
    RealMatrix m = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            m(i, j) = v + 0.1 * static_cast<double>(i == j);
    return m;
}

double max_abs_diff(const RealMatrix& a, const RealMatrix& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < a.shape(0); ++i)
        for (std::size_t j = 0; j < a.shape(1); ++j)
            worst = std::max(worst, std::abs(a(i, j) - b(i, j)));
    return worst;
}

} // namespace

// A1: the K=5 row of Table I, bitwise.
TEST(XlBomd, CoefficientTableK5Exact) {
    const XlBomd xl(5);
    EXPECT_EQ(xl.order(), 5);
    EXPECT_DOUBLE_EQ(xl.kappa(), 1.82);
    EXPECT_DOUBLE_EQ(xl.alpha(), 0.018);
    const std::vector<double> want = {-6.0, 14.0, -8.0, -3.0, 4.0, -1.0};
    ASSERT_EQ(xl.c().size(), want.size());
    for (std::size_t k = 0; k < want.size(); ++k)
        EXPECT_DOUBLE_EQ(xl.c()[k], want[k]) << "c[" << k << "]";
}

// A1: every Table-I row kills constant and linear history bias exactly.
TEST(XlBomd, CoefficientSumsVanishForAllOrders) {
    for (int k = 3; k <= 9; ++k) {
        const XlBomd xl(k);
        ASSERT_EQ(static_cast<int>(xl.c().size()), k + 1) << "K=" << k;
        double sum0 = 0.0, sum1 = 0.0;
        for (std::size_t j = 0; j < xl.c().size(); ++j) {
            sum0 += xl.c()[j];
            sum1 += static_cast<double>(j) * xl.c()[j];
        }
        EXPECT_DOUBLE_EQ(sum0, 0.0) << "sum c_k, K=" << k;
        EXPECT_DOUBLE_EQ(sum1, 0.0) << "sum k*c_k, K=" << k;
        EXPECT_GT(xl.kappa(), 0.0);
        EXPECT_LT(xl.kappa(), 2.0);
        EXPECT_GT(xl.alpha(), 0.0);
    }
}

// A2: only the published K range exists — no silent fallback row.
TEST(XlBomd, OrderOutsideTableThrows) {
    EXPECT_THROW(XlBomd(2), std::invalid_argument);
    EXPECT_THROW(XlBomd(10), std::invalid_argument);
}

// A2: an empty density is a programming error, not a seed.
TEST(XlBomd, AdvanceWithEmptyDensityThrows) {
    XlBomd xl(5);
    EXPECT_THROW(xl.advance(RealMatrix{}), std::invalid_argument);
}

// A2: shape changes after seeding must be surfaced (basis discontinuity).
TEST(XlBomd, AdvanceWithMismatchedShapeThrows) {
    XlBomd xl(5);
    xl.advance(constant_matrix(4, 0.3));
    EXPECT_THROW(xl.advance(constant_matrix(5, 0.3)), std::invalid_argument);
}

// A3: aux() before any density is a logic error.
TEST(XlBomd, AuxBeforeSeedThrows) {
    const XlBomd xl(5);
    EXPECT_FALSE(xl.ready());
    EXPECT_THROW(xl.aux(), std::logic_error);
}

// A3: the first advance seeds a stationary history — aux == D.
TEST(XlBomd, FirstAdvanceSeedsAuxToScfDensity) {
    XlBomd xl(5);
    const RealMatrix d = constant_matrix(4, 0.7);
    xl.advance(d);
    ASSERT_TRUE(xl.ready());
    EXPECT_DOUBLE_EQ(max_abs_diff(xl.aux(), d), 0.0);
}

// A3: a constant SCF density is the recursion's fixed point.
TEST(XlBomd, ConstantDensityIsFixedPoint) {
    XlBomd xl(5);
    const RealMatrix d = constant_matrix(6, -0.25);
    for (int n = 0; n < 21; ++n)
        xl.advance(d);
    EXPECT_LE(max_abs_diff(xl.aux(), d), 1e-12);
}

// A3: reset() drops the history (restart / discontinuity policy).
TEST(XlBomd, ResetDropsHistory) {
    XlBomd xl(5);
    xl.advance(constant_matrix(4, 0.4));
    xl.reset();
    EXPECT_FALSE(xl.ready());
    EXPECT_THROW(xl.aux(), std::logic_error);
}

// A4: bounded response to an oscillatory driving density (no secular
// growth from the dissipative Verlet recursion at K=5).
TEST(XlBomd, BoundedUnderOscillatoryDriving) {
    XlBomd xl(5);
    const int n_ao = 5;
    const RealMatrix d0 = constant_matrix(n_ao, 0.5);
    RealMatrix b = xt::zeros<double>({n_ao, n_ao});
    for (int i = 0; i < n_ao; ++i)
        for (int j = 0; j < n_ao; ++j)
            b(i, j) = ((i + j) % 2 == 0) ? 1.0 : -0.5; // symmetric, |b|max=1
    double worst = 0.0;
    for (int n = 0; n < 200; ++n) {
        const RealMatrix d =
            d0 + 0.01 * std::sin(0.3 * static_cast<double>(n)) * b;
        xl.advance(d);
        worst = std::max(worst, max_abs_diff(xl.aux(), d));
    }
    EXPECT_LE(worst, 0.1) << "auxiliary density lost the SCF density";
}

// ── multi-k XL (complex P_k) ────────────────────────────────────────────

using ATV::cpu::electronic::XlBomdCplx;
using ATV::cpu::linalg::CplxMatrix;
using Cplx = std::complex<double>;

namespace {

CplxMatrix constant_cplx(int n, Cplx v) {
    CplxMatrix m = xt::zeros<Cplx>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            m(i, j) = v + Cplx(0.1 * static_cast<double>(i == j), 0.0);
    // Hermitianize: (M + M†)/2 so fixed-point tests stay physical.
    CplxMatrix h = xt::zeros<Cplx>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            h(i, j) = 0.5 * (m(i, j) + std::conj(m(j, i)));
    return h;
}

double max_abs_diff_cplx(const CplxMatrix& a, const CplxMatrix& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < a.shape(0); ++i)
        for (std::size_t j = 0; j < a.shape(1); ++j)
            worst = std::max(worst, std::abs(a(i, j) - b(i, j)));
    return worst;
}

} // namespace

// A2: XlBomdCplx shares Table I and has a constant-density fixed point.
TEST(XlBomdCplx, ConstantHermitianIsFixedPoint) {
    XlBomdCplx xl(5);
    EXPECT_EQ(xl.order(), 5);
    EXPECT_DOUBLE_EQ(xl.kappa(), 1.82);
    EXPECT_DOUBLE_EQ(xl.alpha(), 0.018);
    const CplxMatrix d = constant_cplx(4, Cplx(0.3, 0.0));
    for (int n = 0; n < 15; ++n)
        xl.advance(d);
    EXPECT_LE(max_abs_diff_cplx(xl.aux(), d), 1e-12);
}

// A2: shape mismatch after seeding throws.
TEST(XlBomdCplx, AdvanceWithMismatchedShapeThrows) {
    XlBomdCplx xl(5);
    xl.advance(constant_cplx(3, Cplx(0.2, 0.0)));
    EXPECT_THROW(xl.advance(constant_cplx(4, Cplx(0.2, 0.0))),
                 std::invalid_argument);
}
