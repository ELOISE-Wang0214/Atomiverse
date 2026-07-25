/// @file test_berry_polarization.cpp
/// @brief Unit tests for Berry-phase string geometry, link-det phase math
///        and quantum wrapping. No SCF e2e (that lives in
///        regressions/pw/berry-polarization).

#include <array>
#include <cmath>
#include <complex>
#include <set>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/planewave/berry_polarization.h"
#include "atomiverse/cpu/planewave/internal/pw_subspace.h"

using namespace ATV::cpu::pw;

namespace {
constexpr double kPi = 3.14159265358979323846;
}

TEST(BerryStrings, LayoutStringMajorAndSpacing) {
    BerryPolarization berry;
    berry.k_perp = 2;
    berry.k_string = 4;

    const auto kf = berry.string_kpoints(/*dir=*/2);
    ASSERT_EQ(kf.size(), 2u * 2u * 4u);

    // String-major: within each block of 4, the transverse components are
    // constant and the dir component walks j/J.
    for (int s = 0; s < 4; ++s) {
        for (int j = 0; j < 4; ++j) {
            const auto& f = kf[static_cast<std::size_t>(s * 4 + j)];
            EXPECT_NEAR(f[2], static_cast<double>(j) / 4.0, 1e-12);
            EXPECT_NEAR(f[0], kf[static_cast<std::size_t>(s * 4)][0], 1e-12);
            EXPECT_NEAR(f[1], kf[static_cast<std::size_t>(s * 4)][1], 1e-12);
            for (int d = 0; d < 3; ++d) {
                EXPECT_GE(f[static_cast<std::size_t>(d)], 0.0);
                EXPECT_LT(f[static_cast<std::size_t>(d)], 1.0);
            }
        }
    }

    // Transverse offsets are the full i/N grid: 4 distinct (f0, f1) pairs.
    std::set<std::pair<double, double>> perp;
    for (int s = 0; s < 4; ++s)
        perp.insert({kf[static_cast<std::size_t>(s * 4)][0],
                     kf[static_cast<std::size_t>(s * 4)][1]});
    EXPECT_EQ(perp.size(), 4u);
    EXPECT_TRUE(perp.count({0.0, 0.0}) == 1);
    EXPECT_TRUE(perp.count({0.5, 0.5}) == 1);
}

TEST(BerryStrings, DirectionAxisSelects) {
    BerryPolarization berry;
    berry.k_perp = 1;
    berry.k_string = 3;
    for (int dir = 0; dir < 3; ++dir) {
        const auto kf = berry.string_kpoints(dir);
        ASSERT_EQ(kf.size(), 3u);
        for (int j = 0; j < 3; ++j) {
            for (int d = 0; d < 3; ++d) {
                const double want =
                    (d == dir) ? static_cast<double>(j) / 3.0 : 0.0;
                EXPECT_NEAR(kf[static_cast<std::size_t>(j)]
                              [static_cast<std::size_t>(d)],
                            want, 1e-12);
            }
        }
    }
}

TEST(BerryStrings, InvalidArgumentsThrow) {
    BerryPolarization berry;
    EXPECT_THROW(berry.string_kpoints(-1), std::invalid_argument);
    EXPECT_THROW(berry.string_kpoints(3), std::invalid_argument);
    berry.k_perp = 0;
    EXPECT_THROW(berry.string_kpoints(0), std::invalid_argument);
    berry.k_perp = 2;
    berry.k_string = 1;
    EXPECT_THROW(berry.string_kpoints(0), std::invalid_argument);
}

TEST(BerryPhaseMath, ProductOfUnitPhases) {
    // Three links of e^{i pi/6} → phase pi/2.
    const std::complex<double> w = std::polar(1.0, kPi / 6.0);
    const double phi = phase_from_link_dets({w, w, w});
    EXPECT_NEAR(phi, kPi / 2.0, 1e-12);
}

TEST(BerryPhaseMath, ModulusIsNormalizedAway) {
    // Moduli != 1 must not change the phase.
    const std::vector<std::complex<double>> dets = {
        std::polar(0.3, 0.4), std::polar(7.0, -1.1), std::polar(0.05, 0.9)};
    EXPECT_NEAR(phase_from_link_dets(dets), 0.4 - 1.1 + 0.9, 1e-12);
}

TEST(BerryPhaseMath, WrapsIntoMinusPiPi) {
    // Total 3*(2 pi/3) = 2 pi → wraps to 0.
    const std::complex<double> w = std::polar(1.0, 2.0 * kPi / 3.0);
    EXPECT_NEAR(phase_from_link_dets({w, w, w}), 0.0, 1e-12);
    // Total 4*(3 pi/4) = 3 pi → wraps to pi (boundary is +pi, not -pi).
    const std::complex<double> v = std::polar(1.0, 3.0 * kPi / 4.0);
    EXPECT_NEAR(std::abs(phase_from_link_dets({v, v, v, v})), kPi, 1e-12);
}

TEST(BerryPhaseMath, ThrowsOnDeadLinkOrEmpty) {
    EXPECT_THROW(phase_from_link_dets({}), std::invalid_argument);
    EXPECT_THROW(phase_from_link_dets({{1.0, 0.0}, {0.0, 0.0}}),
                 std::invalid_argument);
    EXPECT_THROW(phase_from_link_dets({{5e-13, 0.0}}), std::invalid_argument);
}

TEST(BerryWrap, CenteredInterval) {
    EXPECT_NEAR(wrap_centered(0.3, 1.0), 0.3, 1e-12);
    EXPECT_NEAR(wrap_centered(0.7, 1.0), -0.3, 1e-12);
    EXPECT_NEAR(wrap_centered(-0.7, 1.0), 0.3, 1e-12);
    EXPECT_NEAR(wrap_centered(2.0, 2.0), 0.0, 1e-12);
    EXPECT_NEAR(wrap_centered(1.9, 2.0), -0.1, 1e-12);
    // Boundary: exactly -q/2 maps to +q/2.
    EXPECT_NEAR(wrap_centered(-1.0, 2.0), 1.0, 1e-12);
    EXPECT_NEAR(wrap_centered(1.0, 2.0), 1.0, 1e-12);
}

TEST(BerryWrap, ThrowsOnNonPositiveQuantum) {
    EXPECT_THROW(wrap_centered(0.1, 0.0), std::invalid_argument);
    EXPECT_THROW(wrap_centered(0.1, -1.0), std::invalid_argument);
}

TEST(JacobiHermitian, TwoByTwoKnownSpectrum) {
    // [[1, i], [-i, 1]] has eigenvalues {0, 2}.
    using C = std::complex<ATV::AccReal>;
    std::vector<C> A = {C(1, 0), C(0, 1), C(0, -1), C(1, 0)};
    ATV::AccReal eig[2];
    std::vector<C> V(4);
    ATV::cpu::pw::internal::jacobi_eigh_hermitian(2, A.data(), eig, V.data());
    EXPECT_NEAR(static_cast<double>(eig[0]), 0.0, 1e-12);
    EXPECT_NEAR(static_cast<double>(eig[1]), 2.0, 1e-12);
}

TEST(JacobiHermitian, ResidualOrderAndUnitarity) {
    using C = std::complex<ATV::AccReal>;
    const int n = 5;
    // Deterministic pseudo-random Hermitian matrix.
    std::vector<C> A0(static_cast<std::size_t>(n * n));
    unsigned rng = 42u;
    auto next = [&rng]() {
        rng = rng * 1664525u + 1013904223u;
        return (static_cast<double>(rng % 10007u) / 10007.0) - 0.5;
    };
    for (int i = 0; i < n; ++i) {
        A0[static_cast<std::size_t>(i * n + i)] = C(next() * 4.0, 0);
        for (int j = i + 1; j < n; ++j) {
            const C z(next(), next());
            A0[static_cast<std::size_t>(i * n + j)] = z;
            A0[static_cast<std::size_t>(j * n + i)] = std::conj(z);
        }
    }
    std::vector<C> A = A0;
    std::vector<ATV::AccReal> eig(static_cast<std::size_t>(n));
    std::vector<C> V(static_cast<std::size_t>(n * n));
    ATV::cpu::pw::internal::jacobi_eigh_hermitian(n, A.data(), eig.data(),
                                                  V.data());

    // Ascending order.
    for (int k = 1; k < n; ++k)
        EXPECT_LE(eig[static_cast<std::size_t>(k - 1)],
                  eig[static_cast<std::size_t>(k)]);

    // Columns orthonormal and residual ||A0 v - lambda v|| ~ 0.
    for (int k = 0; k < n; ++k) {
        double nrm = 0;
        for (int i = 0; i < n; ++i)
            nrm += std::norm(V[static_cast<std::size_t>(i * n + k)]);
        EXPECT_NEAR(nrm, 1.0, 1e-12);
        double res = 0;
        for (int i = 0; i < n; ++i) {
            C av(0, 0);
            for (int j = 0; j < n; ++j)
                av += A0[static_cast<std::size_t>(i * n + j)]
                      * V[static_cast<std::size_t>(j * n + k)];
            av -= eig[static_cast<std::size_t>(k)]
                  * V[static_cast<std::size_t>(i * n + k)];
            res += std::norm(av);
        }
        EXPECT_LT(std::sqrt(res), 1e-10);
    }
    // Distinct columns orthogonal.
    for (int k = 0; k < n; ++k)
        for (int l = k + 1; l < n; ++l) {
            C dot(0, 0);
            for (int i = 0; i < n; ++i)
                dot += std::conj(V[static_cast<std::size_t>(i * n + k)])
                       * V[static_cast<std::size_t>(i * n + l)];
            EXPECT_LT(std::abs(dot), 1e-11);
        }
}
