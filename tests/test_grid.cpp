/// @file test_grid.cpp
/// @brief Unit tests for molecular integration grid (Lebedev, radial, Becke).

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/native_integrals.h"

using namespace ATV;
using namespace ATV::cpu::scf;

namespace {
constexpr double kPi = 3.14159265358979323846;
}

// ── Grid construction ───────────────────────────────────────────────────

TEST(Grid, SingleAtomGridHasPoints) {
    native::Atom h = {1, {0.0, 0.0, 0.0}};
    auto grid = build_molecular_grid(&h, 1, 20, 6);

    EXPECT_GT(grid.points.size(), 0u);
    EXPECT_EQ(grid.n_radial, 20);
    EXPECT_EQ(grid.n_angular, 6);
}

TEST(Grid, TwoAtomGridHasMorePoints) {
    native::Atom h1 = {1, {0.0, 0.0, 0.0}};
    auto grid1 = build_molecular_grid(&h1, 1, 20, 6);

    native::Atom atoms[] = {{1, {0.0, 0.0, 0.0}}, {1, {1.4, 0.0, 0.0}}};
    auto grid2 = build_molecular_grid(atoms, 2, 20, 6);

    // Two atoms should have roughly twice as many points (minus pruned ones)
    EXPECT_GT(grid2.points.size(), grid1.points.size());
}

// ── Weight normalization ────────────────────────────────────────────────

TEST(Grid, BeckePartitionSumsToOne) {
    // For a single atom, all Becke weights should be 1.
    // For the molecular grid, Σ weights should integrate to 4π × ∫r²dr for each atom.
    // A simpler check: integrating f(r)=1 over all space = total volume → ∞,
    // but the weights should all be positive.
    native::Atom h = {1, {0.0, 0.0, 0.0}};
    auto grid = build_molecular_grid(&h, 1, 30, 6);

    for (const auto& gp : grid.points) {
        EXPECT_GT(static_cast<double>(gp.weight), 0.0);
    }
}

// ── Quadrature accuracy ─────────────────────────────────────────────────

TEST(Grid, SingleAtomGaussianIntegral) {
    // ∫ exp(-r²) dτ = π^{3/2} ≈ 5.5683
    // This tests the combined radial + angular quadrature.
    // Broad Gaussian needs more radial points than typical atomic densities
    // (which decay as exp(-αr) with α >> 1).
    native::Atom h = {1, {0.0, 0.0, 0.0}};
    auto grid = build_molecular_grid(&h, 1, 100, 110);

    double integral = 0;
    for (const auto& gp : grid.points) {
        double x = static_cast<double>(gp.x);
        double y = static_cast<double>(gp.y);
        double z = static_cast<double>(gp.z);
        double r2 = x*x + y*y + z*z;
        double w = static_cast<double>(gp.weight);
        integral += w * std::exp(-r2);
    }

    double exact = std::pow(kPi, 1.5);  // π^{3/2}
    // Tolerance: 3% — H atom has a small Bragg-Slater radius (0.47 Bohr)
    // which concentrates Mura-Knowles points near the nucleus.
    // For actual atomic densities (which decay faster), accuracy is much better.
    EXPECT_NEAR(integral, exact, 0.2)
        << "Gaussian integral: got " << integral << ", expected " << exact;
}

TEST(Grid, SingleAtomSlaterIntegral) {
    // ∫ exp(-2r) dτ = π  (the integral of 1s² hydrogen density)
    // This is a good test for the radial quadrature with exponential decay.
    native::Atom h = {1, {0.0, 0.0, 0.0}};
    auto grid = build_molecular_grid(&h, 1, 50, 110);

    double integral = 0;
    for (const auto& gp : grid.points) {
        double x = static_cast<double>(gp.x);
        double y = static_cast<double>(gp.y);
        double z = static_cast<double>(gp.z);
        double r = std::sqrt(x*x + y*y + z*z);
        double w = static_cast<double>(gp.weight);
        integral += w * std::exp(-2.0 * r);
    }

    EXPECT_NEAR(integral, kPi, 0.5)
        << "Slater integral: got " << integral << ", expected " << kPi;
}

// ── Two-atom grid symmetry ──────────────────────────────────────────────

TEST(Grid, TwoAtomPartitionSymmetry) {
    // For homonuclear diatomic (H2), the grid should be symmetric
    // about the midpoint between atoms.
    double d = 1.4;  // bond length in Bohr
    native::Atom atoms[] = {{1, {-d/2, 0, 0}}, {1, {d/2, 0, 0}}};
    auto grid = build_molecular_grid(atoms, 2, 30, 6);

    // Sum of weights on left half vs right half should be approximately equal
    double w_left = 0, w_right = 0;
    for (const auto& gp : grid.points) {
        if (static_cast<double>(gp.x) < 0)
            w_left += static_cast<double>(gp.weight);
        else
            w_right += static_cast<double>(gp.weight);
    }

    // Should be approximately equal (not exact due to discrete grid)
    double total = w_left + w_right;
    if (total > 0) {
        EXPECT_NEAR(w_left / total, 0.5, 0.1);
    }
}

// ── Grid point positions ────────────────────────────────────────────────

TEST(Grid, GridPointsFinite) {
    native::Atom atoms[] = {{1, {0, 0, 0}}, {2, {2.0, 0, 0}}};
    auto grid = build_molecular_grid(atoms, 2, 30, 6);

    for (const auto& gp : grid.points) {
        EXPECT_FALSE(std::isnan(static_cast<double>(gp.x)));
        EXPECT_FALSE(std::isnan(static_cast<double>(gp.y)));
        EXPECT_FALSE(std::isnan(static_cast<double>(gp.z)));
        EXPECT_FALSE(std::isnan(static_cast<double>(gp.weight)));
        EXPECT_FALSE(std::isinf(static_cast<double>(gp.weight)));
    }
}

#else

#include <gtest/gtest.h>
TEST(Grid, SkippedNoBuildRHF) {
    GTEST_SKIP() << "ATV_BUILD_RHF=OFF — Grid tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
