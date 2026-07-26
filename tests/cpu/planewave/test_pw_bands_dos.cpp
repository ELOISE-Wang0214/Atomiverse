/// @file test_pw_bands_dos.cpp
/// @brief Unit tests for band-path geometry, gap finder, Gaussian DOS.
///        No full SCF e2e (that lives in regressions/pw/bands-dos-si).

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/planewave/bands_dos.h"

using namespace ATV::cpu::pw;

TEST(BandsPath, DiamondPathIncludesGammaAndX) {
    Bands bands;
    bands.set_diamond_path(4);
    ASSERT_GE(bands.kfrac.size(), 5u);

    auto near = [](const std::array<double, 3>& a,
                   const std::array<double, 3>& b, double tol = 1e-12) {
        return std::abs(a[0] - b[0]) < tol && std::abs(a[1] - b[1]) < tol
               && std::abs(a[2] - b[2]) < tol;
    };

    const std::array<double, 3> G{0, 0, 0};
    const std::array<double, 3> X{0.5, 0.0, 0.5};
    bool has_G = false, has_X = false;
    for (const auto& k : bands.kfrac) {
        if (near(k, G))
            has_G = true;
        if (near(k, X))
            has_X = true;
    }
    EXPECT_TRUE(has_G);
    EXPECT_TRUE(has_X);

    // n_seg=4 → 4 samples per segment; 5 segments with shared joins:
    // start + 5*(4-1) = 1 + 15 = 16.
    EXPECT_EQ(bands.kfrac.size(), 16u);
}

TEST(GapFromSpectrum, SyntheticIndirectGap) {
    // 2 k-points, 3 bands; n_occ=2.
    // k0: -1.0, -0.5, 0.8   (VBM candidate -0.5)
    // k1: -0.9, -0.2, 0.3   (VBM -0.2, CBM 0.3) → VBM at k1, CBM at k1
    // Make CBM at different k:
    // k0: -1, -0.1, 1.0
    // k1: -1, -0.5, 0.2
    std::vector<std::vector<double>> evals = {
        {-1.0, -0.1, 1.0},
        {-1.0, -0.5, 0.2},
    };
    double e_vbm = 0, e_cbm = 0;
    int vk = -1, ck = -1, vb = -1, cb = -1;
    const double gap =
        gap_from_spectrum(evals, /*n_occ=*/2, e_vbm, e_cbm, vk, ck, vb, cb);
    EXPECT_NEAR(e_vbm, -0.1, 1e-12);
    EXPECT_NEAR(e_cbm, 0.2, 1e-12);
    EXPECT_NEAR(gap, 0.3, 1e-12);
    EXPECT_EQ(vk, 0);
    EXPECT_EQ(ck, 1);
    EXPECT_EQ(vb, 1);
    EXPECT_EQ(cb, 2);
}

TEST(GaussianDos, SinglePeakIntegratesToWeightTimesTwo) {
    // One eigenvalue at 0 with weight 1 → ∫ D ≈ 2 (closed-shell factor).
    std::vector<double> evals = {0.0};
    std::vector<double> weights = {1.0};
    std::vector<double> energy, dos;
    const double sigma = 0.05;
    gaussian_dos(evals, weights, /*e_min=*/-1.0, /*e_max=*/1.0, /*n_bins=*/2001,
                 sigma, energy, dos);

    double integ = 0;
    const double de = energy[1] - energy[0];
    for (std::size_t i = 0; i + 1 < energy.size(); ++i)
        integ += 0.5 * (dos[i] + dos[i + 1]) * de;
    EXPECT_NEAR(integ, 2.0, 1e-3);
}

TEST(BandsResult, DefaultIsRestricted) {
    BandsResult br;
    EXPECT_EQ(br.nspin, 1);
    EXPECT_TRUE(br.eigenvalues_beta_ry.empty());
    EXPECT_NEAR(br.magnetization, 0.0, 1e-15);
    EXPECT_EQ(br.vbm_spin, 0);
    EXPECT_EQ(br.cbm_spin, 0);
}
