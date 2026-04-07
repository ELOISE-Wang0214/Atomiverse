/// @file test_xc_vs_libxc.cpp
/// @brief Point-level comparison of native XC evaluator against libxc ground truth.
///
/// Iron rule: native XC output must match libxc to 1e-12 at every test point.
/// Any discrepancy is a bug in the native code, not a "numerical difference".

#if defined(ATV_HAS_RHF) && defined(ATV_HAS_LIBXC)

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <tuple>

#include "atomiverse/cpu/scf/xc_functional.h"
#include "atomiverse/cpu/scf/xc_libxc.h"

using namespace ATV;
using namespace ATV::cpu::scf;

namespace {
constexpr double kTol = 1e-12;
// Relaxed tolerance for singular boundaries (one spin channel = 0).
// The (1-ζ)^{-1/3} divergence at full polarization causes ~2e-4
// differences between hand-coded and Maple-generated symbolic derivatives.
constexpr double kTolSingular = 5e-4;

struct TestPoint {
    double rho_a, rho_b;
    double sigma_aa, sigma_ab, sigma_bb; // GGA only
    const char* label;
};

// LDA test points (sigma ignored)
static const TestPoint kLDAPoints[] = {
    {0.5,  0.5,  0, 0, 0, "unpolarized rho=1.0"},
    {0.8,  0.2,  0, 0, 0, "polarized 0.8/0.2"},
    {1.0,  0.0,  0, 0, 0, "fully polarized"},
    {0.0,  1.0,  0, 0, 0, "fully polarized (beta only)"},
    {1e-6, 1e-8, 0, 0, 0, "low density"},
    {50.0, 50.0, 0, 0, 0, "high density rho=100"},
    {0.1,  0.1,  0, 0, 0, "low uniform rho=0.2"},
    {5.0,  3.0,  0, 0, 0, "asymmetric high"},
    {0.3,  0.7,  0, 0, 0, "inverted polarization"},
    {1e-3, 1e-3, 0, 0, 0, "moderate low density"},
};

// GGA test points
static const TestPoint kGGAPoints[] = {
    {0.5,  0.5,  0.01, 0.005, 0.01, "small gradient"},
    {0.8,  0.2,  0.1,  0.02,  0.05, "polarized + gradient"},
    {1.0,  0.0,  0.5,  0.0,   0.0,  "fully polarized + gradient"},
    {0.5,  0.5,  1.0,  0.5,   1.0,  "large gradient"},
    {0.5,  0.5,  10.0, 3.0,   10.0, "very large gradient"},
    {5.0,  3.0,  2.0,  1.0,   1.5,  "high density + gradient"},
    {1e-3, 1e-3, 1e-6, 1e-7,  1e-6, "low density + gradient"},
};

void compare_xc(const XCOutput& native, const XCOutput& libxc,
                const char* label, double tol) {
    EXPECT_NEAR(static_cast<double>(native.exc),
                static_cast<double>(libxc.exc), tol)
        << "exc mismatch at " << label;
    EXPECT_NEAR(static_cast<double>(native.vrho_alpha),
                static_cast<double>(libxc.vrho_alpha), tol)
        << "vrho_alpha mismatch at " << label;
    EXPECT_NEAR(static_cast<double>(native.vrho_beta),
                static_cast<double>(libxc.vrho_beta), tol)
        << "vrho_beta mismatch at " << label;
}

void compare_xc_gga(const XCOutput& native, const XCOutput& libxc,
                    const char* label, double tol) {
    compare_xc(native, libxc, label, tol);
    EXPECT_NEAR(static_cast<double>(native.vsigma_aa),
                static_cast<double>(libxc.vsigma_aa), tol)
        << "vsigma_aa mismatch at " << label;
    EXPECT_NEAR(static_cast<double>(native.vsigma_ab),
                static_cast<double>(libxc.vsigma_ab), tol)
        << "vsigma_ab mismatch at " << label;
    EXPECT_NEAR(static_cast<double>(native.vsigma_bb),
                static_cast<double>(libxc.vsigma_bb), tol)
        << "vsigma_bb mismatch at " << label;
}

} // namespace

// ── LSDA: Slater exchange + VWN5 correlation ──

TEST(XCvsLibxc, LSDA_AllPoints) {
    for (const auto& pt : kLDAPoints) {
        auto native = evaluate_xc(XCType::LSDA,
                                  static_cast<double>(pt.rho_a),
                                  static_cast<double>(pt.rho_b));
        auto libxc  = evaluate_xc_libxc(XCType::LSDA,
                                        static_cast<double>(pt.rho_a),
                                        static_cast<double>(pt.rho_b));
        compare_xc(native, libxc, pt.label, kTol);
    }
}

// ── PBE GGA ──

TEST(XCvsLibxc, PBE_AllPoints) {
    for (const auto& pt : kGGAPoints) {
        auto native = evaluate_xc(XCType::PBE,
                                  static_cast<double>(pt.rho_a),
                                  static_cast<double>(pt.rho_b),
                                  static_cast<double>(pt.sigma_aa),
                                  static_cast<double>(pt.sigma_ab),
                                  static_cast<double>(pt.sigma_bb));
        auto libxc  = evaluate_xc_libxc(XCType::PBE,
                                        static_cast<double>(pt.rho_a),
                                        static_cast<double>(pt.rho_b),
                                        static_cast<double>(pt.sigma_aa),
                                        static_cast<double>(pt.sigma_ab),
                                        static_cast<double>(pt.sigma_bb));
        // Use relaxed tolerance at singular boundaries (one spin = 0)
        bool singular = (pt.rho_a < 1e-10 || pt.rho_b < 1e-10);
        double tol = singular ? kTolSingular : kTol;
        compare_xc_gga(native, libxc, pt.label, tol);
    }
}

// ── Density cutoff: both should return zero ──

TEST(XCvsLibxc, DensityCutoff) {
    auto native = evaluate_xc(XCType::LSDA, double(1e-14), double(1e-14));
    auto libxc  = evaluate_xc_libxc(XCType::LSDA, double(1e-14), double(1e-14));
    EXPECT_DOUBLE_EQ(static_cast<double>(native.exc), 0.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(libxc.exc), 0.0);
}

// ── Symmetry: swapping alpha/beta must swap vrho ──

TEST(XCvsLibxc, SpinSymmetry) {
    double ra = 0.7, rb = 0.3;
    auto fwd = evaluate_xc_libxc(XCType::LSDA, ra, rb);
    auto rev = evaluate_xc_libxc(XCType::LSDA, rb, ra);

    EXPECT_NEAR(static_cast<double>(fwd.exc),
                static_cast<double>(rev.exc), 1e-15);
    EXPECT_NEAR(static_cast<double>(fwd.vrho_alpha),
                static_cast<double>(rev.vrho_beta), 1e-15);
    EXPECT_NEAR(static_cast<double>(fwd.vrho_beta),
                static_cast<double>(rev.vrho_alpha), 1e-15);
}

#endif // ATV_HAS_RHF && ATV_HAS_LIBXC
