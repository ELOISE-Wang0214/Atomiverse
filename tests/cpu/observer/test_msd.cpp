/// @file tests/cpu/observer/test_msd.cpp
/// @brief MSD via the molrs FFI bridge (cpu::MolRec::compute_msd).
///
/// Validates the molrs-first analysis chain: stage_raw -> commit_frame ->
/// molrs compute::MSD -> curve back to C++. All analysis math is in molrs; this
/// asserts the FFI round-trip reproduces the analytic MSD.

#include "atomiverse/cpu/io/molrec.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

using ATV::cpu::MolRec;

namespace {
std::array<double, 9> cubic(double L) {
    return {L, 0, 0, 0, L, 0, 0, 0, L};
}
} // namespace

// Two atoms with known linear displacement: atom 0 moves +a per frame along x,
// atom 1 moves +b per frame along y. Direct-mode MSD relative to frame 0:
//   MSD(t) = <|r(t) - r(0)|^2>_i = ( (a t)^2 + (b t)^2 ) / 2 = t^2 (a^2+b^2)/2.
TEST(CpuMSD, LinearDisplacementMatchesAnalytic) {
    MolRec rec;
    const double a = 1.0, b = 2.0, L = 100.0; // box large -> no PBC crossing
    const int n_frames = 4;
    const std::vector<std::int32_t> type{1, 1};

    for (int t = 0; t < n_frames; ++t) {
        std::vector<double> x{a * t, 0.0};
        std::vector<double> y{0.0, b * t};
        std::vector<double> z{0.0, 0.0};
        rec.stage_raw(type, x, y, z, cubic(L));
        rec.commit_frame();
    }
    ASSERT_EQ(rec.n_frames(), n_frames);

    std::vector<double> msd = rec.compute_msd();
    ASSERT_EQ(msd.size(), static_cast<std::size_t>(n_frames));

    for (int t = 0; t < n_frames; ++t) {
        const double expected = static_cast<double>(t) * t * (a * a + b * b) / 2.0;
        EXPECT_NEAR(msd[t], expected, 1e-9) << "frame t=" << t;
    }
    EXPECT_DOUBLE_EQ(msd[0], 0.0); // reference frame
}

// Fewer than two frames -> empty curve, never a panic/abort across the FFI.
TEST(CpuMSD, SingleFrameReturnsEmpty) {
    MolRec rec;
    std::vector<double> x{0.0}, y{0.0}, z{0.0};
    rec.stage_raw(std::vector<std::int32_t>{1}, x, y, z, cubic(10.0));
    rec.commit_frame();
    EXPECT_TRUE(rec.compute_msd().empty());
}

TEST(CpuMSD, NoFramesReturnsEmpty) {
    MolRec rec;
    EXPECT_TRUE(rec.compute_msd().empty());
}

// Einstein-diffusion bridge (molrs EinsteinDiffusion + LinearFit). molrs owns
// the numerical correctness (its own Rust tests); here we validate the FFI
// round-trip: a valid trajectory yields a finite D, degenerate input yields NaN.
TEST(CpuMSD, DiffusionBridgeFiniteForValidTrajectory) {
    MolRec rec;
    const std::vector<std::int32_t> type{1, 1};
    for (int t = 0; t < 6; ++t) {
        std::vector<double> x{1.0 * t, 0.0};
        std::vector<double> y{0.0, 2.0 * t};
        std::vector<double> z{0.0, 0.0};
        rec.stage_raw(type, x, y, z, cubic(100.0));
        rec.commit_frame();
    }
    const double d = rec.compute_diffusion(/*dt=*/1.0, /*dims=*/3, 0.2, 0.8);
    EXPECT_TRUE(std::isfinite(d)) << "D must be finite for a valid trajectory";
}

TEST(CpuMSD, DiffusionEmptyIsNaN) {
    MolRec rec;
    EXPECT_TRUE(std::isnan(rec.compute_diffusion(1.0)));
}
