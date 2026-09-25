/// @file tests/test_molrs/test_observe.cpp
/// @brief Drift fence for molrs trajectory-analysis FFI via ATV wrappers.

#include <array>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/observe/trajectory_metrics.h"

using ATV::cpu::compute::MSD;
using ATV::cpu::compute::MSDAccumulator;
using ATV::cpu::compute::RDF;
using ATV::cpu::compute::RDFAccumulator;
using ATV::cpu::compute::VACF;
using ATV::cpu::compute::VACFAccumulator;

namespace {

void push_xyz(std::vector<double>& buf, double x, double y, double z) {
    buf.push_back(x);
    buf.push_back(y);
    buf.push_back(z);
}

TEST(MolrsObserve, MsdBatchTwoFramesFinite) {
    std::vector<double> pos;
    push_xyz(pos, 0, 0, 0);
    push_xyz(pos, 1, 0, 0);
    auto curve = MSD{}.compute(pos, /*n_frames=*/2, /*n_dof=*/3);
    ASSERT_EQ(curve.size(), 2u);
    EXPECT_DOUBLE_EQ(curve[0], 0.0);
    EXPECT_NEAR(curve[1], 1.0, 1e-12);
}

TEST(MolrsObserve, MsdAccumulatorMatchesBatch) {
    MSDAccumulator acc(/*window=*/0);
    EXPECT_TRUE(acc.accumulate({0, 0, 0}));
    EXPECT_TRUE(acc.accumulate({1, 0, 0}));
    EXPECT_EQ(acc.n_frames(), 2);
    auto curve = acc.direct_curve();
    ASSERT_EQ(curve.size(), 2u);
    EXPECT_NEAR(curve[1], 1.0, 1e-12);
}

TEST(MolrsObserve, VacfBatchTwoFramesNonEmpty) {
    std::vector<double> vel;
    push_xyz(vel, 1, 0, 0);
    push_xyz(vel, 1, 0, 0);
    auto curve = VACF(/*dt=*/1.0, /*resolution=*/2).compute(vel, 2, 3);
    EXPECT_FALSE(curve.empty());
}

TEST(MolrsObserve, VacfAccumulatorFinalize) {
    VACFAccumulator acc(/*resolution=*/2);
    EXPECT_TRUE(acc.accumulate({1, 0, 0}));
    EXPECT_TRUE(acc.accumulate({1, 0, 0}));
    auto curve = acc.finalize();
    EXPECT_FALSE(curve.empty());
}

TEST(MolrsObserve, RdfBatchOneFrameFiniteBins) {
    // Two particles in a large box; one bin curve comes back.
    std::vector<double> pos;
    push_xyz(pos, 0, 0, 0);
    push_xyz(pos, 1.5, 0, 0);
    const std::vector<double> box = {20, 0, 0, 0, 20, 0, 0, 0, 20};
    auto g = RDF(/*n_bins=*/5, /*r_max=*/5.0, /*r_min=*/0.0)
                 .compute(pos, box, /*n_frames=*/1, /*n_atoms=*/2);
    ASSERT_EQ(g.size(), 5u);
    for (double v : g)
        EXPECT_TRUE(std::isfinite(v));
}

TEST(MolrsObserve, RdfAccumulatorFinalize) {
    RDFAccumulator acc(/*n_bins=*/4, /*r_max=*/4.0, /*r_min=*/0.0);
    const std::array<double, 9> cell = {20, 0, 0, 0, 20, 0, 0, 0, 20};
    EXPECT_TRUE(acc.accumulate({0, 0, 0, 1.0, 0, 0}, cell));
    EXPECT_EQ(acc.n_frames(), 1);
    auto g = acc.finalize();
    ASSERT_EQ(g.size(), 4u);
}

} // namespace
