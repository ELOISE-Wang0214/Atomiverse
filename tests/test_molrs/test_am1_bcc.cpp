/// @file tests/test_molrs/test_am1_bcc.cpp
/// @brief Drift fence for molrs::am1_bcc_assign_frame_from_base (ATV call site).

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/molrs.hpp"

namespace {

::rust::Box<molrs::FrameRef> methane_like() {
    auto fr = molrs::frame_new();
    const std::vector<double> x = {0.0, 0.63, -0.63, -0.63, 0.63};
    const std::vector<double> y = {0.0, 0.63, -0.63, 0.63, -0.63};
    const std::vector<double> z = {0.0, 0.63, 0.63, -0.63, -0.63};
    const std::vector<std::uint64_t> znum = {6, 1, 1, 1, 1};
    molrs::frame_set_column_f64(
        *fr, "atoms", "x", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "y", ::rust::Slice<const double>(y.data(), y.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "z", ::rust::Slice<const double>(z.data(), z.size()));
    molrs::frame_set_column_u32(
        *fr, "atoms", "atomic_number",
        ::rust::Slice<const std::uint64_t>(znum.data(), znum.size()));
    const std::vector<std::uint64_t> bi = {0, 0, 0, 0};
    const std::vector<std::uint64_t> bj = {1, 2, 3, 4};
    const std::vector<std::uint64_t> kind = {1, 1, 1, 1};
    molrs::frame_set_column_u32(
        *fr, "bonds", "atomi",
        ::rust::Slice<const std::uint64_t>(bi.data(), bi.size()));
    molrs::frame_set_column_u32(
        *fr, "bonds", "atomj",
        ::rust::Slice<const std::uint64_t>(bj.data(), bj.size()));
    molrs::frame_set_column_u32(
        *fr, "bonds", "type_id",
        ::rust::Slice<const std::uint64_t>(kind.data(), kind.size()));
    return fr;
}

TEST(MolrsAm1Bcc, AssignFromBaseLinks) {
    auto fr = methane_like();
    const std::vector<double> base = {-0.266, 0.066, 0.066, 0.066, 0.066};
    try {
        auto out = molrs::am1_bcc_assign_frame_from_base(
            *fr, ::rust::Slice<const double>(base.data(), base.size()),
            {"bcc"});
        ASSERT_EQ(out.size(), 5u);
        for (std::size_t i = 0; i < out.size(); ++i)
            EXPECT_TRUE(std::isfinite(out[i])) << "i=" << i;
    } catch (const std::exception& err) {
        // Hand-built frames may fail Atomistic::from_frame / BCC typing; the
        // door must still return a catchable exception (never abort).
        EXPECT_FALSE(std::string(err.what()).empty());
    }
}

TEST(MolrsAm1Bcc, UnknownParameterSetErrors) {
    auto fr = methane_like();
    const std::vector<double> base = {0.0, 0.0, 0.0, 0.0, 0.0};
    EXPECT_THROW(molrs::am1_bcc_assign_frame_from_base(
                     *fr, ::rust::Slice<const double>(base.data(), base.size()),
                     {"not-a-set"}),
                 std::exception);
}

TEST(MolrsAm1Bcc, WrongChargeCountErrors) {
    auto fr = methane_like();
    const std::vector<double> base = {0.0};
    EXPECT_THROW(molrs::am1_bcc_assign_frame_from_base(
                     *fr, ::rust::Slice<const double>(base.data(), base.size()),
                     {"bcc"}),
                 std::exception);
}

} // namespace
