/// @file tests/cpu/gpw/test_cpu_gpw_cp2k_io.cpp
/// @brief Acceptance ac-001 (native-xc / gpw spec): GTH potential parser parses
///        CP2K data files exactly for Ar, H, O.
///
/// Reference values are read by hand from the CP2K data files shipped with the
/// repo at ATV_CP2K_DATA_DIR (set by CMake). H/O come from GTH_POTENTIALS
/// (GTH-PBE), Ar from POTENTIAL (GTH-PADE, used by tests/QS/regtest-gpw-1/Ar.inp).

#include <string>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/internal/cp2k_io.h"

using ATV::cpu::gpw::GthPotential;
using ATV::cpu::gpw::parse_gth_potential;

namespace {
const std::string kDataDir = ATV_CP2K_DATA_DIR; // compile-time absolute path
const std::string kGthPot = kDataDir + "/GTH_POTENTIALS";
const std::string kPotential = kDataDir + "/POTENTIAL";
} // namespace

// H GTH-PBE-q1: 1 valence electron, local-only (no nonlocal projectors).
TEST(GthPotentialParser, HydrogenPbeQ1) {
    GthPotential p = parse_gth_potential(kGthPot, "H", "GTH-PBE-q1");
    EXPECT_EQ(p.zion, 1);
    EXPECT_DOUBLE_EQ(p.rloc, 0.20000000);
    ASSERT_EQ(p.cloc.size(), 2u);
    EXPECT_DOUBLE_EQ(p.cloc[0], -4.17890044);
    EXPECT_DOUBLE_EQ(p.cloc[1], 0.72446331);
    EXPECT_TRUE(p.channels.empty());
}

// O GTH-PBE-q6: 6 valence electrons, two channels (l=0 with one projector,
// l=1 declared with zero projectors).
TEST(GthPotentialParser, OxygenPbeQ6) {
    GthPotential p = parse_gth_potential(kGthPot, "O", "GTH-PBE-q6");
    EXPECT_EQ(p.zion, 6);
    EXPECT_DOUBLE_EQ(p.rloc, 0.24455430);
    ASSERT_EQ(p.cloc.size(), 2u);
    EXPECT_DOUBLE_EQ(p.cloc[0], -16.66721480);
    EXPECT_DOUBLE_EQ(p.cloc[1], 2.48731132);
    ASSERT_EQ(p.channels.size(), 2u);

    EXPECT_EQ(p.channels[0].l, 0);
    EXPECT_DOUBLE_EQ(p.channels[0].r, 0.22095592);
    ASSERT_EQ(p.channels[0].nproj, 1);
    EXPECT_DOUBLE_EQ(p.channels[0].hij(0, 0), 18.33745811);

    EXPECT_EQ(p.channels[1].l, 1);
    EXPECT_DOUBLE_EQ(p.channels[1].r, 0.21133247);
    EXPECT_EQ(p.channels[1].nproj, 0);
}

// Ar GTH-PADE-q8 (from POTENTIAL): 8 valence electrons, l=0 has a 2x2 symmetric
// h matrix spread over two lines, l=1 has a single projector.
TEST(GthPotentialParser, ArgonPadeQ8) {
    GthPotential p = parse_gth_potential(kPotential, "Ar", "GTH-PADE-q8");
    EXPECT_EQ(p.zion, 8);
    EXPECT_DOUBLE_EQ(p.rloc, 0.40000000);
    ASSERT_EQ(p.cloc.size(), 1u);
    EXPECT_DOUBLE_EQ(p.cloc[0], -7.10000000);
    ASSERT_EQ(p.channels.size(), 2u);

    const auto& c0 = p.channels[0];
    EXPECT_DOUBLE_EQ(c0.r, 0.31738081);
    ASSERT_EQ(c0.nproj, 2);
    EXPECT_DOUBLE_EQ(c0.hij(0, 0), 10.24948699);
    EXPECT_DOUBLE_EQ(c0.hij(0, 1), -2.16984522);
    EXPECT_DOUBLE_EQ(c0.hij(1, 0), -2.16984522); // symmetric fill
    EXPECT_DOUBLE_EQ(c0.hij(1, 1), 5.60251627);

    const auto& c1 = p.channels[1];
    EXPECT_DOUBLE_EQ(c1.r, 0.35161921);
    ASSERT_EQ(c1.nproj, 1);
    EXPECT_DOUBLE_EQ(c1.hij(0, 0), 4.97880101);
}

// Missing entries must throw, not silently return garbage.
TEST(GthPotentialParser, MissingEntryThrows) {
    EXPECT_THROW(parse_gth_potential(kGthPot, "H", "GTH-NOPE-q9"), std::runtime_error);
}
