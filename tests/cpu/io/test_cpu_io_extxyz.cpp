/// @file tests/cpu/io/test_cpu_io_extxyz.cpp
/// @brief Tests for cpu::load_extxyz -- (ext)XYZ loading through molrs.
///
/// Covers gpw-edl-aimd acceptance criteria:
///   - ac-001: load the real 1214-atom final_frame.xyz; correct atom count,
///             box from the Lattice= header, and per-atom type = Z for the
///             five elements {H, O, Na, Cl, Zn}. Plus a small deterministic
///             inline fixture (no external dependency).
///   - ac-002: an ExtXYZ file lacking a Lattice= header throws a named
///             exception -- no fallback box.

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/io/extxyz.h"

using ATV::cpu::Frame;
using ATV::cpu::load_extxyz;

namespace {

// Write `content` to a unique temp file and return its path.
std::string write_temp(const std::string& tag, const std::string& content) {
    const std::string path =
        (std::filesystem::temp_directory_path() / ("atv_extxyz_" + tag + ".xyz"))
            .string();
    std::ofstream out(path);
    out << content;
    out.close();
    return path;
}

std::map<int, int> type_histogram(const Frame& f) {
    std::map<int, int> hist;
    for (int i = 0; i < f.atoms.size(); ++i) hist[f.atoms.type[i]]++;
    return hist;
}

} // namespace

// ── ac-001: small deterministic inline fixture ─────────────────────────────
TEST(CpuIoExtxyz, InlineFixtureLoadsAtomsBoxAndZ) {
    // 3 atoms (O, H, H) in an orthorhombic box via a Lattice= header.
    const std::string xyz =
        "3\n"
        "Lattice=\"5.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 7.0\" "
        "Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n"
        "O 0.0 0.0 0.0\n"
        "H 0.9 0.0 0.0\n"
        "H -0.3 0.85 0.0\n";
    const std::string path = write_temp("inline", xyz);

    Frame f = load_extxyz(path);
    std::remove(path.c_str());

    ASSERT_EQ(f.atoms.size(), 3);
    EXPECT_EQ(f.atoms.type[0], 8); // O
    EXPECT_EQ(f.atoms.type[1], 1); // H
    EXPECT_EQ(f.atoms.type[2], 1); // H

    EXPECT_NEAR(f.box.length(0), 5.0, 1e-6);
    EXPECT_NEAR(f.box.length(1), 6.0, 1e-6);
    EXPECT_NEAR(f.box.length(2), 7.0, 1e-6);

    EXPECT_NEAR(f.atoms.x[1], 0.9, 1e-9);
    EXPECT_NEAR(f.atoms.y[2], 0.85, 1e-9);
}

// ── ac-002: missing Lattice header throws, no fallback box ─────────────────
TEST(CpuIoExtxyz, MissingLatticeThrows) {
    const std::string xyz =
        "2\n"
        "Properties=species:S:1:pos:R:3\n"
        "O 0.0 0.0 0.0\n"
        "H 0.9 0.0 0.0\n";
    const std::string path = write_temp("nolattice", xyz);

    EXPECT_THROW(load_extxyz(path), std::runtime_error);
    std::remove(path.c_str());
}

// ── ac-001: the real 1214-atom ZnO electric-double-layer frame ─────────────
#ifdef ATV_FINAL_FRAME_XYZ
TEST(CpuIoExtxyz, FinalFrameZnoEdlLoadsFully) {
    Frame f = load_extxyz(ATV_FINAL_FRAME_XYZ);

    ASSERT_EQ(f.atoms.size(), 1214);

    // Box from Lattice="13.128 0 0  0 10.582 0  0 0 47.05377225".
    EXPECT_NEAR(f.box.length(0), 13.128, 1e-6);
    EXPECT_NEAR(f.box.length(1), 10.582, 1e-6);
    EXPECT_NEAR(f.box.length(2), 47.05377225, 1e-6);

    // Per-element Z counts: H=540, O=462, Na=10, Cl=10, Zn=192.
    auto hist = type_histogram(f);
    EXPECT_EQ(hist.size(), 5u);
    EXPECT_EQ(hist[1], 540);  // H
    EXPECT_EQ(hist[8], 462);  // O
    EXPECT_EQ(hist[11], 10);  // Na
    EXPECT_EQ(hist[17], 10);  // Cl
    EXPECT_EQ(hist[30], 192); // Zn
}
#endif
