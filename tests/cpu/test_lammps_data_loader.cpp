/// @file tests/cpu/test_lammps_data_loader.cpp
/// @brief Smoke test for `cpu::load_lammps_data` against the MSU PE deck.
///
/// Validates that the molrs LAMMPS reader, exposed via the cxx bridge,
/// returns expected counts and a sane geometry. Numbers come straight
/// from `msu_pe.data`'s header:
///     10000 atoms
///      9990 bonds
///      9980 angles
///      9970 dihedrals
///     box: 0..80.0586 along each axis (orthogonal)
///     mass: 14.02 (single type)

#include <gtest/gtest.h>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/lammps_data.h"

#include <cmath>
#include <string>

namespace {

#ifndef ATV_SOURCE_DIR
#define ATV_SOURCE_DIR "."
#endif

const std::string kDataPath = std::string(ATV_SOURCE_DIR) +
    "/benchmarks/md/reference/polyethylene/msu_pe.data";

TEST(LammpsDataLoader, MsuPolyethyleneCountsAndBox) {
    ATV::cpu::Frame frame = ATV::cpu::load_lammps_data(kDataPath);

    EXPECT_EQ(frame.atoms.size(), 10000);
    EXPECT_EQ(frame.bonds.size(), 9990);
    EXPECT_EQ(frame.angles.size(), 9980);
    EXPECT_EQ(frame.dihedrals.size(), 9970);

    // Orthogonal 80.0586 Å cube.
    EXPECT_NEAR(frame.box(0, 0), 80.0586, 1e-3);
    EXPECT_NEAR(frame.box(1, 1), 80.0586, 1e-3);
    EXPECT_NEAR(frame.box(2, 2), 80.0586, 1e-3);
    EXPECT_NEAR(frame.box(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(frame.box(0, 2), 0.0, 1e-12);
    EXPECT_NEAR(frame.box(1, 2), 0.0, 1e-12);
}

TEST(LammpsDataLoader, MsuPolyethyleneAtomDetails) {
    ATV::cpu::Frame frame = ATV::cpu::load_lammps_data(kDataPath);

    // Spot-check the first atom: "1  1  1  8.6550  61.6668  5.4094"
    EXPECT_NEAR(frame.atoms.x[0], 8.6550, 1e-4);
    EXPECT_NEAR(frame.atoms.y[0], 61.6668, 1e-4);
    EXPECT_NEAR(frame.atoms.z[0], 5.4094, 1e-4);

    // Single type → 0-indexed should be 0 for every atom.
    for (int i = 0; i < frame.atoms.size(); ++i) {
        EXPECT_EQ(frame.atoms.type[i], 0);
        EXPECT_NEAR(frame.atoms.mass[i], 14.02, 1e-9);
    }
}

TEST(LammpsDataLoader, BondsAreZeroIndexed) {
    ATV::cpu::Frame frame = ATV::cpu::load_lammps_data(kDataPath);

    // LAMMPS file is 1-indexed; we should see 0-indexed pairs.
    // First bond row in msu_pe.data is "1 1 1 2", so after subtraction (i,j)=(0,1).
    EXPECT_EQ(frame.bonds.atomi[0], 0);
    EXPECT_EQ(frame.bonds.atomj[0], 1);
    // Atom indices must all fall inside [0, n_atoms).
    const int n = frame.atoms.size();
    for (int b = 0; b < frame.bonds.size(); ++b) {
        EXPECT_GE(frame.bonds.atomi[b], 0);
        EXPECT_LT(frame.bonds.atomi[b], n);
        EXPECT_GE(frame.bonds.atomj[b], 0);
        EXPECT_LT(frame.bonds.atomj[b], n);
    }
}

} // namespace
