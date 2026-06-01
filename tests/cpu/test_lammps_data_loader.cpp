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
#include <cstdint>
#include <string>

namespace {

#ifndef ATV_SOURCE_DIR
#define ATV_SOURCE_DIR "."
#endif

const std::string kDataPath =
    std::string(ATV_SOURCE_DIR)
    + "/benchmarks/md/reference/polyethylene/msu_pe.data";

TEST(LammpsDataLoader, MsuPolyethyleneCountsAndBox) {
    ATV::cpu::Frame frame = ATV::cpu::load_lammps_data(kDataPath);

    EXPECT_EQ(frame.atoms.size(), 10000);

    // Topology lives in the LAMMPS per-atom layout (num_bond[i] /
    // bond_atom[i][m]). Each bond/angle/dihedral is stored exactly once on
    // its Newton owner (smaller-tag atom), so summing the per-atom counts
    // reproduces the deck header totals.
    long total_bonds = 0, total_angles = 0, total_dihedrals = 0;
    for (int i = 0; i < frame.atoms.size(); ++i) {
        total_bonds += frame.atoms.num_bond[i];
        total_angles += frame.atoms.num_angle[i];
        total_dihedrals += frame.atoms.num_dihedral[i];
    }
    EXPECT_EQ(total_bonds, 9990);
    EXPECT_EQ(total_angles, 9980);
    EXPECT_EQ(total_dihedrals, 9970);

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
    const int n = frame.atoms.size();

    // LAMMPS file is 1-indexed; the loader stores 0-indexed rows. The first
    // bond row "1 1 1 2" makes atom 0 the Newton owner with partner row 1.
    ASSERT_GE(frame.atoms.num_bond[0], 1);
    EXPECT_EQ(frame.atoms.bond_atom(0, 0), 1);

    // Every stored bond partner must resolve to a valid local row [0, n).
    // (The owner row is the loop index, in range by construction.)
    for (int i = 0; i < n; ++i) {
        for (int m = 0; m < frame.atoms.num_bond[i]; ++m) {
            const std::int64_t partner = frame.atoms.bond_atom(i, m);
            EXPECT_GE(partner, 0);
            EXPECT_LT(partner, n);
        }
    }
}

} // namespace
