/// @file tests/cpu/test_atom_variable_n.cpp
/// @brief Unit tests for the variable-N Atom primitives GCMC lifts out of
///        comm/: `Atom::copy_row` (the promoted migration row-move) and
///        `Atom::acquire_tag` (the global next-free-tag allocator).

#include <gtest/gtest.h>

#include <random>
#include <set>

#include "atomiverse/cpu/frame/atom.h"
#include "atomiverse/cpu/frame/frame.h"

using namespace ATV;

namespace {

/// @brief Fill row `i` with distinct, recognizable values keyed by `seed`.
void fill_row(cpu::Atom& a, int i, double seed) {
    a.x[i] = seed + 0.1;
    a.y[i] = seed + 0.2;
    a.z[i] = seed + 0.3;
    a.vx[i] = seed + 0.4;
    a.vy[i] = seed + 0.5;
    a.vz[i] = seed + 0.6;
    a.fx[i] = seed + 0.7;
    a.fy[i] = seed + 0.8;
    a.fz[i] = seed + 0.9;
    a.mass[i] = seed + 1.0;
    a.charge[i] = seed + 1.1;
    a.type[i] = static_cast<int>(seed);
    a.tag[i] = static_cast<TagInt>(seed) + 100;
    a.image(i, 0) = 1;
    a.image(i, 1) = -1;
    a.image(i, 2) = 2;
}

} // namespace

TEST(AtomCopyRow, CopiesEverySoAField) {
    cpu::Frame frame(4);
    auto& a = frame.atoms;
    fill_row(a, 3, 30.0);
    fill_row(a, 1, 10.0);

    a.copy_row(1, 3);

    EXPECT_EQ(a.x[1], a.x[3]);
    EXPECT_EQ(a.y[1], a.y[3]);
    EXPECT_EQ(a.z[1], a.z[3]);
    EXPECT_EQ(a.vx[1], a.vx[3]);
    EXPECT_EQ(a.vy[1], a.vy[3]);
    EXPECT_EQ(a.vz[1], a.vz[3]);
    EXPECT_EQ(a.fx[1], a.fx[3]);
    EXPECT_EQ(a.fy[1], a.fy[3]);
    EXPECT_EQ(a.fz[1], a.fz[3]);
    EXPECT_EQ(a.mass[1], a.mass[3]);
    EXPECT_EQ(a.charge[1], a.charge[3]);
    EXPECT_EQ(a.type[1], a.type[3]);
    EXPECT_EQ(a.tag[1], a.tag[3]);
    for (int d = 0; d < 3; ++d)
        EXPECT_EQ(a.image(1, d), a.image(3, d));
}

TEST(AtomCopyRow, TopologyFollowsTheRow) {
    cpu::Frame frame(4);
    auto& a = frame.atoms;
    // Tags are identity by construction; bond 3-0 is owned by atom 0
    // (smaller tag), angle 0-3-2 by its center atom 3.
    a.map_clear();
    a.map_set(0, 4);
    a.add_bond(/*ti=*/3, /*tj=*/0, /*type=*/7);
    a.add_angle(/*ti=*/0, /*tj=*/3, /*tk=*/2, /*type=*/5);

    a.copy_row(1, 3);

    EXPECT_EQ(a.num_angle(1), a.num_angle(3));
    EXPECT_EQ(a.angle_atom(1, 0), a.angle_atom(3, 0));
    EXPECT_EQ(a.angle_atom(1, 1), a.angle_atom(3, 1));
    EXPECT_EQ(a.angle_type(1, 0), a.angle_type(3, 0));

    a.copy_row(2, 0); // bond owner row moves too
    EXPECT_EQ(a.num_bond(2), 1);
    EXPECT_EQ(a.bond_atom(2, 0), a.bond_atom(0, 0));
    EXPECT_EQ(a.bond_type(2, 0), 7);
}

TEST(AtomCopyRow, SameRowIsANoOp) {
    cpu::Frame frame(2);
    auto& a = frame.atoms;
    fill_row(a, 0, 5.0);
    a.copy_row(0, 0);
    EXPECT_EQ(a.x[0], 5.1);
    EXPECT_EQ(a.tag[0], 105);
}

TEST(AtomGrow, ChunkedGrowPreservesEveryLiveRow) {
    // Regression pin for a real data-loss bug: grow() trusted
    // xt::xtensor::resize to keep the prefix, and it does NOT — the first
    // GCMC insertion past the initial capacity (a chunked 64 → 4160 grow)
    // zeroed every live atom, collapsed them onto the origin, and froze the
    // Markov chain on a NaN energy. grow() now reallocates through
    // set_capacity's explicit copy.
    cpu::Frame frame(45); // ctor floors capacity at 64
    auto& a = frame.atoms;
    for (int i = 0; i < 45; ++i) {
        fill_row(a, i, 10.0 * i);
        a.tag[i] = 1000 + i;
    }
    const int nmax0 = a.nmax;
    a.grow(a.nmax + 1); // forces the chunked reallocation
    ASSERT_GT(a.nmax, nmax0);
    for (int i = 0; i < 45; ++i) {
        EXPECT_EQ(a.x[i], 10.0 * i + 0.1) << "row " << i;
        EXPECT_EQ(a.vz[i], 10.0 * i + 0.6);
        EXPECT_EQ(a.mass[i], 10.0 * i + 1.0);
        EXPECT_EQ(a.tag[i], 1000 + i);
    }
    // The new tail is zeroed, tags identity-seeded (the grow contract).
    EXPECT_EQ(a.x[50], 0.0);
    EXPECT_EQ(a.tag[50], 50);
}

TEST(AtomAcquireTag, SeedsPastTheLiveMaximum) {
    cpu::Frame frame(3); // identity tags 0,1,2
    frame.atoms.tag[1] = 41; // a restart-style sparse tag
    EXPECT_EQ(frame.atoms.acquire_tag(), 42);
    EXPECT_EQ(frame.atoms.acquire_tag(), 43); // monotonic thereafter
}

TEST(AtomAcquireTag, UniqueAcrossRandomInsertDeleteCycles) {
    // ac-004: 1000 random insert/delete cycles never produce a duplicate
    // tag. Delete = compact-last-into-victim (the GCMC deletion primitive);
    // insert = write a fresh row at nlocal with acquire_tag().
    cpu::Frame frame(8);
    auto& a = frame.atoms;
    std::mt19937_64 rng(2026u);

    for (int cycle = 0; cycle < 1000; ++cycle) {
        std::uniform_int_distribution<int> coin(0, 1);
        if (a.nlocal > 1 && coin(rng)) {
            std::uniform_int_distribution<int> pick(0, a.nlocal - 1);
            const int victim = pick(rng);
            a.copy_row(victim, a.nlocal - 1);
            --a.nlocal;
        } else {
            if (a.nlocal + 1 > a.nmax)
                a.grow(a.nlocal + 1);
            a.tag[a.nlocal] = a.acquire_tag();
            ++a.nlocal;
        }
        std::set<TagInt> seen;
        for (int i = 0; i < a.nlocal; ++i)
            EXPECT_TRUE(seen.insert(a.tag[i]).second)
                << "duplicate tag " << a.tag[i] << " at cycle " << cycle;
    }
}
