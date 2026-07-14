/// @file tests/cpu/mc/test_transaction_rollback.cpp
/// @brief ac-002 of mc-gcmc-npt: `mc::Transaction` is a Frame patch, and
///        rollback restores the frame EXACTLY (bitwise, field-by-field
///        against a pre-move copy) for every verb — a rejected displacement,
///        insertion (no orphan row), deletion (full resurrection, topology
///        included), and volume move (box AND coordinates).

#include <gtest/gtest.h>

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/mc/transaction.h"

using namespace ATV;
using Verb = cpu::mc::Transaction::Verb;

namespace {

cpu::Frame make_frame(int n) {
    cpu::Frame frame(n);
    auto& a = frame.atoms;
    for (int i = 0; i < n; ++i) {
        a.x[i] = 0.11 + i;
        a.y[i] = 0.22 + i;
        a.z[i] = 0.33 + i;
        a.vx[i] = 1.0 + 0.1 * i;
        a.vy[i] = 2.0 + 0.1 * i;
        a.vz[i] = 3.0 + 0.1 * i;
        a.fx[i] = -1.0 - i;
        a.fy[i] = -2.0 - i;
        a.fz[i] = -3.0 - i;
        a.mass[i] = 39.948;
        a.charge[i] = 0.5 - 0.25 * i;
        a.type[i] = 1 + (i % 2);
        a.image(i, 0) = static_cast<std::int8_t>(i % 3 - 1);
    }
    // A bond and an angle so a deleted row carries real topology.
    a.map_clear();
    a.map_set(0, n);
    a.add_bond(/*ti=*/0, /*tj=*/1, /*type=*/3);
    a.add_angle(/*ti=*/0, /*tj=*/1, /*tk=*/2, /*type=*/4);
    frame.box.set_cell(Mat3::diagonal(10.0, 12.0, 14.0));
    frame.scalars.energy = -1.25;
    return frame;
}

/// @brief Bitwise equality of every per-atom field and the box, rows
///        `[0, a.nlocal)`.
void expect_frames_identical(const cpu::Frame& f, const cpu::Frame& ref) {
    ASSERT_EQ(f.atoms.nlocal, ref.atoms.nlocal);
    const auto& a = f.atoms;
    const auto& b = ref.atoms;
    for (int i = 0; i < b.nlocal; ++i) {
        EXPECT_EQ(a.x[i], b.x[i]) << "row " << i;
        EXPECT_EQ(a.y[i], b.y[i]);
        EXPECT_EQ(a.z[i], b.z[i]);
        EXPECT_EQ(a.vx[i], b.vx[i]);
        EXPECT_EQ(a.vy[i], b.vy[i]);
        EXPECT_EQ(a.vz[i], b.vz[i]);
        EXPECT_EQ(a.fx[i], b.fx[i]);
        EXPECT_EQ(a.fy[i], b.fy[i]);
        EXPECT_EQ(a.fz[i], b.fz[i]);
        EXPECT_EQ(a.mass[i], b.mass[i]);
        EXPECT_EQ(a.charge[i], b.charge[i]);
        EXPECT_EQ(a.type[i], b.type[i]);
        EXPECT_EQ(a.tag[i], b.tag[i]) << "row " << i;
        for (int d = 0; d < 3; ++d)
            EXPECT_EQ(a.image(i, d), b.image(i, d));
        EXPECT_EQ(a.num_bond(i), b.num_bond(i)) << "row " << i;
        for (int m = 0; m < b.num_bond(i); ++m) {
            EXPECT_EQ(a.bond_atom(i, m), b.bond_atom(i, m));
            EXPECT_EQ(a.bond_type(i, m), b.bond_type(i, m));
        }
        EXPECT_EQ(a.num_angle(i), b.num_angle(i));
        for (int m = 0; m < b.num_angle(i); ++m) {
            EXPECT_EQ(a.angle_atom(i, 2 * m + 0), b.angle_atom(i, 2 * m + 0));
            EXPECT_EQ(a.angle_atom(i, 2 * m + 1), b.angle_atom(i, 2 * m + 1));
            EXPECT_EQ(a.angle_type(i, m), b.angle_type(i, m));
        }
        for (int d = 0; d < 3; ++d)
            EXPECT_EQ(a.nspecial(i, d), b.nspecial(i, d));
        for (int s = 0; s < b.nspecial(i, 2); ++s)
            EXPECT_EQ(a.special(i, s), b.special(i, s));
    }
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            EXPECT_EQ(f.box.cell()(r, c), ref.box.cell()(r, c));
}

} // namespace

TEST(TransactionRollback, RejectedDisplacementRestoresTheRow) {
    cpu::Frame frame = make_frame(4);
    const cpu::Frame before = frame;

    cpu::mc::Transaction t;
    t.begin(Verb::Displace, frame);
    t.record_touch(frame, 2);
    frame.atoms.x[2] += 0.7;
    frame.atoms.y[2] -= 0.4;
    frame.atoms.z[2] += 0.1;

    t.rollback(frame);
    expect_frames_identical(frame, before);
}

TEST(TransactionRollback, RejectedInsertionLeavesNoOrphanRow) {
    cpu::Frame frame = make_frame(3);
    const cpu::Frame before = frame;

    cpu::mc::Transaction t;
    t.begin(Verb::Insert, frame);
    auto& a = frame.atoms;
    if (a.nlocal + 1 > a.nmax)
        a.grow(a.nlocal + 1);
    a.x[a.nlocal] = 5.0;
    a.y[a.nlocal] = 6.0;
    a.z[a.nlocal] = 7.0;
    a.mass[a.nlocal] = 39.948;
    a.type[a.nlocal] = 1;
    a.tag[a.nlocal] = a.acquire_tag();
    ++a.nlocal;
    ASSERT_EQ(frame.atoms.nlocal, 4);

    t.rollback(frame);
    expect_frames_identical(frame, before);
}

TEST(TransactionRollback, RejectedDeletionResurrectsTheAtomIntact) {
    cpu::Frame frame = make_frame(5);
    const cpu::Frame before = frame;

    // GCMC deletion primitive: snapshot the victim, compact last into its
    // slot, shrink.
    const int victim = 1; // owns bond 0-1's partner entry via atom 0; carries
                          // angle topology as the center of 0-1-2
    cpu::mc::Transaction t;
    t.begin(Verb::Delete, frame);
    t.record_touch(frame, victim);
    frame.atoms.copy_row(victim, frame.atoms.nlocal - 1);
    --frame.atoms.nlocal;
    ASSERT_EQ(frame.atoms.nlocal, 4);

    t.rollback(frame);
    expect_frames_identical(frame, before);
}

TEST(TransactionRollback, RejectedVolumeMoveRestoresBoxAndCoordinates) {
    cpu::Frame frame = make_frame(4);
    const cpu::Frame before = frame;

    cpu::mc::Transaction t;
    t.begin(Verb::Box, frame);
    t.record_box(frame);
    const double s = 1.03;
    for (int i = 0; i < frame.atoms.nlocal; ++i) {
        t.record_touch(frame, i);
        frame.atoms.x[i] *= s;
        frame.atoms.y[i] *= s;
        frame.atoms.z[i] *= s;
    }
    frame.box.scale_isotropic(s);

    t.rollback(frame);
    // s*x/s != x in floating point — only a stored-row splice passes this.
    expect_frames_identical(frame, before);
}
