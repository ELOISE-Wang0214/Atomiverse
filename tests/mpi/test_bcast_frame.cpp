/// @file tests/mpi/test_bcast_frame.cpp
/// @brief cpu::comm::bcast_frame roundtrip tests
///        (mpi-launcher-embedded-python-01).
///
/// Launched under `mpirun -n N`. Runs at N=1, 2, 4 per CTest config.
/// bcast_frame is a collective on MPI_COMM_WORLD: the root rank packs every
/// owned atom row via comm::pack_exchange, broadcasts a POD header
/// {n_atoms, payload_doubles, has_bonded, cpu::Box, SystemScalars} plus the
/// payload, and
/// every non-root rank reconstructs a full replica (all rows owned,
/// nghost = 0, tag→local map rebuilt, box + scalars overwritten). The wire is
/// pure byte transit — every assertion below is EXACT (==, zero tolerance).
///
/// Verifies:
///   - Multi-rank roundtrip: every schema field, per-atom bonded topology,
///     box, and SystemScalars replicate exactly; replica invariants
///     (nlocal == n, nghost == 0, has_bonded(), map usable); root untouched.
///   - World size 1: strict no-op — the frame is field-for-field unchanged.
///   - 0-atom frame: header-only broadcast completes; box/scalars still land.
///   - No-topology frame: has_bonded() stays false on the replica.
///   - Exchanged extension columns: both-sides-registered columns ride the
///     payload exactly; a non-root-only extra exchanged column makes the
///     non-root decode consume more doubles than the root announced →
///     std::runtime_error on that rank (collectives complete before decode,
///     so all ranks stay in lockstep — no deadlock).
///
/// Not testable here: the MPI-uninitialized degrade path (this binary's main
/// initializes MPI before any test runs); it is covered by the non-MPI-build
/// compile path.

#include <gtest/gtest.h>

#include <mpi.h>

#include <cstdint>
#include <stdexcept>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/comm/bcast_frame.h"
#include "atomiverse/cpu/frame/frame.h"

using namespace ATV;

namespace {

int world_size() {
    int s = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &s);
    return s;
}

int world_rank() {
    int r = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    return r;
}

/// Noncontiguous, nonmonotonic global tags for the 8-atom reference frame.
constexpr TagInt kRefTags[8] = {7, 101, 13, 55, 2, 89, 34, 68};

/// @brief Deterministic 8-atom reference frame: every core field distinct and
///        nonzero, bonded topology (bond/angle/dihedral/special), an
///        orthorhombic non-cubic box, and nonzero SystemScalars. Constructed
///        from literals only, so two independent constructions are bitwise
///        identical — each rank builds its own expected copy locally.
cpu::Frame make_reference_frame() {
    constexpr int n = 8;
    cpu::Frame f(n);
    for (int i = 0; i < n; ++i) {
        f.atoms.x[i] = StorageReal(0.5 + 1.25 * i);
        f.atoms.y[i] = StorageReal(-1.0 + 0.5 * i);
        f.atoms.z[i] = StorageReal(2.0 + 0.75 * i);
        f.atoms.vx[i] = StorageReal(0.01 * (i + 1));
        f.atoms.vy[i] = StorageReal(-0.02 * (i + 1));
        f.atoms.vz[i] = StorageReal(0.03 * (i + 1));
        f.atoms.fx[i] = AccReal(1.5 * i - 3.0);
        f.atoms.fy[i] = AccReal(0.25 * i + 0.125);
        f.atoms.fz[i] = AccReal(-0.5 * i + 1.0);
        f.atoms.mass[i] = StorageReal(12.011 + i);
        f.atoms.charge[i] =
            StorageReal((i % 2 == 0 ? 0.42 : -0.42) + 0.01 * i); // nonzero
        f.atoms.type[i] = 1 + (i % 3);
        f.atoms.tag[i] = kRefTags[i]; // noncontiguous
        f.atoms.image(i, 0) = static_cast<std::int8_t>(i - 4);
        f.atoms.image(i, 1) = static_cast<std::int8_t>(2 * i - 7);
        f.atoms.image(i, 2) = static_cast<std::int8_t>(-i);
    }
    // add_* prerequisites: nlocal set (Frame ctor), tags set, map usable.
    f.atoms.map_clear();
    f.atoms.map_set(0, n);
    // Chain over rows 0-1-2-3 by global tag: 7 - 101 - 13 - 55.
    f.atoms.add_bond(kRefTags[0], kRefTags[1], 1);
    f.atoms.add_bond(kRefTags[1], kRefTags[2], 2);
    f.atoms.add_bond(kRefTags[2], kRefTags[3], 3);
    f.atoms.add_angle(kRefTags[0], kRefTags[1], kRefTags[2], 1);
    f.atoms.add_angle(kRefTags[1], kRefTags[2], kRefTags[3], 2);
    f.atoms.add_dihedral(kRefTags[0], kRefTags[1], kRefTags[2], kRefTags[3],
                         1);
    f.atoms.build_special_list();

    f.box = Mat3::diagonal(12.0, 15.5, 19.25); // orthorhombic, non-cubic

    f.scalars.energy = -123.5;
    f.scalars.energy_pair = -100.25;
    f.scalars.energy_bond = -13.5;
    f.scalars.energy_angle = -6.25;
    f.scalars.energy_dihedral = -3.5;
    f.scalars.kinetic = 45.75;
    f.scalars.conserved = -77.75;
    for (int k = 0; k < 6; ++k) {
        f.scalars.stress[k] = 0.5 * (k + 1);
        f.scalars.stress_kinetic[k] = -0.25 * (k + 1);
    }
    return f;
}

/// @brief Assert `got.box` == `want.box` element-for-element (exact).
void expect_box_equal(const cpu::Frame& got, const cpu::Frame& want) {
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            EXPECT_EQ(got.box(r, c), want.box(r, c))
                << "box(" << r << ", " << c << ")";
    EXPECT_EQ(got.box.is_orthorhombic(), want.box.is_orthorhombic());
}

/// @brief Assert every SystemScalars field equal (exact).
void expect_scalars_equal(const cpu::Frame& got, const cpu::Frame& want) {
    EXPECT_EQ(got.scalars.energy, want.scalars.energy);
    EXPECT_EQ(got.scalars.energy_pair, want.scalars.energy_pair);
    EXPECT_EQ(got.scalars.energy_bond, want.scalars.energy_bond);
    EXPECT_EQ(got.scalars.energy_angle, want.scalars.energy_angle);
    EXPECT_EQ(got.scalars.energy_dihedral, want.scalars.energy_dihedral);
    EXPECT_EQ(got.scalars.kinetic, want.scalars.kinetic);
    EXPECT_EQ(got.scalars.conserved, want.scalars.conserved);
    for (int k = 0; k < 6; ++k) {
        EXPECT_EQ(got.scalars.stress[k], want.scalars.stress[k]) << "k=" << k;
        EXPECT_EQ(got.scalars.stress_kinetic[k],
                  want.scalars.stress_kinetic[k])
            << "k=" << k;
    }
}

/// @brief Assert every per-atom field of rows [0, want.atoms.nlocal) equal
///        (exact ==; the payload is byte transit — no tolerance is legal),
///        plus the replica invariants (counts, has_bonded, tag→local map).
///        Topology strides may legitimately differ (the receiver grows
///        defensively), so wide arrays are compared through the per-row
///        counts, never by full stride.
void expect_frame_equal(const cpu::Frame& got, const cpu::Frame& want) {
    const int n = want.atoms.nlocal;
    ASSERT_EQ(got.atoms.nlocal, n);
    EXPECT_EQ(got.atoms.nghost, 0);
    EXPECT_EQ(got.atoms.has_bonded(), want.atoms.has_bonded());

    for (int i = 0; i < n; ++i) {
        SCOPED_TRACE("atom row " + std::to_string(i));
        EXPECT_EQ(got.atoms.x[i], want.atoms.x[i]);
        EXPECT_EQ(got.atoms.y[i], want.atoms.y[i]);
        EXPECT_EQ(got.atoms.z[i], want.atoms.z[i]);
        EXPECT_EQ(got.atoms.vx[i], want.atoms.vx[i]);
        EXPECT_EQ(got.atoms.vy[i], want.atoms.vy[i]);
        EXPECT_EQ(got.atoms.vz[i], want.atoms.vz[i]);
        EXPECT_EQ(got.atoms.fx[i], want.atoms.fx[i]);
        EXPECT_EQ(got.atoms.fy[i], want.atoms.fy[i]);
        EXPECT_EQ(got.atoms.fz[i], want.atoms.fz[i]);
        EXPECT_EQ(got.atoms.mass[i], want.atoms.mass[i]);
        EXPECT_EQ(got.atoms.charge[i], want.atoms.charge[i]);
        EXPECT_EQ(got.atoms.type[i], want.atoms.type[i]);
        EXPECT_EQ(got.atoms.tag[i], want.atoms.tag[i]);
        for (int d = 0; d < 3; ++d)
            EXPECT_EQ(got.atoms.image(i, d), want.atoms.image(i, d))
                << "image axis " << d;

        // tag → local map is rebuilt and usable on the replica.
        ASSERT_EQ(got.atoms.map.count(want.atoms.tag[i]), 1u);
        EXPECT_EQ(got.atoms.map.at(want.atoms.tag[i]), i);

        // Per-atom bonded topology, through the counts (stride-agnostic).
        ASSERT_EQ(got.atoms.num_bond[i], want.atoms.num_bond[i]);
        for (int m = 0; m < want.atoms.num_bond[i]; ++m) {
            EXPECT_EQ(got.atoms.bond_atom(i, m), want.atoms.bond_atom(i, m));
            EXPECT_EQ(got.atoms.bond_type(i, m), want.atoms.bond_type(i, m));
        }
        ASSERT_EQ(got.atoms.num_angle[i], want.atoms.num_angle[i]);
        for (int m = 0; m < want.atoms.num_angle[i]; ++m) {
            EXPECT_EQ(got.atoms.angle_atom(i, 2 * m + 0),
                      want.atoms.angle_atom(i, 2 * m + 0));
            EXPECT_EQ(got.atoms.angle_atom(i, 2 * m + 1),
                      want.atoms.angle_atom(i, 2 * m + 1));
            EXPECT_EQ(got.atoms.angle_type(i, m), want.atoms.angle_type(i, m));
        }
        ASSERT_EQ(got.atoms.num_dihedral[i], want.atoms.num_dihedral[i]);
        for (int m = 0; m < want.atoms.num_dihedral[i]; ++m) {
            for (int k = 0; k < 3; ++k)
                EXPECT_EQ(got.atoms.dihedral_atom(i, 3 * m + k),
                          want.atoms.dihedral_atom(i, 3 * m + k));
            EXPECT_EQ(got.atoms.dihedral_type(i, m),
                      want.atoms.dihedral_type(i, m));
        }
        for (int k = 0; k < 3; ++k)
            ASSERT_EQ(got.atoms.nspecial(i, k), want.atoms.nspecial(i, k))
                << "nspecial tier " << k;
        for (int s = 0; s < want.atoms.nspecial(i, 2); ++s)
            EXPECT_EQ(got.atoms.special(i, s), want.atoms.special(i, s))
                << "special slot " << s;
    }

    expect_box_equal(got, want);
    expect_scalars_equal(got, want);
}

} // namespace

// ── (a) Multi-rank roundtrip: every field exact, invariants, root untouched ─
//
// Root builds the reference frame; non-root ranks start from Frame(0). After
// the broadcast EVERY rank (root included — which proves the root frame is
// untouched) must equal an independently constructed reference bit-for-bit.
// Also runs at N=1, where it degenerates to the root-untouched assertion.

TEST(MpiBcastFrame, RoundtripReplicatesEveryFieldExactly) {
    const cpu::Frame want = make_reference_frame();
    cpu::Frame frame = (world_rank() == 0) ? make_reference_frame()
                                           : cpu::Frame(0);

    cpu::comm::bcast_frame(frame, /*root=*/0);

    EXPECT_TRUE(frame.atoms.has_bonded());
    expect_frame_equal(frame, want);
}

// ── (b) World size 1: strict no-op degrade ──────────────────────────────────
//
// Only meaningful at N=1 (at N>1 non-root frames are legitimately rewritten);
// skipping is world-uniform, so lockstep is preserved.

TEST(MpiBcastFrame, SerialWorldDegradesToNoOp) {
    if (world_size() != 1)
        GTEST_SKIP() << "serial-degrade contract is a world-size-1 property";

    const cpu::Frame want = make_reference_frame();
    cpu::Frame frame = make_reference_frame();

    cpu::comm::bcast_frame(frame, /*root=*/0);

    expect_frame_equal(frame, want); // field-for-field identical to pre-call
}

// ── (c) 0-atom frame: header-only broadcast ─────────────────────────────────
//
// No payload rows; the header (box + scalars) must still land on non-root.

TEST(MpiBcastFrame, ZeroAtomFrameBroadcastsHeaderOnly) {
    cpu::Frame frame(0);
    if (world_rank() == 0) {
        frame.box = Mat3::diagonal(6.0, 7.5, 9.0);
        frame.scalars.energy = -42.125;
        frame.scalars.kinetic = 3.5;
    }

    EXPECT_NO_THROW(cpu::comm::bcast_frame(frame, /*root=*/0));

    EXPECT_EQ(frame.atoms.nlocal, 0);
    EXPECT_EQ(frame.atoms.nghost, 0);
    EXPECT_EQ(frame.box(0, 0), 6.0);
    EXPECT_EQ(frame.box(1, 1), 7.5);
    EXPECT_EQ(frame.box(2, 2), 9.0);
    EXPECT_EQ(frame.scalars.energy, -42.125);
    EXPECT_EQ(frame.scalars.kinetic, 3.5);
}

// ── (d) No-topology frame: replica stays atom_style atomic ─────────────────

TEST(MpiBcastFrame, NoTopologyFrameStaysAtomic) {
    constexpr int n = 4;
    cpu::Frame frame = (world_rank() == 0) ? cpu::Frame(n) : cpu::Frame(0);
    if (world_rank() == 0) {
        for (int i = 0; i < n; ++i) {
            frame.atoms.x[i] = StorageReal(1.0 + i);
            frame.atoms.y[i] = StorageReal(2.0 + i);
            frame.atoms.z[i] = StorageReal(3.0 + i);
            frame.atoms.mass[i] = StorageReal(4.0);
            frame.atoms.tag[i] = TagInt(20 + 3 * i);
        }
        frame.box = Mat3::cubic(10.0);
    }

    cpu::comm::bcast_frame(frame, /*root=*/0);

    EXPECT_EQ(frame.atoms.nlocal, n);
    EXPECT_FALSE(frame.atoms.has_bonded());
    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(frame.atoms.x[i], StorageReal(1.0 + i));
        EXPECT_EQ(frame.atoms.tag[i], TagInt(20 + 3 * i));
        EXPECT_EQ(frame.atoms.num_bond[i], 0);
    }
}

// ── (e-i) Exchanged extension column rides the roundtrip exactly ────────────
//
// Precondition of bcast_frame: exchanged extension-column registration must
// match on every rank. Here it does — the column's values must transit
// bit-exactly alongside the schema fields.

TEST(MpiBcastFrame, ExchangedExtensionColumnRidesTheRoundtrip) {
    constexpr int n = 3;
    cpu::Frame frame = (world_rank() == 0) ? cpu::Frame(n) : cpu::Frame(0);
    frame.insert("aux", ColDType::Storage, 2, /*migrates=*/true,
                 /*exchanged=*/true); // every rank, identical registration
    if (world_rank() == 0) {
        auto aux = frame.span(ColKey<StorageReal>{"atoms", "aux"});
        for (int i = 0; i < n; ++i) {
            frame.atoms.x[i] = StorageReal(0.5 * i);
            frame.atoms.tag[i] = TagInt(i + 1);
            aux(i, 0) = StorageReal(10.0 + i);
            aux(i, 1) = StorageReal(-20.0 - i);
        }
        frame.box = Mat3::cubic(8.0);
    }

    cpu::comm::bcast_frame(frame, /*root=*/0);

    ASSERT_EQ(frame.atoms.nlocal, n);
    auto aux = frame.span(ColKey<StorageReal>{"atoms", "aux"});
    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(frame.atoms.x[i], StorageReal(0.5 * i));
        EXPECT_EQ(aux(i, 0), StorageReal(10.0 + i)) << "row " << i;
        EXPECT_EQ(aux(i, 1), StorageReal(-20.0 - i)) << "row " << i;
    }
}

// ── (e-ii) Mismatched exchanged registration throws on non-root ─────────────
//
// ONLY non-root ranks register an extra exchanged column, so their decode
// consumes more doubles than the root's payload_doubles announces →
// std::runtime_error on those ranks. Every rank still calls bcast_frame in
// lockstep (the collectives complete before decode, so the root and any
// passing rank return normally while the mismatched ranks throw — no
// deadlock). At N=1 there is no non-root rank and the call is a no-op.

TEST(MpiBcastFrame, MismatchedExchangedRegistrationThrowsOnNonRoot) {
    constexpr int n = 4;
    cpu::Frame frame = (world_rank() == 0) ? cpu::Frame(n) : cpu::Frame(0);
    if (world_rank() == 0) {
        for (int i = 0; i < n; ++i) {
            frame.atoms.x[i] = StorageReal(i);
            frame.atoms.tag[i] = TagInt(i);
        }
        frame.box = Mat3::cubic(5.0);
    } else {
        // The registration divergence under test: non-root only.
        frame.insert("rogue", ColDType::Storage, 1, /*migrates=*/true,
                     /*exchanged=*/true);
    }

    if (world_rank() == 0) {
        EXPECT_NO_THROW(cpu::comm::bcast_frame(frame, /*root=*/0));
    } else {
        EXPECT_THROW(cpu::comm::bcast_frame(frame, /*root=*/0),
                     std::runtime_error);
    }
}

// ── Custom main that initializes MPI before/after GoogleTest ────────────
//
// See test_topology.cpp for the OpenMPI 5.0.x macOS teardown rationale.

int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    ::testing::InitGoogleTest(&argc, argv);
    const int rc = RUN_ALL_TESTS();
    MPI_Finalize();
    return rc;
}
