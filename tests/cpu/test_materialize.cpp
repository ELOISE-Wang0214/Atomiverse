/// @file tests/cpu/test_materialize.cpp
/// @brief Round-trip + boundary tests for cpu::materialize / cpu::write_back.
///
/// Covers the rust-frame-lowering-02 acceptance criteria:
///   - ac-002: molrs.Frame -> cpu::Frame -> molrs.Frame round-trips
///             positions / velocities / forces / charge / topology / box /
///             system scalars element-wise.
///   - ac-003: required fields fail loud (named exception); optional fields
///             default to zero without throwing.
///   - ac-004: a cpu::Frame produced by cpu::materialize is runnable under
///             cpu::md::MD and yields a finite energy.

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/materialize.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/molrs.hpp"

using ATV::cpu::Frame;
using ATV::cpu::materialize;
using ATV::cpu::write_back;

namespace {

// Build a fully-populated molrs.Frame with three atoms, one bond, one
// angle, one dihedral, and a cubic box.
::rust::Box<molrs::FrameRef> make_full_frame() {
    auto fr = molrs::frame_new();

    const std::vector<double> x = {0.0, 1.5, 3.0};
    const std::vector<double> y = {0.1, 1.6, 3.1};
    const std::vector<double> z = {0.2, 1.7, 3.2};
    const std::vector<double> vx = {0.01, 0.02, 0.03};
    const std::vector<double> vy = {-0.01, -0.02, -0.03};
    const std::vector<double> vz = {0.04, 0.05, 0.06};
    const std::vector<double> fx = {1.0, 2.0, 3.0};
    const std::vector<double> fy = {-1.0, -2.0, -3.0};
    const std::vector<double> fz = {0.5, 0.6, 0.7};
    const std::vector<double> mass = {12.011, 1.008, 15.999};
    const std::vector<double> charge = {-0.5, 0.25, 0.25};
    const std::vector<std::int32_t> type = {1, 2, 3};
    const std::vector<std::int32_t> tag = {10, 11, 12};

    molrs::frame_set_column_f64(
        *fr, "atoms", "x", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "y", ::rust::Slice<const double>(y.data(), y.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "z", ::rust::Slice<const double>(z.data(), z.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "vx", ::rust::Slice<const double>(vx.data(), vx.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "vy", ::rust::Slice<const double>(vy.data(), vy.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "vz", ::rust::Slice<const double>(vz.data(), vz.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "fx", ::rust::Slice<const double>(fx.data(), fx.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "fy", ::rust::Slice<const double>(fy.data(), fy.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "fz", ::rust::Slice<const double>(fz.data(), fz.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "mass",
        ::rust::Slice<const double>(mass.data(), mass.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "charge",
        ::rust::Slice<const double>(charge.data(), charge.size()));
    molrs::frame_set_column_i32(
        *fr, "atoms", "type",
        ::rust::Slice<const std::int32_t>(type.data(), type.size()));
    molrs::frame_set_column_i32(
        *fr, "atoms", "tag",
        ::rust::Slice<const std::int32_t>(tag.data(), tag.size()));

    // Topology: one bond, one angle, one dihedral. molrs topology columns
    // reference atoms by tag (materialize scatters them via the tag->local
    // map), so use the atom tags {10,11,12} above -- not row indices.
    const std::vector<std::int32_t> bi = {10}, bj = {11}, bt = {1};
    molrs::frame_set_column_i32(
        *fr, "bonds", "atomi",
        ::rust::Slice<const std::int32_t>(bi.data(), bi.size()));
    molrs::frame_set_column_i32(
        *fr, "bonds", "atomj",
        ::rust::Slice<const std::int32_t>(bj.data(), bj.size()));
    molrs::frame_set_column_i32(
        *fr, "bonds", "type",
        ::rust::Slice<const std::int32_t>(bt.data(), bt.size()));

    const std::vector<std::int32_t> ai = {10}, aj = {11}, ak = {12}, at = {1};
    molrs::frame_set_column_i32(
        *fr, "angles", "atomi",
        ::rust::Slice<const std::int32_t>(ai.data(), ai.size()));
    molrs::frame_set_column_i32(
        *fr, "angles", "atomj",
        ::rust::Slice<const std::int32_t>(aj.data(), aj.size()));
    molrs::frame_set_column_i32(
        *fr, "angles", "atomk",
        ::rust::Slice<const std::int32_t>(ak.data(), ak.size()));
    molrs::frame_set_column_i32(
        *fr, "angles", "type",
        ::rust::Slice<const std::int32_t>(at.data(), at.size()));

    const std::vector<std::int32_t> di = {10}, dj = {11}, dk = {12}, dl = {12},
                                    dt = {1};
    molrs::frame_set_column_i32(
        *fr, "dihedrals", "atomi",
        ::rust::Slice<const std::int32_t>(di.data(), di.size()));
    molrs::frame_set_column_i32(
        *fr, "dihedrals", "atomj",
        ::rust::Slice<const std::int32_t>(dj.data(), dj.size()));
    molrs::frame_set_column_i32(
        *fr, "dihedrals", "atomk",
        ::rust::Slice<const std::int32_t>(dk.data(), dk.size()));
    molrs::frame_set_column_i32(
        *fr, "dihedrals", "atoml",
        ::rust::Slice<const std::int32_t>(dl.data(), dl.size()));
    molrs::frame_set_column_i32(
        *fr, "dihedrals", "type",
        ::rust::Slice<const std::int32_t>(dt.data(), dt.size()));

    // Cubic box, row-major 3x3 H.
    const std::vector<double> h = {20.0, 0.0, 0.0, 0.0, 20.0,
                                   0.0,  0.0, 0.0, 20.0};
    molrs::frame_set_simbox(*fr,
                            ::rust::Slice<const double>(h.data(), h.size()));

    return fr;
}

} // namespace

// ── ac-002: round-trip materialize -> write_back -> read back ───────────

TEST(Materialize, RoundTrip) {
    auto fr = make_full_frame();

    Frame cf = materialize(*fr);

    ASSERT_EQ(cf.atoms.size(), 3);
    int n_bond = 0, n_angle = 0, n_dihedral = 0;
    for (int i = 0; i < cf.atoms.nlocal; ++i) {
        n_bond += cf.atoms.num_bond[i];
        n_angle += cf.atoms.num_angle[i];
        n_dihedral += cf.atoms.num_dihedral[i];
    }
    ASSERT_EQ(n_bond, 1);
    ASSERT_EQ(n_angle, 1);
    ASSERT_EQ(n_dihedral, 1);

    // Mutate the cpu::Frame so write-back has something distinct to carry.
    cf.atoms.x[1] = 7.25;
    cf.atoms.vy[2] = -9.5;
    cf.atoms.fz[0] = 42.0;
    cf.atoms.charge[0] = -0.875;
    cf.scalars.energy = -123.5;
    cf.scalars.kinetic = 64.25;
    cf.scalars.conserved = -59.25;
    for (int k = 0; k < 6; ++k) {
        cf.scalars.stress[k] = 1.0 + k;
        cf.scalars.stress_kinetic[k] = 10.0 + k;
    }

    write_back(cf, *fr);

    // Read every per-atom column back out of the molrs.Frame.
    auto rx = molrs::frame_column_f64(*fr, "atoms", "x");
    auto ry = molrs::frame_column_f64(*fr, "atoms", "y");
    auto rz = molrs::frame_column_f64(*fr, "atoms", "z");
    auto rvx = molrs::frame_column_f64(*fr, "atoms", "vx");
    auto rvy = molrs::frame_column_f64(*fr, "atoms", "vy");
    auto rvz = molrs::frame_column_f64(*fr, "atoms", "vz");
    auto rfx = molrs::frame_column_f64(*fr, "atoms", "fx");
    auto rfy = molrs::frame_column_f64(*fr, "atoms", "fy");
    auto rfz = molrs::frame_column_f64(*fr, "atoms", "fz");
    auto rq = molrs::frame_column_f64(*fr, "atoms", "charge");

    ASSERT_EQ(rx.size(), 3u);
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(rx[i], cf.atoms.x[i]);
        EXPECT_DOUBLE_EQ(ry[i], cf.atoms.y[i]);
        EXPECT_DOUBLE_EQ(rz[i], cf.atoms.z[i]);
        EXPECT_DOUBLE_EQ(rvx[i], cf.atoms.vx[i]);
        EXPECT_DOUBLE_EQ(rvy[i], cf.atoms.vy[i]);
        EXPECT_DOUBLE_EQ(rvz[i], cf.atoms.vz[i]);
        EXPECT_DOUBLE_EQ(rfx[i], cf.atoms.fx[i]);
        EXPECT_DOUBLE_EQ(rfy[i], cf.atoms.fy[i]);
        EXPECT_DOUBLE_EQ(rfz[i], cf.atoms.fz[i]);
        EXPECT_DOUBLE_EQ(rq[i], cf.atoms.charge[i]);
    }

    // Topology survives the materialize (write_back leaves it untouched).
    auto bi = molrs::frame_column_i32(*fr, "bonds", "atomi");
    auto bj = molrs::frame_column_i32(*fr, "bonds", "atomj");
    ASSERT_EQ(bi.size(), 1u);
    EXPECT_EQ(bi[0], 10);
    EXPECT_EQ(bj[0], 11);

    // cpu side: topology scattered into the per-atom layout, keyed by tag.
    // bond (tags 10-11): smaller tag 10 -> row 0 owns; partner tag 11.
    ASSERT_EQ(cf.atoms.num_bond[0], 1);
    EXPECT_EQ(cf.atoms.bond_atom(0, 0), 11);
    EXPECT_EQ(cf.atoms.bond_type(0, 0), 1);
    // angle (10-11-12): center tag 11 -> row 1 owns; k-partner tag is 12.
    ASSERT_EQ(cf.atoms.num_angle[1], 1);
    EXPECT_EQ(cf.atoms.angle_atom(1, 1), 12);
    // dihedral (10-11-12-12): atom-j tag 11 -> row 1 owns; l-partner tag 12.
    ASSERT_EQ(cf.atoms.num_dihedral[1], 1);
    EXPECT_EQ(cf.atoms.dihedral_atom(1, 2), 12);

    // Box round-trips.
    auto h = molrs::frame_simbox(*fr);
    ASSERT_EQ(h.size(), 9u);
    EXPECT_DOUBLE_EQ(h[0], 20.0);
    EXPECT_DOUBLE_EQ(h[4], 20.0);
    EXPECT_DOUBLE_EQ(h[8], 20.0);
    EXPECT_DOUBLE_EQ(cf.box.length(0), 20.0);
    EXPECT_DOUBLE_EQ(cf.box.length(1), 20.0);
    EXPECT_DOUBLE_EQ(cf.box.length(2), 20.0);

    // System scalars round-trip through the dedicated "scalars" block.
    auto se = molrs::frame_column_f64(*fr, "scalars", "energy");
    auto sk = molrs::frame_column_f64(*fr, "scalars", "kinetic");
    auto sc = molrs::frame_column_f64(*fr, "scalars", "conserved");
    ASSERT_EQ(se.size(), 1u);
    EXPECT_DOUBLE_EQ(se[0], -123.5);
    EXPECT_DOUBLE_EQ(sk[0], 64.25);
    EXPECT_DOUBLE_EQ(sc[0], -59.25);

    auto ss = molrs::frame_column_f64(*fr, "stress", "value");
    auto ssk = molrs::frame_column_f64(*fr, "stress_kinetic", "value");
    ASSERT_EQ(ss.size(), 6u);
    ASSERT_EQ(ssk.size(), 6u);
    for (int k = 0; k < 6; ++k) {
        EXPECT_DOUBLE_EQ(ss[k], 1.0 + k);
        EXPECT_DOUBLE_EQ(ssk[k], 10.0 + k);
    }
}

// ── ac-003: required fields fail loud; optional fields default to zero ──

TEST(Materialize, MissingPositionsThrows) {
    auto fr = molrs::frame_new();
    // type present, x/y/z absent.
    const std::vector<std::int32_t> type = {1, 2};
    molrs::frame_set_column_i32(
        *fr, "atoms", "type",
        ::rust::Slice<const std::int32_t>(type.data(), type.size()));

    try {
        Frame cf = materialize(*fr);
        FAIL() << "materialize must throw when atoms.x is absent";
    } catch (const std::exception& e) {
        EXPECT_NE(std::string(e.what()).find("x"), std::string::npos)
            << "exception message must name the missing field: " << e.what();
    }
}

TEST(Materialize, MissingTypeThrows) {
    auto fr = molrs::frame_new();
    // positions present, type absent.
    const std::vector<double> x = {0.0, 1.0};
    molrs::frame_set_column_f64(
        *fr, "atoms", "x", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "y", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "z", ::rust::Slice<const double>(x.data(), x.size()));

    try {
        Frame cf = materialize(*fr);
        FAIL() << "materialize must throw when atoms.type is absent";
    } catch (const std::exception& e) {
        EXPECT_NE(std::string(e.what()).find("type"), std::string::npos)
            << "exception message must name the missing field: " << e.what();
    }
}

TEST(Materialize, MissingOptionalFieldsDefaultToZero) {
    auto fr = molrs::frame_new();
    // Only required fields: x, y, z, type. No velocities, no charge.
    const std::vector<double> x = {0.0, 1.0, 2.0};
    const std::vector<std::int32_t> type = {1, 1, 1};
    molrs::frame_set_column_f64(
        *fr, "atoms", "x", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "y", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "z", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_i32(
        *fr, "atoms", "type",
        ::rust::Slice<const std::int32_t>(type.data(), type.size()));

    Frame cf = materialize(*fr); // must NOT throw

    ASSERT_EQ(cf.atoms.size(), 3);
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(cf.atoms.vx[i], 0.0);
        EXPECT_DOUBLE_EQ(cf.atoms.vy[i], 0.0);
        EXPECT_DOUBLE_EQ(cf.atoms.vz[i], 0.0);
        EXPECT_DOUBLE_EQ(cf.atoms.fx[i], 0.0);
        EXPECT_DOUBLE_EQ(cf.atoms.charge[i], 0.0);
    }
    // No topology -> per-atom counts all zero.
    for (int i = 0; i < cf.atoms.nlocal; ++i) {
        EXPECT_EQ(cf.atoms.num_bond[i], 0);
        EXPECT_EQ(cf.atoms.num_angle[i], 0);
        EXPECT_EQ(cf.atoms.num_dihedral[i], 0);
    }
}

// ── ac-004: a materialized cpu::Frame is runnable under cpu::md::MD ─────

TEST(Materialize, MaterializedFrameIsRunnable) {
    auto fr = molrs::frame_new();

    // Two argon-like atoms separated by ~sigma in a roomy cubic box.
    const std::vector<double> x = {0.0, 3.4};
    const std::vector<double> y = {0.0, 0.0};
    const std::vector<double> z = {0.0, 0.0};
    const std::vector<double> mass = {39.948, 39.948};
    const std::vector<std::int32_t> type = {1, 1};
    molrs::frame_set_column_f64(
        *fr, "atoms", "x", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "y", ::rust::Slice<const double>(y.data(), y.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "z", ::rust::Slice<const double>(z.data(), z.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "mass",
        ::rust::Slice<const double>(mass.data(), mass.size()));
    molrs::frame_set_column_i32(
        *fr, "atoms", "type",
        ::rust::Slice<const std::int32_t>(type.data(), type.size()));
    const std::vector<double> h = {30.0, 0.0, 0.0, 0.0, 30.0,
                                   0.0,  0.0, 0.0, 30.0};
    molrs::frame_set_simbox(*fr,
                            ::rust::Slice<const double>(h.data(), h.size()));

    Frame cf = materialize(*fr);
    ASSERT_EQ(cf.atoms.size(), 2);

    ATV::cpu::md::MD sim;
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{6.8, /*skin=*/0.5});
    sim.set_pair(ATV::cpu::PairLJCut{0.0104, 3.4, 6.8});
    sim.set_integrator(ATV::cpu::md::VelocityVerlet{/*dt=*/0.5});
    sim.build();
    sim.run(cf, 10);

    EXPECT_TRUE(std::isfinite(cf.scalars.energy))
        << "MD run on a materialized frame must leave a finite energy";
}
