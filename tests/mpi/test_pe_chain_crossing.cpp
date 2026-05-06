/// @file tests/mpi/test_pe_chain_crossing.cpp
/// @brief Phase 0.3 final -- a single PE chain that genuinely crosses a
///        rank boundary, with bonded forces on both sides.
///
/// Prior multi-rank PE tests rigged the geometry so each chain sat fully
/// inside one rank's slab. That config never exercised the cross-boundary
/// part of the bonded force pipeline. This test fixes that:
///
///   * 1 chain of 20 united-atom CH2 sites laid along x.
///   * Box 30 Å × 30 Å × 30 Å, dims = (2, 1, 1) → slab_x = 15 Å.
///   * Chain atoms span x ≈ 3 → 27.2 Å; rank 0 owns atoms whose x < 15,
///     rank 1 owns the rest. Roughly half-and-half.
///   * HaloExchangeForward cutoff = 5 Å so a dihedral that spans 3 bond
///     lengths can resolve every atom on at least one rank as either
///     owner or ghost.
///
/// Asserts the global energy stays finite and bounded over a 200-step NVE
/// run -- the chain's own bonded forces should hold it together across
/// the rank boundary.

#include <gtest/gtest.h>

#include <mpi.h>

#include <cmath>
#include <vector>

#include "atomiverse/cpu/bonded/harmonic_angle.h"
#include "atomiverse/cpu/bonded/harmonic_bond.h"
#include "atomiverse/cpu/bonded/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"

using namespace ATV;

namespace {

inline constexpr double kKcalPerMolToEv = 1.0 / 23.060547830619026;
inline constexpr double kBondK_kcal  = 268.0;
inline constexpr double kBondR0      = 1.529;
inline constexpr double kAngleK_kcal = 58.35;
inline constexpr double kAngleTheta0 = 112.7 * kPI / 180.0;
inline constexpr double kV1_kcal = 1.411;
inline constexpr double kV2_kcal = -0.271;
inline constexpr double kV3_kcal = 3.145;

inline double bond_K_eV()  { return kBondK_kcal  * kKcalPerMolToEv; }
inline double angle_K_eV() { return kAngleK_kcal * kKcalPerMolToEv; }
inline std::vector<double> dihedral_Cn_eV() {
    const double v1 = kV1_kcal * kKcalPerMolToEv;
    const double v2 = kV2_kcal * kKcalPerMolToEv;
    const double v3 = kV3_kcal * kKcalPerMolToEv;
    return {
        v1 / 2.0 + v2 + v3 / 2.0,
        v1 / 2.0 - 1.5 * v3,
        -v2,
         2.0 * v3,
         0.0,
    };
}

// Chain laid along x: atom positions form an all-trans planar zigzag,
// with the chain axis pointing in +x. The y/z zigzag gives the angle/
// dihedral non-degenerate geometry; chain centre is placed so it spans
// the x = box_x/2 partition boundary.
cpu::Frame make_chain_along_x(int n_carbons, double box_x, double box_yz) {
    cpu::Frame frame(n_carbons);

    const double half_apex = (kPI - kAngleTheta0) / 2.0;
    const double dx_step = kBondR0 * std::cos(half_apex);
    const double dy_step = kBondR0 * std::sin(half_apex);

    const double chain_len = (n_carbons - 1) * dx_step;
    const double x_start = (box_x - chain_len) / 2.0;
    const double y0 = box_yz / 2.0;
    const double z0 = box_yz / 2.0;

    for (int i = 0; i < n_carbons; ++i) {
        frame.atoms.x[i] = x_start + i * dx_step;
        frame.atoms.y[i] = y0 + (i % 2 == 0 ? 0.0 : dy_step);
        frame.atoms.z[i] = z0;
        frame.atoms.mass[i] = 14.027;
        frame.atoms.type[i] = 0;
    }

    const int n_b = n_carbons - 1;
    frame.bonds.atomi = xt::zeros<int>({n_b});
    frame.bonds.atomj = xt::zeros<int>({n_b});
    frame.bonds.type  = xt::zeros<int>({n_b});
    for (int k = 0; k < n_b; ++k) {
        frame.bonds.atomi[k] = k;
        frame.bonds.atomj[k] = k + 1;
    }
    const int n_a = n_carbons - 2;
    frame.angles.atomi = xt::zeros<int>({n_a});
    frame.angles.atomj = xt::zeros<int>({n_a});
    frame.angles.atomk = xt::zeros<int>({n_a});
    frame.angles.type  = xt::zeros<int>({n_a});
    for (int k = 0; k < n_a; ++k) {
        frame.angles.atomi[k] = k;
        frame.angles.atomj[k] = k + 1;
        frame.angles.atomk[k] = k + 2;
    }
    const int n_d = n_carbons - 3;
    frame.dihedrals.atomi = xt::zeros<int>({n_d});
    frame.dihedrals.atomj = xt::zeros<int>({n_d});
    frame.dihedrals.atomk = xt::zeros<int>({n_d});
    frame.dihedrals.atoml = xt::zeros<int>({n_d});
    frame.dihedrals.type  = xt::zeros<int>({n_d});
    for (int k = 0; k < n_d; ++k) {
        frame.dihedrals.atomi[k] = k;
        frame.dihedrals.atomj[k] = k + 1;
        frame.dihedrals.atomk[k] = k + 2;
        frame.dihedrals.atoml[k] = k + 3;
    }

    frame.box = Mat3::diagonal(box_x, box_yz, box_yz);
    return frame;
}

} // namespace

TEST(PeChainCrossing, BondedNveHandlesChainAcrossBoundary) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 2) {
        GTEST_SKIP() << "2-rank specific test; got world_size=" << world_size;
    }

    auto frame = make_chain_along_x(20, 30.0, 30.0);

    // Pair cutoff = 5 Å on its own covers the dihedral span (3 × 1.529 ≈
    // 4.6 Å). The spec-A.5 auto-extension added in `build()` makes this
    // redundant here; the explicit 5.0 preserves the historical baseline
    // so a future regression in the auto-extension cannot silently mask a
    // cross-rank bonded miss through this test.
    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(ATV::md::Phase::pre_force)
        .add(cpu::md::mpi::HaloExchangeForward{5.0}, -20);
    sim.on(ATV::md::Phase::post_force)
        .add(cpu::md::mpi::HaloExchangeReverse{}, -1000);
    sim.on(ATV::md::Phase::force)
        .add(cpu::bonded::HarmonicBond{{bond_K_eV()}, {kBondR0}})
        .add(cpu::bonded::HarmonicAngle{{angle_K_eV()}, {kAngleTheta0}})
        .add(cpu::bonded::MultiHarmonicDihedral{{dihedral_Cn_eV()}});
    sim.build();

    sim.run(frame, 1);
    const double e0 = frame.energy + frame.kinetic;   // already global via ReduceEnergy
    ASSERT_FALSE(std::isnan(e0)) << "step-1 global energy is NaN";

    sim.run(frame, 199);
    const double e1 = frame.energy + frame.kinetic;
    EXPECT_FALSE(std::isnan(e1)) << "global energy went NaN";
    EXPECT_LT(std::abs(e1), 1e3)
        << "global energy " << e1 << " is wildly large -- cross-boundary "
        << "bonded resolve / halo skin probably wrong (E0=" << e0 << ")";

    // No NaN in any owner position.
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    for (int i = 0; i < frame.atoms.size(); ++i) {
        ASSERT_FALSE(std::isnan(frame.atoms.x[i]))
            << "rank " << my_rank << " atom " << i << " x went NaN";
    }
}

// Same chain + same bonded FF, but the user passes a PAIR cutoff well
// below the dihedral span (0.1 Å). The auto-extension in build() scans
// bonds / angles / dihedrals and must add enough ghost skin on top so
// every cross-rank tuple still resolves. Without the refactor this run
// would drop bonded tuples at the boundary and blow up within a handful
// of steps.
TEST(PeChainCrossing, BondedSkinAutoExtendsWithSmallPairCutoff) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 2) {
        GTEST_SKIP() << "2-rank specific test; got world_size=" << world_size;
    }

    auto frame = make_chain_along_x(20, 30.0, 30.0);

    const double tiny_pair_cutoff = 0.1;  // << bond length 1.529 Å
    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(ATV::md::Phase::pre_force)
        .add(cpu::md::mpi::HaloExchangeForward{tiny_pair_cutoff}, -20);
    sim.on(ATV::md::Phase::post_force)
        .add(cpu::md::mpi::HaloExchangeReverse{}, -1000);
    sim.on(ATV::md::Phase::force)
        .add(cpu::bonded::HarmonicBond{{bond_K_eV()}, {kBondR0}})
        .add(cpu::bonded::HarmonicAngle{{angle_K_eV()}, {kAngleTheta0}})
        .add(cpu::bonded::MultiHarmonicDihedral{{dihedral_Cn_eV()}});
    sim.build();

    sim.run(frame, 50);
    const double e = frame.energy + frame.kinetic;
    EXPECT_FALSE(std::isnan(e)) << "global energy went NaN";
    EXPECT_LT(std::abs(e), 1e3)
        << "global energy " << e << " blew up -- auto-extension failed "
        << "to cover the cross-rank dihedral span";

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    for (int i = 0; i < frame.atoms.size(); ++i) {
        ASSERT_FALSE(std::isnan(frame.atoms.x[i]))
            << "rank " << my_rank << " atom " << i << " x went NaN";
    }
}

int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    ::testing::InitGoogleTest(&argc, argv);
    const int rc = RUN_ALL_TESTS();
    MPI_Finalize();
    return rc;
}
