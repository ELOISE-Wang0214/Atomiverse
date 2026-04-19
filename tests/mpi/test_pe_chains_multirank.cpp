/// @file tests/mpi/test_pe_chains_multirank.cpp
/// @brief Phase 0.3 MVP -- multi-rank PE bonded forces.
///
/// Two PE chains (C20H42 each) placed entirely inside two adjacent ranks'
/// slabs along x, so no bond/angle/dihedral straddles the slab boundary
/// and the partition-time topology filter (msu-pe-validation Phase 0.3
/// MVP) can keep every owned tuple in its rank's local subset. Bonds
/// using ghost atoms are deferred to the Phase 0.3-extended milestone.
///
/// Drives a short bonded-only NVE run and asserts the global energy stays
/// finite (no NaN, no runaway). Cross-validation against the single-rank
/// run on the same chains is left to a downstream commit because Phase 2
/// has no momentum-zero initialisation across ranks.

#include <gtest/gtest.h>

#include <mpi.h>

#include <cmath>
#include <vector>

#include "atomiverse/cpu/bonded/harmonic_angle.h"
#include "atomiverse/cpu/bonded/harmonic_bond.h"
#include "atomiverse/cpu/bonded/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
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

// Build a frame with two C20H42 chains, each sitting fully inside one
// half of the box along x. With 2 MPI ranks and dims = (2,1,1) the
// partition assigns chain 0 (centred at x = box_x/4) to rank 0 and chain 1
// (centred at x = 3·box_x/4) to rank 1.
cpu::Frame make_two_chain_box(double box_x, double box_yz) {
    const int n_carbons = 20;
    const int n_atoms   = 2 * n_carbons;
    cpu::Frame frame(n_atoms);

    const double half_apex = (kPI - kAngleTheta0) / 2.0;
    const double dy_step = kBondR0 * std::cos(half_apex);   // along y now
    const double dz_step = kBondR0 * std::sin(half_apex);   // small z zigzag

    auto place_chain = [&](int base, double cx) {
        // Chain laid along y, centred at (cx, box_yz/2, box_yz/2).
        const double y0 = (box_yz - (n_carbons - 1) * dy_step) / 2.0;
        for (int i = 0; i < n_carbons; ++i) {
            frame.atoms.x[base + i] = cx;
            frame.atoms.y[base + i] = y0 + i * dy_step;
            frame.atoms.z[base + i] = box_yz / 2.0 + (i % 2 == 0 ? 0.0 : dz_step);
            frame.atoms.mass[base + i] = 14.027;
            frame.atoms.type[base + i] = 0;
        }
    };
    place_chain(0,           box_x * 0.25);   // owned by rank 0 (x in [0, box_x/2))
    place_chain(n_carbons,   box_x * 0.75);   // owned by rank 1 (x in [box_x/2, box_x))

    // Tiny opposite-sign kick per chain so we exercise the integrator
    // without driving atoms across the slab boundary in 200 steps.
    for (int c = 0; c < 2; ++c) {
        for (int i = 0; i < n_carbons; ++i) {
            const double sign = (i % 2 == 0) ? +1.0 : -1.0;
            frame.atoms.vz[c * n_carbons + i] = sign * 1e-3;
        }
    }

    // Pre-allocate to TOTAL size for both chains, then fill all entries.
    // (xt::xtensor::resize wipes the existing data on grow, so per-chain
    // append doesn't work without manually preserving previous writes.)
    const int n_chains = 2;
    const int n_b_per = n_carbons - 1;
    const int n_a_per = n_carbons - 2;
    const int n_d_per = n_carbons - 3;
    const int n_b_tot = n_chains * n_b_per;
    const int n_a_tot = n_chains * n_a_per;
    const int n_d_tot = n_chains * n_d_per;

    frame.bonds.atomi = xt::zeros<int>({n_b_tot});
    frame.bonds.atomj = xt::zeros<int>({n_b_tot});
    frame.bonds.type  = xt::zeros<int>({n_b_tot});
    frame.angles.atomi = xt::zeros<int>({n_a_tot});
    frame.angles.atomj = xt::zeros<int>({n_a_tot});
    frame.angles.atomk = xt::zeros<int>({n_a_tot});
    frame.angles.type  = xt::zeros<int>({n_a_tot});
    frame.dihedrals.atomi = xt::zeros<int>({n_d_tot});
    frame.dihedrals.atomj = xt::zeros<int>({n_d_tot});
    frame.dihedrals.atomk = xt::zeros<int>({n_d_tot});
    frame.dihedrals.atoml = xt::zeros<int>({n_d_tot});
    frame.dihedrals.type  = xt::zeros<int>({n_d_tot});

    for (int c = 0; c < n_chains; ++c) {
        const int base = c * n_carbons;
        for (int k = 0; k < n_b_per; ++k) {
            const int b = c * n_b_per + k;
            frame.bonds.atomi[b] = base + k;
            frame.bonds.atomj[b] = base + k + 1;
        }
        for (int k = 0; k < n_a_per; ++k) {
            const int a = c * n_a_per + k;
            frame.angles.atomi[a] = base + k;
            frame.angles.atomj[a] = base + k + 1;
            frame.angles.atomk[a] = base + k + 2;
        }
        for (int k = 0; k < n_d_per; ++k) {
            const int d = c * n_d_per + k;
            frame.dihedrals.atomi[d] = base + k;
            frame.dihedrals.atomj[d] = base + k + 1;
            frame.dihedrals.atomk[d] = base + k + 2;
            frame.dihedrals.atoml[d] = base + k + 3;
        }
    }

    frame.box = Mat3{};
    frame.box(0, 0) = box_x;
    frame.box(1, 1) = box_yz;
    frame.box(2, 2) = box_yz;
    return frame;
}

double total_energy_allreduced(const cpu::Frame& frame) {
    double local = frame.energy + frame.kinetic;
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

} // namespace

TEST(PeChainsMultiRank, BondedOnlyNveStable) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 2) {
        GTEST_SKIP() << "2-rank specific test; got world_size=" << world_size;
    }

    auto frame = make_two_chain_box(40.0, 30.0);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(ATV::md::Phase::force)
        .add(cpu::bonded::HarmonicBond{{bond_K_eV()}, {kBondR0}})
        .add(cpu::bonded::HarmonicAngle{{angle_K_eV()}, {kAngleTheta0}})
        .add(cpu::bonded::MultiHarmonicDihedral{{dihedral_Cn_eV()}});
    sim.build();

    sim.run(frame, 200);

    const double e_global = total_energy_allreduced(frame);
    EXPECT_FALSE(std::isnan(e_global)) << "global energy went NaN";
    EXPECT_LT(std::abs(e_global), 1e6)
        << "global energy " << e_global
        << " is wildly large -- topology filter or bonded kernel is broken";

    // Per-rank: owners must remain finite.
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
