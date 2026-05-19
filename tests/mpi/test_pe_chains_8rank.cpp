/// @file tests/mpi/test_pe_chains_8rank.cpp
/// @brief Phase 0.3 / Phase B.2 (msu-pe-validation, internal-consistency
///        variant) -- 8-rank PE bonded NVE conservation.
///
/// Builds 8 PE chains, one centred inside each rank's slab in a 2×2×2
/// Cartesian partition, each chain laid along y so it fits inside a single
/// `slab_y` extent without straddling a rank boundary. Drives 1000
/// VelocityVerlet steps with HarmonicBond + HarmonicAngle +
/// MultiHarmonicDihedral and asserts the **globally** allreduced total
/// energy drifts by less than 1e-3 over the run.
///
/// Per-chain dynamics are decoupled (no LJ, no inter-chain bonded), so
/// global conservation is a sufficient check that:
///   1. the topology filter places every chain's bonded set on its rank,
///   2. multi-rank integration + bonded force assembly does not introduce
///      drift beyond a single-rank baseline,
///   3. all 8 ranks compute consistently (otherwise rank-local cancellation
///      would mask per-rank drift in the allreduce sum).

#include <gtest/gtest.h>

#include <mpi.h>

#include <cmath>
#include <vector>

#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
#include "atomiverse/cpu/potentials/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"

using namespace ATV;

namespace {

inline constexpr double kKcalPerMolToEv = 1.0 / 23.060547830619026;
inline constexpr double kBondK_kcal  = 268.0;
inline constexpr double kBondR0      = 1.53;
inline constexpr double kAngleK_kcal = 60.0;
inline constexpr double kAngleTheta0 = 109.5 * kPI / 180.0;

inline double bond_K_eV()  { return kBondK_kcal  * kKcalPerMolToEv; }
inline double angle_K_eV() { return kAngleK_kcal * kKcalPerMolToEv; }
inline std::vector<double> dihedral_Cn_eV() {
    return {
        1.73    * kKcalPerMolToEv,   // C0  (MSU paper, Tschopp 2010)
        -4.49   * kKcalPerMolToEv,   // C1
        0.776   * kKcalPerMolToEv,   // C2
        6.99    * kKcalPerMolToEv,   // C3
        0.0,                          // C4
    };
}

// Build a 2×2×2 box with one C18H38 chain per rank slab.
//
// Box edges all equal `box_len`, slab edges `box_len/2`. Chain laid along
// y, length (n_carbons - 1) * dy_step. The chain is centred in (x, z) at
// the slab centre, and its y-span is positioned with `y_margin` slack on
// both sides so no atom sits within `y_margin` of the slab boundary.
// Atom and topology indexing is contiguous chain-by-chain in (cx, cy, cz)
// row-major order.
cpu::Frame make_8chain_box(double box_len, int n_carbons, double y_margin) {
    const int n_chains = 8;
    const int n_atoms  = n_chains * n_carbons;
    cpu::Frame frame(n_atoms);

    const double half_apex = (kPI - kAngleTheta0) / 2.0;
    const double dy_step   = kBondR0 * std::cos(half_apex);
    const double dz_step   = kBondR0 * std::sin(half_apex);

    const double slab = box_len / 2.0;
    const double chain_len = (n_carbons - 1) * dy_step;
    if (chain_len + 2 * y_margin > slab) {
        // Will trip the partition check; surface a clear error early.
        throw std::runtime_error("chain length > slab - 2 margin -- shrink chain or grow box");
    }

    int idx = 0;
    for (int cx = 0; cx < 2; ++cx) {
        for (int cy = 0; cy < 2; ++cy) {
            for (int cz = 0; cz < 2; ++cz) {
                const double sub_x_lo = cx * slab;
                const double sub_y_lo = cy * slab;
                const double sub_z_lo = cz * slab;
                const double cx0 = sub_x_lo + slab / 2.0;
                const double cz0 = sub_z_lo + slab / 2.0;
                const double y0  = sub_y_lo + y_margin;
                for (int i = 0; i < n_carbons; ++i) {
                    frame.atoms.x[idx] = cx0;
                    frame.atoms.y[idx] = y0 + i * dy_step;
                    frame.atoms.z[idx] = cz0 + (i % 2 == 0 ? 0.0 : dz_step);
                    frame.atoms.mass[idx] = 14.027;
                    frame.atoms.type[idx] = 0;
                    // Tiny deterministic transverse kick along x; sums to
                    // zero on each chain so per-chain COM does not drift.
                    const double sign = (i % 2 == 0) ? +1.0 : -1.0;
                    frame.atoms.vx[idx] = sign * 1e-3;
                    ++idx;
                }
            }
        }
    }

    // Pre-allocate to TOTAL topology size and write every entry exactly
    // once. (xt::xtensor::resize() does not preserve previous data on grow,
    // so an append-style construction would silently zero earlier chains.)
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

    frame.box = Mat3::cubic(box_len);
    return frame;
}

double total_energy_global(const cpu::Frame& frame) {
    // ReduceEnergy already wrote the global PE/KE back into frame, so
    // every rank has the same value here.
    return frame.energy + frame.kinetic;
}

} // namespace

TEST(PeChains8Rank, BondedNveConservation) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 8) {
        GTEST_SKIP() << "8-rank specific test; got world_size=" << world_size;
    }

    // box 48 → slab 24 → chain length ~21.6 Å (18 atoms × 1.272 dy_step) +
    // 1 Å margin both sides ≤ 24. Comfortable.
    auto frame = make_8chain_box(48.0, 18, 1.0);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(ATV::md::Phase::force)
        .add(HarmonicBond{{bond_K_eV()}, {kBondR0}})
        .add(HarmonicAngle{{angle_K_eV()}, {kAngleTheta0}})
        .add(MultiHarmonicDihedral{{dihedral_Cn_eV()}});
    sim.build();

    sim.run(frame, 1);
    const double e0 = total_energy_global(frame);
    ASSERT_GT(std::abs(e0), 0.0);
    ASSERT_FALSE(std::isnan(e0));

    sim.run(frame, 999);
    const double e1 = total_energy_global(frame);
    const double drift = std::abs(e1 - e0) / std::abs(e0);
    EXPECT_LT(drift, 1e-3)
        << "8-rank bonded NVE drift " << drift
        << " exceeds 1e-3 over 1000 steps (E0=" << e0 << ", E1=" << e1 << ")";

    // No NaN explosion on any rank.
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
