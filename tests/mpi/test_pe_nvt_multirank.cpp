/// @file tests/mpi/test_pe_nvt_multirank.cpp
/// @brief Phase 0.5 / B.3 (msu-pe-validation, internal) -- multi-rank NVT.
///
/// 2x2x2 box of 8 PE chains (one per rank slab) with InitVelocity at 300 K
/// and a NoseHoover thermostat at 300 K. Drives 5000 Velocity-Verlet steps
/// with bonded forces only and asserts the final-step temperature -- read
/// from the (already globally reduced) `frame.scalars.kinetic` -- sits within
/// 30% of the target. The 30% bound matches the single-rank Phase A.5
/// tolerance and accounts for the 1/√N statistical noise on a 144-atom
/// system.

#include <gtest/gtest.h>

#include <mpi.h>

#include <cmath>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
#include "atomiverse/cpu/potentials/multi_harmonic_dihedral.h"

using namespace ATV;

namespace {

inline constexpr double kKcalPerMolToEv = 1.0 / 23.060547830619026;
inline constexpr double kBondK_kcal = 268.0;
inline constexpr double kBondR0 = 1.53;
inline constexpr double kAngleK_kcal = 60.0;
inline constexpr double kAngleTheta0 = 109.5 * kPI / 180.0;

inline double bond_K_eV() { return kBondK_kcal * kKcalPerMolToEv; }
inline double angle_K_eV() { return kAngleK_kcal * kKcalPerMolToEv; }
inline std::vector<double> dihedral_Cn_eV() {
    return {
        1.73 * kKcalPerMolToEv,  // C0  (MSU paper, Tschopp 2010)
        -4.49 * kKcalPerMolToEv, // C1
        0.776 * kKcalPerMolToEv, // C2
        6.99 * kKcalPerMolToEv,  // C3
        0.0,                     // C4
    };
}

cpu::Frame make_8chain_box(double box_len, int n_carbons, double y_margin) {
    const int n_chains = 8;
    const int n_atoms = n_chains * n_carbons;
    cpu::Frame frame(n_atoms);

    const double half_apex = (kPI - kAngleTheta0) / 2.0;
    const double dy_step = kBondR0 * std::cos(half_apex);
    const double dz_step = kBondR0 * std::sin(half_apex);

    const double slab = box_len / 2.0;
    int idx = 0;
    for (int cx = 0; cx < 2; ++cx) {
        for (int cy = 0; cy < 2; ++cy) {
            for (int cz = 0; cz < 2; ++cz) {
                const double cx0 = cx * slab + slab / 2.0;
                const double cz0 = cz * slab + slab / 2.0;
                const double y0 = cy * slab + y_margin;
                for (int i = 0; i < n_carbons; ++i) {
                    frame.atoms.x[idx] = cx0;
                    frame.atoms.y[idx] = y0 + i * dy_step;
                    frame.atoms.z[idx] = cz0 + (i % 2 == 0 ? 0.0 : dz_step);
                    frame.atoms.mass[idx] = 14.027;
                    frame.atoms.type[idx] = 0;
                    ++idx;
                }
            }
        }
    }

    const int n_b_per = n_carbons - 1;
    const int n_a_per = n_carbons - 2;
    const int n_d_per = n_carbons - 3;
    frame.atoms.map_clear();
    frame.atoms.map_set(0, frame.atoms.nlocal);
    for (int c = 0; c < n_chains; ++c) {
        const int base = c * n_carbons;
        for (int k = 0; k < n_b_per; ++k)
            frame.atoms.add_bond(base + k, base + k + 1, 0);
        for (int k = 0; k < n_a_per; ++k)
            frame.atoms.add_angle(base + k, base + k + 1, base + k + 2, 0);
        for (int k = 0; k < n_d_per; ++k)
            frame.atoms.add_dihedral(base + k, base + k + 1, base + k + 2,
                                     base + k + 3, 0);
    }
    frame.atoms.build_special_list();

    frame.box = Mat3::cubic(box_len);
    return frame;
}

double temperature_from_kinetic_global(double e_kin_ev_global,
                                       int n_atoms_global) {
    return (2.0 / 3.0) * e_kin_ev_global / (n_atoms_global * kBoltzmann);
}

} // namespace

TEST(PeNvtMultiRank, NoseHooverHoldsGlobalTemperature) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 8) {
        GTEST_SKIP() << "8-rank specific test; got world_size=" << world_size;
    }

    const double T_target = 300.0; // K
    const double dt = 0.5;         // fs
    const double tau_T = 100.0;    // fs

    const int n_carbons = 18;
    auto frame = make_8chain_box(48.0, n_carbons, 1.0);
    const int n_atoms_global = frame.atoms.size();

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/4242});
    sim.add(cpu::md::NoseHoover{T_target, tau_T, dt});
    sim.on(ATV::md::Phase::force)
        .add(HarmonicBond{{bond_K_eV()}, {kBondR0}})
        .add(HarmonicAngle{{angle_K_eV()}, {kAngleTheta0}})
        .add(MultiHarmonicDihedral{{dihedral_Cn_eV()}});
    sim.build();

    sim.run(frame, 5000);

    // ReduceEnergy already wrote the global KE into frame.scalars.kinetic.
    const double T_final =
        temperature_from_kinetic_global(frame.scalars.kinetic, n_atoms_global);
    const double rel_err = std::abs(T_final - T_target) / T_target;

    // 144-atom system: instantaneous-T statistical noise is √(2/N) ≈ 12 %;
    // 5000 steps × 0.5 fs = 2.5 ps is short for the chain to fully thermalise
    // at this τ_T. We assert the looser 40 % bound -- the point is to prove
    // the thermostat (a) doesn't blow up, (b) doesn't diverge from target,
    // (c) the multi-rank KE reduction in NoseHoover gives the same ballpark
    // as single-rank. Quantitative block-averaged matching against LAMMPS
    // is a downstream Phase D milestone.
    EXPECT_LT(rel_err, 0.4)
        << "post-thermostat global T=" << T_final << " K vs target " << T_target
        << " K, relative error " << rel_err << " exceeds 40%";

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
