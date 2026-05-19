/// @file tests/mpi/test_pe_msu_pipeline.cpp
/// @brief End-to-end multi-rank PE simulation pipeline.
///
/// Combines every CPU MD component the MSU amorphous-PE protocol needs:
/// bonded forces, LJ with special-bonds exclusions, halo, atom migration,
/// NoseHoover NVT, Berendsen Barostat, and Deform uniaxial tension. Runs
/// under 8 MPI ranks (2×2×2 partition).
///
/// This is **not** the full MSU acceptance test (no LAMMPS cross-check,
/// no 50 000 atoms, no quantitative stress-strain comparison). It is the
/// first integration test that lets the entire multi-rank PE stack run
/// for thousands of steps with every moving part live, and asserts the
/// stack does not blow up. Quantitative validation (Phase D) requires a
/// LAMMPS reference and is tracked under msu-pe-validation §D.

#include <gtest/gtest.h>

#include <mpi.h>

#include <cmath>
#include <vector>

#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
#include "atomiverse/cpu/potentials/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/barostat.h"
#include "atomiverse/cpu/md/deform.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/virial_stress.h"
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

cpu::Frame make_pe_box(int n_chains_per_axis, int n_carbons_per_chain,
                        double box_len, double margin) {
    const int n_chains = n_chains_per_axis * n_chains_per_axis * n_chains_per_axis;
    const int n_atoms  = n_chains * n_carbons_per_chain;
    cpu::Frame frame(n_atoms);

    const double half_apex = (kPI - kAngleTheta0) / 2.0;
    const double dy_step = kBondR0 * std::cos(half_apex);
    const double dz_step = kBondR0 * std::sin(half_apex);

    const double cell = box_len / n_chains_per_axis;
    int idx = 0;
    for (int cx = 0; cx < n_chains_per_axis; ++cx) {
        for (int cy = 0; cy < n_chains_per_axis; ++cy) {
            for (int cz = 0; cz < n_chains_per_axis; ++cz) {
                const double cx0 = cx * cell + cell / 2.0;
                const double cz0 = cz * cell + cell / 2.0;
                const double y0  = cy * cell + margin;
                for (int i = 0; i < n_carbons_per_chain; ++i) {
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

    const int n_b_per = n_carbons_per_chain - 1;
    const int n_a_per = n_carbons_per_chain - 2;
    const int n_d_per = n_carbons_per_chain - 3;
    frame.bonds.atomi = xt::zeros<int>({n_chains * n_b_per});
    frame.bonds.atomj = xt::zeros<int>({n_chains * n_b_per});
    frame.bonds.type  = xt::zeros<int>({n_chains * n_b_per});
    frame.angles.atomi = xt::zeros<int>({n_chains * n_a_per});
    frame.angles.atomj = xt::zeros<int>({n_chains * n_a_per});
    frame.angles.atomk = xt::zeros<int>({n_chains * n_a_per});
    frame.angles.type  = xt::zeros<int>({n_chains * n_a_per});
    frame.dihedrals.atomi = xt::zeros<int>({n_chains * n_d_per});
    frame.dihedrals.atomj = xt::zeros<int>({n_chains * n_d_per});
    frame.dihedrals.atomk = xt::zeros<int>({n_chains * n_d_per});
    frame.dihedrals.atoml = xt::zeros<int>({n_chains * n_d_per});
    frame.dihedrals.type  = xt::zeros<int>({n_chains * n_d_per});
    for (int c = 0; c < n_chains; ++c) {
        const int base = c * n_carbons_per_chain;
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

} // namespace

// ── Full MSU-class pipeline at 8 ranks ──────────────────────────────

TEST(PeMsuPipeline, NvtPlusNptPlusDeform) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 8) {
        GTEST_SKIP() << "8-rank specific test; got world_size=" << world_size;
    }

    // 8 chains × 18 atoms = 144 atoms in a 48 Å box (slab 24 Å, chain
    // length ~21.6 Å fits inside one slab).
    auto frame = make_pe_box(/*chains_per_axis=*/2, /*n_carbons_per_chain=*/18,
                              /*box_len=*/48.0, /*margin=*/1.0);

    const double T_target = 300.0;     // K
    const double P_target = 0.0;       // eV/Å³
    const double dt       = 0.5;       // fs
    const double tau_T    = 100.0;     // fs
    const double tau_P    = 500.0;     // fs (slower than T_target)
    const double strain_rate = 5e-5;   // 1/fs (~5e10 / s -- tension-like)

    const double lj_eps = 0.00512;
    const double lj_sig = 4.01;
    const double lj_cut = 10.5;

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/12345});
    sim.add(cpu::md::NoseHoover{T_target, tau_T, dt});
    sim.on(ATV::md::Phase::pre_force)
        .add(cpu::md::mpi::HaloExchangeForward{lj_cut}, -20);
    sim.on(ATV::md::Phase::post_force)
        .add(cpu::md::mpi::HaloExchangeReverse{}, -1000);
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::NeighborList<>{lj_cut}, -10)
        .add(cpu::md::LJ{lj_eps, lj_sig, lj_cut}.set_special_bonds(0.0, 0.0, 0.0))
        .add(HarmonicBond{{bond_K_eV()}, {kBondR0}})
        .add(HarmonicAngle{{angle_K_eV()}, {kAngleTheta0}})
        .add(MultiHarmonicDihedral{{dihedral_Cn_eV()}});
    sim.add(cpu::md::VirialStress{});
    sim.add(cpu::md::Barostat{P_target, tau_P, dt});
    sim.add(cpu::md::Deform{strain_rate, dt, /*axis=*/0});
    sim.build();

    const double box_x_initial = frame.box(0, 0);
    sim.run(frame, 1000);

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // No NaN explosion under the full bonded + LJ + halo + thermostat +
    // barostat + deform stack.
    for (int i = 0; i < frame.atoms.size(); ++i) {
        ASSERT_FALSE(std::isnan(frame.atoms.x[i]))
            << "rank " << my_rank << " atom " << i << " x went NaN";
    }

    // Box_x reflects the deformation pull modulated by the Berendsen
    // barostat trying to maintain P_target = 0 in the transverse axes.
    // We only assert it grew (Deform's pull is applied first; barostat
    // can damp but not reverse). Quantitative stress-strain analysis
    // belongs to msu-pe-validation §D against a LAMMPS reference.
    EXPECT_GT(frame.box(0, 0), box_x_initial)
        << "Deform pull did not grow box_x (initial=" << box_x_initial
        << ", final=" << frame.box(0, 0) << ")";

    // Box must agree across ranks.
    double bx = frame.box(0, 0);
    double bx_min, bx_max;
    MPI_Allreduce(&bx, &bx_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&bx, &bx_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_DOUBLE_EQ(bx_min, bx_max)
        << "ranks disagree on box_x after combined deform + barostat";

    // Total energy is finite.
    EXPECT_FALSE(std::isnan(frame.energy + frame.kinetic))
        << "global E went NaN after full pipeline";
}

int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    ::testing::InitGoogleTest(&argc, argv);
    const int rc = RUN_ALL_TESTS();
    MPI_Finalize();
    return rc;
}
