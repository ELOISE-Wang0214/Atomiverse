/// @file tests/mpi/test_gpw_kpoint.cpp
/// @brief MPI k-point parallel GPW RKS SCF — serial/parallel parity.
///
/// Launched under `mpirun -n N` (N = 1, 2, 4 per CTest config). The
/// k-point distribution is round-robin over MPI_COMM_WORLD inside
/// `cpu::gpw::GpwDft::run()` multi-k path; density and energy scalars
/// are Allreduced, so every rank must converge to the identical energy.

#include <gtest/gtest.h>

#include <mpi.h>

#include <cmath>
#include <vector>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/molrs.hpp"

using namespace ATV;
using namespace ATV::cpu::gpw;

namespace {

constexpr double kA2B = 1.0 / 0.52917720859;

std::vector<GpwAtom> si_bct(double a_bohr) {
    const double at = a_bohr / std::sqrt(2.0);
    const GaussianBasis bSi(molrs::Element::Si, "DZVP-MOLOPT-GTH");
    const GthPotential pSi(molrs::Element::Si, "GTH-PBE-q4");
    const double frac[4][3] = {
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {0.5, 0.0, 0.25}, {0.0, 0.5, 0.75}};
    std::vector<GpwAtom> atoms;
    for (int i = 0; i < 4; ++i) {
        const double pos[3] = {frac[i][0] * at, frac[i][1] * at,
                               frac[i][2] * a_bohr};
        atoms.emplace_back(pos, bSi, pSi, 14);
    }
    return atoms;
}

GpwDft make_method(double a_bohr) {
    GpwDft m;
    m.cell[0] = m.cell[1] = a_bohr / std::sqrt(2.0);
    m.cell[2] = a_bohr;
    m.cutoff_ry = 60.0; // unit-test budget (seconds)
    m.n_occ = 8;
    m.max_iter = 30;
    m.tol = 1e-6;
    m.xc = XcKind::PBE;
    m.eps_default = 1e-6;
    m.k_mesh[0] = m.k_mesh[1] = m.k_mesh[2] = 2;
    m.k_time_reversal = true;
    return m;
}

} // namespace

TEST(GpwMpiKpoint, AllRanksAgreeOnEnergy) {
    const double a = 5.43 * kA2B;
    auto atoms = si_bct(a);
    GpwDft method = make_method(a);
    method.build(atoms);
    const auto& r = method.run(atoms);

    ASSERT_TRUE(std::isfinite(r.e_total));
    EXPECT_EQ(r.n_kpts, 4);

    int rank = 0, nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int expected_local = 0;
    for (int ik = rank; ik < r.n_kpts; ik += nprocs)
        ++expected_local;
    EXPECT_EQ(r.n_kpts_local, expected_local);

    double local = r.e_total;
    double min_e = 0, max_e = 0;
    MPI_Allreduce(&local, &min_e, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local, &max_e, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_NEAR(min_e, max_e, 1e-12)
        << "ranks diverged: density/energy Allreduce broken";

    int local_iters = r.iters;
    int min_it = 0, max_it = 0;
    MPI_Allreduce(&local_iters, &min_it, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_iters, &max_it, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_EQ(min_it, max_it) << "ranks took different SCF iteration counts";

    int sum_local = 0;
    MPI_Allreduce(&r.n_kpts_local, &sum_local, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    EXPECT_EQ(sum_local, r.n_kpts);
}

TEST(GpwMpiKpoint, FiniteConvergedSmoke) {
    const double a = 5.43 * kA2B;
    auto atoms = si_bct(a);
    GpwDft method = make_method(a);
    method.build(atoms);
    const auto& r = method.run(atoms);
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_TRUE(r.converged);
    EXPECT_GT(r.iters, 0);
}

int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    ::testing::InitGoogleTest(&argc, argv);
    const int rc = RUN_ALL_TESTS();
    MPI_Finalize();
    return rc;
}
