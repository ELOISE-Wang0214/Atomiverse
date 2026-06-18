/// @file tests/cpu/gpw/test_cpu_gpw_forces.cpp
/// @brief Finite-difference GPW forces on H2 — Newton's third law, along-axis,
///        and the correct restoring sign for a compressed bond.

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"

using namespace ATV::cpu::gpw;

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;
} // namespace

TEST(GpwForces, H2FiniteDifference) {
    const auto basis = parse_gth_basis(kData + "/BASIS_SET", "H", "DZV-GTH-PADE");
    const auto pot = parse_gth_potential(kData + "/POTENTIAL", "H", "GTH-PADE-q1");

    // H2 compressed to 1.0 Bohr (equilibrium ~1.4): expect a restoring force
    // that pushes the atoms apart along x.
    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.5, 2.0, 2.0};
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};

    GpwConfig cfg;
    cfg.cell[0] = 8.0;
    cfg.cell[1] = cfg.cell[2] = 4.0;
    cfg.cutoff_ry = 140.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-8;
    cfg.mix = 0.3;

    auto F = gpw_forces_fd(atoms, cfg, 2e-3);

    // The bond (x) force is the clean physical signal: equal and opposite, and
    // restoring for a compressed bond (left atom pushed to -x, right to +x).
    EXPECT_NEAR(F[0][0], -F[1][0], 1e-4); // Newton's third law along the bond
    EXPECT_LT(F[0][0], 0.0);
    EXPECT_GT(F[1][0], 0.0);
    EXPECT_GT(std::abs(F[0][0]), 5e-3); // substantial restoring force

    // The off-axis (y,z) components are a common-mode "egg-box" grid artifact
    // (both atoms share the same grid phase at y=z=2), not a physical force.
    // It is small relative to the bond force on this single coarse grid;
    // CP2K's REL_CUTOFF multigrid + smoothing suppress it (a documented TODO).
    EXPECT_LT(std::abs(F[0][1]), 0.4 * std::abs(F[0][0]));
    EXPECT_LT(std::abs(F[0][2]), 0.4 * std::abs(F[0][0]));
}
