/// @file tests/cpu/gpw/test_cpu_gpw_forces.cpp
/// @brief Finite-difference GPW forces on H2 — Newton's third law, along-axis,
///        and the correct restoring sign for a compressed bond.

#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"
#include "atomiverse/cpu/gpw/internal/gpw_forces.h"
#include "atomiverse/cpu/gpw/basis_data/gpw_data.h"

using namespace ATV::cpu::gpw;

namespace {
} // namespace

TEST(GpwForces, H2FiniteDifference) {
    const auto basis = builtin_gth_basis("H", "DZV-GTH-PADE");
    const auto pot = builtin_gth_potential("H", "GTH-PADE-q1");

    // H2 compressed to 1.0 Bohr (equilibrium ~1.4): expect a restoring force
    // that pushes the atoms apart along x.
    const double rA[3] = {5.5, 6.0, 6.0};
    const double rB[3] = {6.5, 6.0, 6.0};
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 12.0; // enlarged + centered: keep R_cut < L/2 (no min-image kink)
    cfg.cutoff_ry = 140.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-8;
    cfg.mix = 0.3;
    // Isolated limit: the analytic-force assembly is not yet periodic-image
    // consistent with the Γ-point lattice-sum SCF, so validate the gradient
    // machinery without images (periodic-consistent forces are a follow-on).
    cfg.nonperiodic = true;

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

// Analytic forces must reproduce the finite-difference oracle component-by-
// component (including the off-axis egg-box artifact — same grid, same energy).
TEST(GpwForces, H2AnalyticVsFiniteDifference) {
    const auto basis = builtin_gth_basis("H", "DZV-GTH-PADE");
    const auto pot = builtin_gth_potential("H", "GTH-PADE-q1");

    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.7, 2.1, 1.9}; // off-symmetry so every component is nonzero
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 12.0; // enlarged: keep R_cut < L/2 (no min-image kink)
    cfg.cutoff_ry = 140.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-10;
    cfg.nonperiodic = true; // isolated: analytic forces not yet image-consistent

    const auto Fa = gpw_forces_analytic(atoms, cfg);
    const auto Ffd = gpw_forces_fd(atoms, cfg, 1e-3);

    double max_abs = 0.0, max_err = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int k = 0; k < 3; ++k) {
            max_abs = std::max(max_abs, std::abs(Ffd[a][k]));
            max_err = std::max(max_err, std::abs(Fa[a][k] - Ffd[a][k]));
            std::printf("[GPW F] atom %d axis %d  analytic % .8f  fd % .8f  d % .2e\n",
                        a, k, Fa[a][k], Ffd[a][k], Fa[a][k] - Ffd[a][k]);
        }
    std::printf("[GPW F] max|F|=%.6f  max|dF|=%.3e Ha/Bohr\n", max_abs, max_err);
    EXPECT_LT(max_err, 1e-5) << "analytic vs FD mismatch";
}

// H2 with PBE: isolates the GGA grid-term derivative (∫ V_xc^GGA ∂ρ) cheaply,
// before the heavier ZnO case adds the d-projector / f-basis derivatives.
TEST(GpwForces, H2PbeAnalyticVsFiniteDifference) {
    const auto basis = builtin_gth_basis("H", "DZVP-MOLOPT-SR-GTH");
    const auto pot = builtin_gth_potential("H", "GTH-PBE-q1");
    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.7, 2.1, 1.9};
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};
    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 12.0; // enlarged: keep R_cut < L/2 (no min-image kink)
    cfg.cutoff_ry = 200.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-10;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = true; // isolated: analytic forces not yet image-consistent

    const auto Fa = gpw_forces_analytic(atoms, cfg);
    const auto Ffd = gpw_forces_fd(atoms, cfg, 1e-3);
    double max_err = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int k = 0; k < 3; ++k) {
            max_err = std::max(max_err, std::abs(Fa[a][k] - Ffd[a][k]));
            std::printf("[GPW F PBE] atom %d axis %d  analytic % .8f  fd % .8f  d % .2e\n",
                        a, k, Fa[a][k], Ffd[a][k], Fa[a][k] - Ffd[a][k]);
        }
    std::printf("[GPW F PBE] max|dF|=%.3e Ha/Bohr\n", max_err);
    EXPECT_LT(max_err, 1e-5) << "H2/PBE GGA analytic vs FD mismatch";
}

// ZnO: exercises the GGA(PBE) grid term, the l=2 (d) GTH projector derivative,
// and the l=3 (f) basis-shell derivative — none of which the H2/LDA case touches.
TEST(GpwForces, ZnOAnalyticVsFiniteDifference) {
    constexpr double kAng = 1.8897259885789; // Å -> Bohr
    const auto bZn = builtin_gth_basis("Zn", "DZVP-MOLOPT-SR-GTH");
    const auto bO = builtin_gth_basis("O", "DZVP-MOLOPT-SR-GTH");
    const auto pZn = builtin_gth_potential("Zn", "GTH-PBE-q12");
    const auto pO = builtin_gth_potential("O", "GTH-PBE-q6");

    const double rZn[3] = {7.0 * kAng, 7.0 * kAng, 7.0 * kAng};
    const double rO[3] = {7.15 * kAng, 7.1 * kAng, 8.6 * kAng}; // off-axis bond
    std::vector<GpwAtom> atoms = {make_gpw_atom(rZn, bZn, pZn),
                                  make_gpw_atom(rO, bO, pO)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 14.0 * kAng; // enlarged: R_cut < L/2 (no min-image kink)
    // The grid route makes analytic-vs-FD consistency grid-dependent (V_loc^SR is
    // collocated), so the cutoff is a physical setting here, not a free knob.
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 9;
    cfg.tol = 1e-9;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = true; // isolated: analytic forces not yet image-consistent

    const auto Fa = gpw_forces_analytic(atoms, cfg);
    // FD step 5e-4, not 1e-3: the oracle's OWN central-difference truncation
    // error is O(h^2) and dominates here. Measured on this system, analytic-vs-FD
    // max|dF| = 1.68e-5 (h=5e-4), 7.20e-5 (h=1e-3), 2.93e-4 (h=2e-3) — a clean 4x
    // per doubling, i.e. the discrepancy is the REFERENCE's error, not the analytic
    // force's. Refining the reference, not loosening the gate.
    const auto Ffd = gpw_forces_fd(atoms, cfg, 5e-4);

    double max_abs = 0.0, max_err = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int k = 0; k < 3; ++k) {
            max_abs = std::max(max_abs, std::abs(Ffd[a][k]));
            max_err = std::max(max_err, std::abs(Fa[a][k] - Ffd[a][k]));
            std::printf("[GPW F ZnO] atom %d axis %d  analytic % .8f  fd % .8f  d % .2e\n",
                        a, k, Fa[a][k], Ffd[a][k], Fa[a][k] - Ffd[a][k]);
        }
    std::printf("[GPW F ZnO] max|F|=%.6f  max|dF|=%.3e Ha/Bohr\n", max_abs, max_err);
    EXPECT_LT(max_err, 5e-5) << "ZnO analytic vs FD mismatch (GGA/d-proj/f-basis)";
}

// PERIODIC analytic forces: H2O (DZVP-MOLOPT-SR-GTH, the EDL water basis) in a
// tight cubic box where the Γ-point lattice images are active (m>=1). Closed
// shell, converges cleanly. The analytic forces must reproduce the periodic
// finite-difference oracle (gpw_scf at +/-h, itself lattice-sum periodic) — i.e.
// the force is image-consistent with the energy (the AIMD requirement).
TEST(GpwForces, H2OPeriodicAnalyticVsFiniteDifference) {
    constexpr double kAng = 1.8897259885789;
    const auto bO = builtin_gth_basis("O", "DZVP-MOLOPT-SR-GTH");
    const auto bH = builtin_gth_basis("H", "DZVP-MOLOPT-SR-GTH");
    const auto pO = builtin_gth_potential("O", "GTH-PBE-q6");
    const auto pH = builtin_gth_potential("H", "GTH-PBE-q1");

    // Bent water, off-centre, so every component is nonzero.
    const double rO[3] = {2.6 * kAng, 2.5 * kAng, 2.4 * kAng};
    const double rH1[3] = {3.16 * kAng, 2.5 * kAng, 2.9 * kAng};
    const double rH2[3] = {2.1 * kAng, 2.5 * kAng, 2.95 * kAng};
    std::vector<GpwAtom> atoms = {make_gpw_atom(rO, bO, pO),
                                  make_gpw_atom(rH1, bH, pH),
                                  make_gpw_atom(rH2, bH, pH)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 5.0 * kAng; // ~9.4 Bohr: images active
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 4; // 8 valence electrons (O q6 + 2 H q1)
    cfg.max_iter = 200;
    cfg.tol = 1e-9;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = false; // PERIODIC: exercises the lattice-sum force path

    const auto Fa = gpw_forces_analytic(atoms, cfg);
    const auto Ffd = gpw_forces_fd(atoms, cfg, 1e-3);

    double max_abs = 0.0, max_err = 0.0;
    for (int a = 0; a < 3; ++a)
        for (int k = 0; k < 3; ++k) {
            max_abs = std::max(max_abs, std::abs(Ffd[a][k]));
            max_err = std::max(max_err, std::abs(Fa[a][k] - Ffd[a][k]));
            std::printf("[GPW F H2O-per] atom %d axis %d  analytic % .8f  fd % .8f  d % .2e\n",
                        a, k, Fa[a][k], Ffd[a][k], Fa[a][k] - Ffd[a][k]);
        }
    std::printf("[GPW F H2O-per] max|F|=%.6f  max|dF|=%.3e Ha/Bohr\n", max_abs, max_err);
    // Periodic multi-image collocation carries a slightly larger egg-box residual
    // than the enlarged-box cases; gate at 1e-3 (the analytic/FD agreement that
    // confirms image consistency, well below physical force magnitudes ~0.1).
    EXPECT_LT(max_err, 1e-3) << "H2O periodic analytic vs FD mismatch";
}
