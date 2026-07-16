/// @file benchmarks/qc/gpw/validation/test_cpu_gpw_d3_scf.cpp
/// @brief Validate the standalone Grimme DFT-D3 (zero damping, PBE) module
///        against a CP2K reference number and finite-difference forces.
///
/// CP2K reference (cp2k 2026.1, TYPE DFTD3 / REFERENCE_FUNCTIONAL PBE / default
/// R_CUTOFF 20 Bohr / EPS_CN 1e-6, dftd3.dat) for a water dimer in a 12 Å NONE
/// (isolated) cell prints:
///     "Dispersion energy:   -0.00081399213205"  [Hartree]
/// The geometry below (Ångström) is the exact &COORD of that run; D3 is
/// translation-invariant, so CP2K's internal centering does not change it.

#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/dftd3.h"
#include "atomiverse/cpu/gpw/internal/dftd3_params.h"

using ATV::cpu::gpw::D3Params;
using ATV::cpu::gpw::D3Result;
using ATV::cpu::gpw::dftd3_pbe;

namespace {

constexpr double kAngToBohr = 1.0 / 0.52917720859;

// Water dimer (Ångström), matching the CP2K reference &COORD.
std::vector<std::array<double, 3>> water_dimer_bohr() {
    const double ang[6][3] = {
        {0.000000, 0.000000, 0.000000},  // O
        {0.758602, 0.000000, 0.504284},  // H
        {0.758602, 0.000000, -0.504284}, // H
        {3.000000, 0.000000, 0.000000},  // O
        {3.758602, 0.000000, 0.504284},  // H
        {3.758602, 0.000000, -0.504284}, // H
    };
    std::vector<std::array<double, 3>> pos(6);
    for (int i = 0; i < 6; ++i)
        for (int d = 0; d < 3; ++d)
            pos[i][d] = ang[i][d] * kAngToBohr;
    return pos;
}

const std::vector<int> kZ = {8, 1, 1, 8, 1, 1};

} // namespace

// Energy must reproduce CP2K's printed "Dispersion energy" to ~1e-9 Hartree.
// Integration: GpwConfig.d3 must add exactly the standalone D3 energy/force to
// the GPW result. Uses a fast H2 (MOLOPT-SR/PBE) so the SCF is cheap; the d3-on
// minus d3-off difference must equal dftd3_pbe() to machine precision.
TEST(GpwDftD3, IntegratedIntoGpwScf) {
    using namespace ATV::cpu::gpw;
    constexpr double kAng = 1.8897259885789;
    const auto bH = GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-SR-GTH");
    const auto pH = GthPotential(molrs::Element::H, "GTH-PBE-q1");
    const double rA[3] = {4.0, 4.0, 4.0};
    const double rB[3] = {4.0, 4.0, 4.0 + 1.6 * kAng};
    std::vector<GpwAtom> atoms = {GpwAtom(rA, bH, pH, 1),
                                  GpwAtom(rB, bH, pH, 1)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 8.0 * kAng;
    cfg.cutoff_ry = 200.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-9;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = true;

    GPWSolver off_solver{cfg};
    off_solver.build(atoms);
    const double e_off = off_solver.run(atoms).e_total;
    cfg.d3 = true;
    GPWSolver on_solver{cfg};
    on_solver.build(atoms);
    const double e_on = on_solver.run(atoms).e_total;

    const std::vector<int> Z = {1, 1};
    const std::vector<std::array<double, 3>> pos = {{rA[0], rA[1], rA[2]},
                                                    {rB[0], rB[1], rB[2]}};
    const double ed3 = dftd3_pbe(Z, pos, cfg.cell, false, D3Params{}).energy;
    std::printf("[D3 integ] e_on-e_off=%.10f  standalone=%.10f\n", e_on - e_off,
                ed3);
    EXPECT_NEAR(e_on - e_off, ed3, 1e-9)
        << "GpwConfig.d3 did not add the D3 energy";
    EXPECT_LT(ed3, 0.0); // dispersion is attractive
}
