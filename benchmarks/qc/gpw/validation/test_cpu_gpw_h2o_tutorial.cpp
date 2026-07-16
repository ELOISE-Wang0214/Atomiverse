/// @file tests/cpu/gpw/test_cpu_gpw_h2o_tutorial.cpp
/// @brief GPW M3 — the getting-started h2o.inp tutorial: isolated H2O, PBE,
///        DZVP-MOLOPT-GTH + GTH-PBE-q6/q1, MT Poisson. The headline goal.
///
/// docs/getting-started/h2o.inp: O at (5,5,5) Ang, H at (5.7586,5,5.5043) and
/// (4.2414,5,5.5043) Ang, 10 Ang cubic box, CUTOFF 400, PERIODIC NONE / PSOLVER
/// MT, PBE. 8 valence electrons (n_occ=4). CP2K total energy -17.21233714458516
/// (decomp: self=-44.34715500 core_ham=13.07138998 hartree=18.34339895
/// xc=-4.27997143). Exercises the full non-periodic GGA path: MOLOPT basis with
/// spherical d, GTH-PBE local+nonlocal, GGA-on-grid, and MT isolated Poisson.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"

using namespace ATV::cpu::gpw;

namespace {
constexpr double kAng = 1.8897259885789; // Angstrom -> Bohr
} // namespace

TEST(GpwH2OTutorial, IsolatedPbeWater) {
    const auto basis_O = GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-GTH");
    const auto basis_H = GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH");
    const auto pot_O = GthPotential(molrs::Element::O, "GTH-PBE-q6");
    const auto pot_H = GthPotential(molrs::Element::H, "GTH-PBE-q1");

    const double O[3] = {5.0 * kAng, 5.0 * kAng, 5.0 * kAng};
    const double H1[3] = {5.7586 * kAng, 5.0 * kAng, 5.5043 * kAng};
    const double H2[3] = {4.2414 * kAng, 5.0 * kAng, 5.5043 * kAng};

    std::vector<GpwAtom> atoms = {GpwAtom(O, basis_O, pot_O),
                                  GpwAtom(H1, basis_H, pot_H),
                                  GpwAtom(H2, basis_H, pot_H)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 10.0 * kAng;
    cfg.cutoff_ry = 400.0;
    cfg.n_occ = 4; // 8 valence electrons
    cfg.max_iter = 200;
    cfg.tol = 1e-6;
    cfg.mix = 0.3;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = true; // MT isolated Poisson

    GPWSolver solver{cfg};
    solver.build(atoms);
    const GpwResult r = solver.run(atoms);

    std::printf(
        "[GPW H2O] iters=%d conv=%d  E_total=%.10f  (CP2K -17.2123371446)\n"
        "          kinetic=%.8f TrP*Vpp=%.8f core_ham=%.8f es=%.8f xc=%.8f\n",
        r.iters, (int)r.converged, r.e_total, r.e_kinetic, r.e_loc_sr,
        r.e_kinetic + r.e_loc_sr, r.e_es, r.e_xc);

    EXPECT_TRUE(r.converged) << "iters=" << r.iters;
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_LT(r.e_xc, 0.0);
    const double ref = -17.21233714458516;
    EXPECT_LT(std::abs(r.e_total - ref), 2.0e-3)
        << "E_total=" << r.e_total << " |dE|=" << std::abs(r.e_total - ref);
}
