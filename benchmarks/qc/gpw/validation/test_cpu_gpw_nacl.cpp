/// @file tests/cpu/gpw/test_cpu_gpw_nacl_cp2k.cpp
/// @brief GPW vs CP2K parity for a periodic NaCl diatomic — validates the last
///        two heaviest supported elements: Na GTH-PBE-q9 (semicore, 9 valence e-) and
///        Cl GTH-PBE-q7. Rung 2 of the GPW per-element anchor ladder.
///
/// System: Na at (5,5,5) Å, Cl at (5,5,7.36) Å (2.36 Å bond) in a 10 Å cubic
/// PERIODIC XYZ box; PBE; CUTOFF 600 Ry; GTH-PBE-q9 (Na) / GTH-PBE-q7 (Cl);
/// DZVP-MOLOPT-SR-GTH. 16 valence electrons -> n_occ = 8 (closed shell).
///
/// CP2K 2026.1 reference (local run, same geometry/cell/CUTOFF/PP/basis/XC,
/// PERIODIC XYZ, EPS_SCF 1e-7, converged 25 steps):
/// Total energy = -62.838137517556305 Ha.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"

using ATV::cpu::gpw::GaussianBasis;
using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwConfig;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::GPWSolver;
using ATV::cpu::gpw::GthPotential;
using ATV::cpu::gpw::XcKind;

namespace {
constexpr double kAng = 1.8897259885789; // Angstrom -> Bohr
} // namespace

TEST(GpwNaClCP2K, PeriodicDiatomicPbe) {
    const auto basis_Na =
        GaussianBasis(molrs::Element::Na, "DZVP-MOLOPT-SR-GTH");
    const auto basis_Cl =
        GaussianBasis(molrs::Element::Cl, "DZVP-MOLOPT-SR-GTH");
    const auto pot_Na = GthPotential(molrs::Element::Na, "GTH-PBE-q9");
    const auto pot_Cl = GthPotential(molrs::Element::Cl, "GTH-PBE-q7");

    const double Na[3] = {5.0 * kAng, 5.0 * kAng, 5.0 * kAng};
    const double Cl[3] = {5.0 * kAng, 5.0 * kAng, 7.36 * kAng};

    std::vector<GpwAtom> atoms = {GpwAtom(Na, basis_Na, pot_Na),
                                  GpwAtom(Cl, basis_Cl, pot_Cl)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 10.0 * kAng;
    cfg.cutoff_ry = 600.0;
    cfg.n_occ = 8; // 16 valence electrons (Na q9 + Cl q7)
    cfg.max_iter = 400;
    cfg.tol = 1e-7;
    // Na semicore (2s2p included in q9) is stiff like Zn; linear mixing at 0.1
    // converges cleanly. (DIIS/Pulay would relax this — GPW solver TODO.)
    cfg.mix = 0.1;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = false; // periodic Poisson (PERIODIC XYZ)

    GPWSolver solver{cfg};
    solver.build(atoms);
    const GpwResult r = solver.run(atoms);

    const double ref = -62.838137517556305; // CP2K 2026.1
    std::printf(
        "[GPW NaCl] iters=%d conv=%d  E_total=%.10f  (CP2K %.10f)\n"
        "           kinetic=%.8f loc_sr=%.8f nonlocal=%.8f es=%.8f xc=%.8f\n"
        "           |dE|=%.3e Ha\n",
        r.iters, (int)r.converged, r.e_total, ref, r.e_kinetic, r.e_loc_sr,
        r.e_nonlocal, r.e_es, r.e_xc, std::abs(r.e_total - ref));

    EXPECT_TRUE(r.converged) << "iters=" << r.iters;
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_LT(r.e_xc, 0.0);
    // Single uniform grid vs CP2K CUTOFF/REL_CUTOFF multigrid leaves a
    // discretization residual; gate at the same 5 mHa scale as the sibling rungs.
    const double err = std::abs(r.e_total - ref);
    EXPECT_LT(err, 5.0e-3) << "E_total=" << r.e_total << " ref=" << ref
                           << " |dE|=" << err;
}
