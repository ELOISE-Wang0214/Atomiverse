/// @file tests/cpu/gpw/test_cpu_gpw_zno_cp2k.cpp
/// @brief GPW vs CP2K parity for a periodic ZnO diatomic — the first
///        transition-metal validation (Zn GTH-PBE-q12 d-channel projector +
///        DZVP-MOLOPT-SR-GTH, which carries f-shells, lmax=3). Rung 1 of the
///        gpw-edl-aimd per-element ladder (ac-004).
///
/// System: Zn at (4,4,4) Å, O at (4,4,5.72) Å (1.72 Å bond) in an 8 Å cubic
/// PERIODIC XYZ box; PBE; CUTOFF 600 Ry; GTH-PBE-q12 (Zn) / GTH-PBE-q6 (O);
/// DZVP-MOLOPT-SR-GTH. 18 valence electrons -> n_occ = 9 (closed shell).
///
/// CP2K 2026.1 reference (local run, same geometry/cell/CUTOFF/PP/basis/XC,
/// PERIODIC XYZ, EPS_SCF 1e-7): Total energy = -76.40905670322425 Ha.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwConfig;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::gpw_scf;
using ATV::cpu::gpw::make_gpw_atom;
using ATV::cpu::gpw::parse_gth_basis;
using ATV::cpu::gpw::parse_gth_potential;
using ATV::cpu::gpw::XcKind;

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;
constexpr double kAng = 1.8897259885789; // Angstrom -> Bohr
} // namespace

TEST(GpwZnOCP2K, PeriodicDiatomicPbe) {
    const auto basis_Zn =
        parse_gth_basis(kData + "/BASIS_MOLOPT", "Zn", "DZVP-MOLOPT-SR-GTH");
    const auto basis_O =
        parse_gth_basis(kData + "/BASIS_MOLOPT", "O", "DZVP-MOLOPT-SR-GTH");
    const auto pot_Zn =
        parse_gth_potential(kData + "/GTH_POTENTIALS", "Zn", "GTH-PBE-q12");
    const auto pot_O =
        parse_gth_potential(kData + "/GTH_POTENTIALS", "O", "GTH-PBE-q6");

    const double Zn[3] = {4.0 * kAng, 4.0 * kAng, 4.0 * kAng};
    const double O[3] = {4.0 * kAng, 4.0 * kAng, 5.72 * kAng};

    std::vector<GpwAtom> atoms = {make_gpw_atom(Zn, basis_Zn, pot_Zn),
                                  make_gpw_atom(O, basis_O, pot_O)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 8.0 * kAng;
    cfg.cutoff_ry = 600.0;
    cfg.n_occ = 9; // 18 valence electrons (Zn q12 + O q6)
    cfg.max_iter = 400;
    cfg.tol = 1e-7;
    // Zn is an 18-valence-electron semicore-3d transition metal; plain linear
    // mixing at 0.3 oscillates, at 0.1 it converges cleanly (~157 iters) to
    // 0.1 mHa of CP2K. (A DIIS/Pulay mixer would relax this — GPW solver TODO.)
    cfg.mix = 0.1;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = false; // periodic Poisson (PERIODIC XYZ)

    GpwResult r = gpw_scf(atoms, cfg);

    const double ref = -76.40905670322425; // CP2K 2026.1
    std::printf("[GPW ZnO] iters=%d conv=%d  E_total=%.10f  (CP2K %.10f)\n"
                "          kinetic=%.8f loc_sr=%.8f nonlocal=%.8f es=%.8f xc=%.8f\n"
                "          |dE|=%.3e Ha\n",
                r.iters, (int)r.converged, r.e_total, ref, r.e_kinetic,
                r.e_loc_sr, r.e_nonlocal, r.e_es, r.e_xc,
                std::abs(r.e_total - ref));

    EXPECT_TRUE(r.converged) << "iters=" << r.iters;
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_LT(r.e_xc, 0.0);
    // Single uniform grid vs CP2K CUTOFF/REL_CUTOFF multigrid leaves a
    // discretization residual (the H2/H2O rungs land ~1.4-2 mHa); gate at the
    // same scale used by the sibling CP2K-parity tests.
    const double err = std::abs(r.e_total - ref);
    EXPECT_LT(err, 5.0e-3) << "E_total=" << r.e_total << " ref=" << ref
                           << " |dE|=" << err;
}
