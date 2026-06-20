/// @file tests/cpu/gpw/test_cpu_gpw_diis.cpp
/// @brief Pulay/DIIS acceleration for the GPW SCF (gpw-edl-aimd, M5 prereq).
///
/// Two guarantees:
///   1. DIIS and the legacy linear density mixing converge to the SAME fixed
///      point (energy identical), with DIIS taking no more iterations — on a
///      cheap H2 system (fast both ways).
///   2. On the Zn semicore-3d system that stalls linear mixing (~157 iters at
///      mix=0.1), DIIS converges in far fewer iterations to the CP2K energy.

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

// DIIS must reach the SAME self-consistent energy as linear mixing, in no more
// iterations. Cheap H2 (2 electrons) so both paths are fast.
TEST(GpwDiis, MatchesLinearMixingOnH2) {
    const auto basis = parse_gth_basis(kData + "/BASIS_SET", "H", "DZV-GTH-PADE");
    const auto pot = parse_gth_potential(kData + "/POTENTIAL", "H", "GTH-PADE-q1");
    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.5, 2.0, 2.0};
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};

    GpwConfig base;
    base.cell[0] = 8.0;
    base.cell[1] = 4.0;
    base.cell[2] = 4.0;
    base.cutoff_ry = 140.0;
    base.n_occ = 1;
    base.max_iter = 200;
    base.tol = 1e-8;
    base.mix = 0.3;

    GpwConfig cfg_lin = base;
    cfg_lin.diis_depth = 0; // legacy linear density mixing
    GpwConfig cfg_diis = base;
    cfg_diis.diis_depth = 6; // Pulay DIIS

    GpwResult lin = gpw_scf(atoms, cfg_lin);
    GpwResult dii = gpw_scf(atoms, cfg_diis);

    std::printf("[GPW DIIS H2] linear: iters=%d E=%.10f | diis: iters=%d E=%.10f\n",
                lin.iters, lin.e_total, dii.iters, dii.e_total);

    ASSERT_TRUE(lin.converged);
    ASSERT_TRUE(dii.converged);
    // Same fixed point.
    EXPECT_LT(std::abs(dii.e_total - lin.e_total), 1e-7)
        << "linear=" << lin.e_total << " diis=" << dii.e_total;
    // DIIS is never slower.
    EXPECT_LE(dii.iters, lin.iters);
}

// The motivating case: Zn semicore-3d (GTH-PBE-q12) stalls linear mixing
// (~157 iters at mix=0.1). DIIS converges in far fewer iterations to the same
// CP2K energy. (ZnO diatomic, identical geometry to test_cpu_gpw_zno_cp2k.)
TEST(GpwDiis, ZnOSemicoreConvergesFast) {
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
    cfg.n_occ = 9;
    cfg.max_iter = 200;
    cfg.tol = 1e-7;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = false;
    cfg.diis_depth = 6; // DIIS (default); no mix=0.1 workaround needed

    GpwResult r = gpw_scf(atoms, cfg);

    std::printf("[GPW DIIS ZnO] iters=%d conv=%d  E_total=%.10f  (CP2K -76.4090567)\n",
                r.iters, (int)r.converged, r.e_total);

    EXPECT_TRUE(r.converged) << "iters=" << r.iters;
    // The linear-mixing baseline needed ~157 iters; DIIS must beat it decisively.
    EXPECT_LT(r.iters, 50) << "DIIS should converge the Zn semicore in << 157 iters";
    const double ref = -76.40905670322425; // CP2K 2026.1 (same as zno_cp2k test)
    EXPECT_LT(std::abs(r.e_total - ref), 5.0e-3) << "E_total=" << r.e_total;
}
