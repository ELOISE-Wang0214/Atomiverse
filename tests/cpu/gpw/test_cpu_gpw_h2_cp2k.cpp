/// @file tests/cpu/gpw/test_cpu_gpw_h2_cp2k.cpp
/// @brief GPW M1 (in progress): real CP2K H2 — DZV-GTH-PADE + GTH-PADE-q1, PADE
///        XC — toward the CP2K reference Total energy -1.06345574964568 Ha.
///
/// Mirrors tests/QS/regtest-gpw-1/H2.inp: H at (0,0,0) and (1,0,0) Bohr in an
/// 8x4x4 box, CUTOFF 140 Ry. The CP2K energy decomposition (from a local run) is
///   self=-2.82094792  core_ham=0.80557972  hartree=1.51921160  xc=-0.56729915
/// This test reports how close the assembled GPW energy is; the tolerance is
/// loose while the electrostatic/basis conventions are tuned to bit-parity.

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

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;
} // namespace

TEST(GpwH2CP2K, RealBasisConverges) {
    const auto basis = parse_gth_basis(kData + "/BASIS_SET", "H", "DZV-GTH-PADE");
    const auto pot = parse_gth_potential(kData + "/POTENTIAL", "H", "GTH-PADE-q1");

    // CP2K H2.inp geometry, in the 8x4x4 box (place the molecule away from the
    // origin so it does not straddle the periodic boundary).
    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.5, 2.0, 2.0}; // 1.0 Bohr apart
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};

    GpwConfig cfg;
    cfg.cell[0] = 8.0;
    cfg.cell[1] = 4.0;
    cfg.cell[2] = 4.0;
    cfg.cutoff_ry = 140.0; // CP2K CUTOFF 140
    cfg.n_occ = 1;         // 2 electrons
    cfg.max_iter = 200;
    cfg.tol = 1e-7;
    cfg.mix = 0.3;

    GpwResult r = gpw_scf(atoms, cfg);

    std::printf("[GPW H2] iters=%d conv=%d  E_total=%.10f  (CP2K -1.0634557496)\n"
                "         kinetic=%.8f  TrP*Vpp=%.8f  core_ham=%.8f  es=%.8f  xc=%.8f\n",
                r.iters, (int)r.converged, r.e_total, r.e_kinetic, r.e_loc_sr,
                r.e_kinetic + r.e_loc_sr, r.e_es, r.e_xc);

    // Per-H DZV-GTH-PADE = 2 contracted s functions -> 4 basis functions total.
    EXPECT_TRUE(r.converged) << "iters=" << r.iters;
    EXPECT_TRUE(std::isfinite(r.e_total));
    // Report distance to the CP2K reference (loose gate during convention tuning).
    // First real-data run lands ~1.4 mHa from CP2K; the residual is grid
    // discretization (single uniform grid vs CP2K's CUTOFF/REL_CUTOFF multigrid).
    // Gate at 3 mHa to lock that in; tightening to 1e-6 needs grid-scheme parity.
    const double ref = -1.06345574964568;
    const double err = std::abs(r.e_total - ref);
    EXPECT_LT(err, 3.0e-3) << "E_total=" << r.e_total << " ref=" << ref
                           << " |dE|=" << err;
}
