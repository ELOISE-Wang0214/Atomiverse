/// @file tests/cpu/gpw/test_cpu_gpw_multigrid_scf.cpp
/// @brief End-to-end CP2K REL_CUTOFF multigrid SCF: the product-based multigrid
///        density / KS path (ngrids>1, rel_cutoff>0) must reproduce the
///        single-grid GPW energy to grid-transfer accuracy.
///
/// The single-grid path collocates the whole density on one fine grid; the
/// multigrid path maps each primitive-pair product (by its combined exponent) to
/// the coarsest grid level that resolves it, collocates there, prolongs to the
/// fine grid for Hartree+XC, and restricts the fine potential back to integrate
/// the KS matrix. Adjoint transfer (prolong/restrict, collocate/integrate) keeps
/// the SCF variational, so the total energy agrees with the single grid to the
/// trilinear grid-transfer residual (a few 1e-4 Ha) — not bit-identical.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"

using namespace ATV::cpu::gpw;

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;

// Periodic water (DZVP-MOLOPT-SR-GTH / GTH-PBE) in a ~5 Å cubic box — the EDL
// water basis, Γ-point images active. Closed shell, converges cleanly.
std::vector<GpwAtom> make_water() {
    constexpr double kAng = 1.8897259885789; // Å -> Bohr
    const auto bO = parse_gth_basis(kData + "/BASIS_MOLOPT", "O", "DZVP-MOLOPT-SR-GTH");
    const auto bH = parse_gth_basis(kData + "/BASIS_MOLOPT", "H", "DZVP-MOLOPT-SR-GTH");
    const auto pO = parse_gth_potential(kData + "/GTH_POTENTIALS", "O", "GTH-PBE-q6");
    const auto pH = parse_gth_potential(kData + "/GTH_POTENTIALS", "H", "GTH-PBE-q1");
    const double rO[3] = {2.6 * kAng, 2.5 * kAng, 2.4 * kAng};
    const double rH1[3] = {3.16 * kAng, 2.5 * kAng, 2.9 * kAng};
    const double rH2[3] = {2.1 * kAng, 2.5 * kAng, 2.95 * kAng};
    return {make_gpw_atom(rO, bO, pO), make_gpw_atom(rH1, bH, pH),
            make_gpw_atom(rH2, bH, pH)};
}

GpwConfig base_config() {
    constexpr double kAng = 1.8897259885789;
    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 5.0 * kAng;
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 4; // 8 valence electrons (O q6 + 2 H q1)
    cfg.max_iter = 200;
    cfg.tol = 1e-9;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = false; // PERIODIC: exercises the lattice-sum collocation
    return cfg;
}
} // namespace

// The product-based multigrid (ngrids=4, rel_cutoff=60) reproduces the
// single-grid (ngrids=1) total energy to grid-transfer accuracy.
TEST(GpwMultigridScf, FourLevelsMatchSingleGrid) {
    const std::vector<GpwAtom> atoms = make_water();

    GpwConfig cfg1 = base_config();
    cfg1.ngrids = 1;
    const GpwResult r1 = gpw_scf(atoms, cfg1);

    GpwConfig cfg4 = base_config();
    cfg4.ngrids = 4;
    cfg4.rel_cutoff = 60.0;
    const GpwResult r4 = gpw_scf(atoms, cfg4);

    std::printf("[GPW mg-scf] ngrids=1 E=%.8f (it=%d conv=%d)\n", r1.e_total,
                r1.iters, r1.converged);
    std::printf("[GPW mg-scf] ngrids=4 E=%.8f (it=%d conv=%d)\n", r4.e_total,
                r4.iters, r4.converged);
    std::printf("[GPW mg-scf] dE = %.3e Ha\n", std::abs(r4.e_total - r1.e_total));
    std::printf("[GPW mg-scf] g1: T=%.5f loc=%.5f es=%.5f xc=%.5f\n", r1.e_kinetic,
                r1.e_loc_sr, r1.e_es, r1.e_xc);
    std::printf("[GPW mg-scf] g4: T=%.5f loc=%.5f es=%.5f xc=%.5f\n", r4.e_kinetic,
                r4.e_loc_sr, r4.e_es, r4.e_xc);

    EXPECT_TRUE(r1.converged);
    EXPECT_TRUE(r4.converged);
    // Grid-transfer (trilinear) accuracy: gate at 1e-3 Ha.
    EXPECT_LT(std::abs(r4.e_total - r1.e_total), 1e-3)
        << "multigrid SCF energy diverges from the single-grid reference";
}
