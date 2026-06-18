/// @file tests/cpu/gpw/test_cpu_gpw_solver.cpp
/// @brief End-to-end GPW SCF: the closed-shell loop converges for a minimal H2.
///
/// This exercises the whole GPW chain in one self-consistent loop — basis S/T,
/// GTH local PP, collocation, FFT-Poisson electrostatics (Gaussian core charge),
/// PADE XC, generalized eigensolve, density mixing. The assertion is that it
/// CONVERGES to a finite fixed point; exact CP2K bit-parity is the M1 milestone.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwConfig;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::gpw_scf;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;

namespace {

// One normalized s-type primitive at center c.
GridBasisFunction s_basis(double alpha, const double c[3]) {
    GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = c[0];
    p.center[1] = c[1];
    p.center[2] = c[2];
    p.l[0] = p.l[1] = p.l[2] = 0;
    p.coeff = std::pow(2.0 * alpha / M_PI, 0.75); // L2-normalized s-Gaussian
    return {p};
}

GpwAtom hydrogen(const double pos[3]) {
    GpwAtom a;
    a.z_val = 1.0;
    a.pos[0] = pos[0];
    a.pos[1] = pos[1];
    a.pos[2] = pos[2];
    a.basis = {s_basis(0.5, pos)};
    // GTH-PADE-q1-like local pseudopotential for H (rloc + short-range C coeffs;
    // q1 H has no nonlocal projectors). Values near the GTH H parameterization.
    a.rloc = 0.2;
    a.cloc = {-4.17890044, 0.72446331};
    return a;
}

} // namespace

TEST(GpwSolver, MinimalH2Converges) {
    const double box = 10.0;
    const double rA[3] = {box / 2 - 0.7, box / 2, box / 2};
    const double rB[3] = {box / 2 + 0.7, box / 2, box / 2}; // 1.4 Bohr apart

    std::vector<GpwAtom> atoms = {hydrogen(rA), hydrogen(rB)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = box;
    cfg.cutoff_ry = 250.0;
    cfg.n_occ = 1; // 2 electrons, closed shell
    cfg.max_iter = 100;
    cfg.tol = 1e-7;
    cfg.mix = 0.5;

    GpwResult r = gpw_scf(atoms, cfg);

    EXPECT_TRUE(r.converged) << "SCF did not converge in " << r.iters << " iters";
    EXPECT_LT(r.iters, cfg.max_iter);
    EXPECT_TRUE(std::isfinite(r.e_total)) << "E_total = " << r.e_total;
    EXPECT_TRUE(std::isfinite(r.e_kinetic) && std::isfinite(r.e_es) &&
                std::isfinite(r.e_xc));
    // Kinetic energy of a real density is positive; XC is negative.
    EXPECT_GT(r.e_kinetic, 0.0);
    EXPECT_LT(r.e_xc, 0.0);
}
