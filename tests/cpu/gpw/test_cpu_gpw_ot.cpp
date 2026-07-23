/// @file tests/cpu/gpw/test_cpu_gpw_ot.cpp
/// @brief Unit tests for GPW Orbital Transformation (OT) minimizer.
///
/// Single-function units only (no full CP2K e2e — that lives in
/// regressions/gpw/ot-minimizer-water). Covers:
///   1. fixed-H OT band energy vs generalized diagonalization
///   2. mini-H2 SCF: OT vs DIIS total energy parity
///   3. GPWEvaluator::use_ot product switch (no ATV_GPW_OT env required)

#include <cmath>
#include <cstdlib>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gtest/gtest.h>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/ot.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;
using ATV::cpu::gpw::internal::ot_minimize_fixed;
using ATV::cpu::linalg::RealMatrix;
using ATV::cpu::linalg::diagonalize_fock_generalized;

namespace {

void pin_single_thread() {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
}

GridBasisFunction s_basis(double alpha, const double c[3]) {
    GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = c[0];
    p.center[1] = c[1];
    p.center[2] = c[2];
    p.l[0] = p.l[1] = p.l[2] = 0;
    p.coeff = std::pow(2.0 * alpha / M_PI, 0.75);
    return {p};
}

GpwAtom hydrogen(const double pos[3]) {
    GpwAtom a;
    a.z_val = 1.0;
    a.z = 1;
    a.pos[0] = pos[0];
    a.pos[1] = pos[1];
    a.pos[2] = pos[2];
    a.basis = {s_basis(0.5, pos)};
    a.rloc = 0.2;
    a.cloc = {-4.17890044, 0.72446331};
    return a;
}

std::vector<GpwAtom> small_h2(double box = 8.0) {
    const double rA[3] = {box / 2 - 0.7, box / 2, box / 2};
    const double rB[3] = {box / 2 + 0.7, box / 2, box / 2};
    return {hydrogen(rA), hydrogen(rB)};
}

GpwDft small_cfg(double box = 8.0) {
    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = box;
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 1;
    cfg.max_iter = 100;
    cfg.tol = 1e-7;
    cfg.mix = 0.5;
    return cfg;
}

} // namespace

// Fixed-H OT must recover the occupied band energy of the generalized
// eigenproblem (same subspace as diagonalization).
TEST(GpwOt, FixedHMatchesDiagonalizationBandEnergy) {
    pin_single_thread();
    constexpr int nb = 3;
    constexpr int no = 1;
    RealMatrix H = xt::zeros<double>({nb, nb});
    RealMatrix S = xt::zeros<double>({nb, nb});
    RealMatrix T = xt::zeros<double>({nb, nb});
    // Symmetric H with a clear gap; S = I; T = I (FULL_KINETIC precond).
    H(0, 0) = -1.5;
    H(1, 1) = -0.4;
    H(2, 2) = 0.8;
    H(0, 1) = H(1, 0) = 0.12;
    H(0, 2) = H(2, 0) = 0.05;
    H(1, 2) = H(2, 1) = -0.03;
    for (int i = 0; i < nb; ++i) {
        S(i, i) = 1.0;
        T(i, i) = 1.0;
    }

    auto [eig, C] = diagonalize_fock_generalized(H, S);
    RealMatrix C0 = xt::zeros<double>({nb, no});
    for (int u = 0; u < nb; ++u)
        C0(u, 0) = C(u, 0); // lowest eigenvector as guess (already optimal)

    // Perturb guess so OT has work to do.
    RealMatrix Cguess = C0;
    Cguess(1, 0) += 0.3;
    Cguess(2, 0) -= 0.2;

    const auto ot = ot_minimize_fixed(H, S, T, Cguess, /*tol=*/1e-8,
                                      /*max_iter=*/200);
    ASSERT_TRUE(ot.converged) << "OT fixed-H did not converge in " << ot.iters
                              << " iters";

    const double e_diag = 2.0 * eig(0); // closed-shell band energy
    EXPECT_NEAR(ot.band_energy, e_diag, 1e-6)
        << "OT band=" << ot.band_energy << " diag=" << e_diag;
}

// Mini-H2 SCF: OT and DIIS must agree on the same gapped fixed point.
TEST(GpwOt, MiniH2OtMatchesDiisEnergy) {
    pin_single_thread();
    const auto atoms = small_h2();

    GpwDft diis = small_cfg();
    diis.use_ot = false;
    diis.build(atoms);
    const GpwResult r_diis = diis.run(atoms);
    ASSERT_TRUE(r_diis.converged) << "DIIS failed in " << r_diis.iters
                                  << " iters";

    GpwDft ot = small_cfg();
    ot.use_ot = true;
    ot.max_iter = 200;
    ot.build(atoms);
    const GpwResult r_ot = ot.run(atoms);
    ASSERT_TRUE(r_ot.converged) << "OT failed in " << r_ot.iters << " iters";

    EXPECT_NEAR(r_ot.e_total, r_diis.e_total, 1e-5)
        << "OT=" << r_ot.e_total << " DIIS=" << r_diis.e_total
        << " iters ot/diis=" << r_ot.iters << "/" << r_diis.iters;
}

// Product switch is a public member (default false); ctor sets it without env.
// Full Frame+MOLOPT OT SCF is the e2e regression (≤60 s), not this unit gate.
TEST(GpwOt, EvaluatorUseOtMemberDefaultAndCtor) {
    unsetenv("ATV_GPW_OT");
    const ATV::cpu::gpw::Kind H{
        molrs::Element::H,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};

    ATV::cpu::gpw::GPWEvaluator def({H});
    EXPECT_FALSE(def.use_ot);

    ATV::cpu::gpw::GPWEvaluator on(
        /*kinds*/ {H}, /*cutoff_ry*/ 120.0, /*max_iter*/ 50, /*tol*/ 1e-6,
        /*mix*/ 0.4, ATV::cpu::gpw::XcKind::PADE, /*nonperiodic*/ false,
        /*compute_forces*/ false, /*fd_step*/ 1e-3, /*d3*/ false,
        /*use_ot*/ true);
    EXPECT_TRUE(on.use_ot);
    on.use_ot = false;
    EXPECT_FALSE(on.use_ot);
}
