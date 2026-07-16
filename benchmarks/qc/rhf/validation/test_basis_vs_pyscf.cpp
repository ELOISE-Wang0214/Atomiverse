/// @file tests/cpu/scf/test_basis_vs_pyscf.cpp
/// @brief PySCF parity for multi-basis SCF (Phase 1 Stage B2 acceptance).
///
/// Ground truth is PySCF 2.12.x with the atom geometries checked into
/// this file (both sides pin the same Angstrom→Bohr constant and the
/// same Pople "unnormalized" contraction coefficients, so agreement
/// should be at the numerical-precision level).
///
/// Stage B2a covers 6-31G (l ∈ {0, 1}). 6-31G*/def2-SVP/def2-TZVP
/// land in Stage B2b once the d-shell normalization plumbing exists.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;
using namespace ATV::cpu;
using ATV::cpu::scf::BasisSet;

namespace {

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    // Geometry copied verbatim from the PySCF reference script.
    f.span(ATV::kAtom_type)[0] = 8; // O
    f.span(ATV::kAtom_type)[1] = 1; // H
    f.span(ATV::kAtom_type)[2] = 1; // H
    f.span(ATV::kAtom_x)[0] = 0.000000;
    f.span(ATV::kAtom_y)[0] = 0.000000;
    f.span(ATV::kAtom_z)[0] = 0.117790;
    f.span(ATV::kAtom_x)[1] = 0.000000;
    f.span(ATV::kAtom_y)[1] = 0.755453;
    f.span(ATV::kAtom_z)[1] = -0.471161;
    f.span(ATV::kAtom_x)[2] = 0.000000;
    f.span(ATV::kAtom_y)[2] = -0.755453;
    f.span(ATV::kAtom_z)[2] = -0.471161;
    f.simbox = ATV::cpu::Box::cubic(100.0);
    return f;
}

} // namespace

TEST(BasisVsPySCF, H2O_631G_RHF) {
    const auto frame = make_h2o();
    auto integrals = scf::build_integrals(frame, BasisSet::Pople631G);

    // Structural sanity: basis size is deterministic from geometry + basis.
    EXPECT_EQ(integrals.n_basis, 13);

    // With units::kAngstromInBohr pinned to PySCF's 1/0.52917721092, V_nuc
    // should match the reference to the precision of the checked-in
    // literal (~1e-12).
    EXPECT_NEAR(integrals.V_nuc, 9.189193229310, 1e-11);

    auto config = scf::rhf_config(integrals, /*tol=*/1e-10);
    ATV::cpu::scf::SCFSolver result_solver{config};
    result_solver.build(integrals);
    const auto& result = result_solver.run(integrals);
    ASSERT_TRUE(result.converged);

    // PySCF 2.12.1 reference:
    //   E_RHF = -75.983831120632 Ha
    EXPECT_NEAR(result.energy, -75.983831120632, 1e-7);
}

TEST(BasisVsPySCF, H2O_631Gs_cart_overlap) {
    // Structural sanity: overlap diagonal on Cartesian d follows the
    // libcint convention (<d_xy|d_xy>=1, <d_xx|d_xx>=3).
    const auto frame = make_h2o();
    auto integrals = scf::build_integrals(frame, BasisSet::Pople631Gs);

    // Find O's d AOs in atomiverse shell order. H shells come after O
    // (O has 9 s/p + 6 d = 15 AOs; H2 adds 4 more = 19).
    // AO layout at O (Z=8): [1s, 2s, 3s, 2p{x,y,z}, 3p{x,y,z}, 3d{xx,xy,xz,yy,yz,zz}]
    // d block starts at AO index 9.
    const int d0 = 9;
    EXPECT_EQ(integrals.shells[d0 + 0].l, 2);    // d_xx
    EXPECT_EQ(integrals.shells[d0 + 1].cart, 1); // d_xy
    EXPECT_NEAR(integrals.S(d0 + 0, d0 + 0), 3.0, 1e-10) << "<d_xx|d_xx>";
    EXPECT_NEAR(integrals.S(d0 + 1, d0 + 1), 1.0, 1e-10) << "<d_xy|d_xy>";
    EXPECT_NEAR(integrals.S(d0 + 2, d0 + 2), 1.0, 1e-10) << "<d_xz|d_xz>";
    EXPECT_NEAR(integrals.S(d0 + 3, d0 + 3), 3.0, 1e-10) << "<d_yy|d_yy>";
    EXPECT_NEAR(integrals.S(d0 + 4, d0 + 4), 1.0, 1e-10) << "<d_yz|d_yz>";
    EXPECT_NEAR(integrals.S(d0 + 5, d0 + 5), 3.0, 1e-10) << "<d_zz|d_zz>";
    // Off-diagonal <d_xx|d_yy> = 1 in libcint cart.
    EXPECT_NEAR(integrals.S(d0 + 0, d0 + 3), 1.0, 1e-10) << "<d_xx|d_yy>";
}

TEST(BasisVsPySCF, H2O_631Gs_cart_RHF) {
    const auto frame = make_h2o();
    auto integrals = scf::build_integrals(frame, BasisSet::Pople631Gs);

    // Cartesian d on O (no d on H): 13 (6-31G) + 6 (d on O) = 19.
    EXPECT_EQ(integrals.n_basis, 19);
    EXPECT_NEAR(integrals.V_nuc, 9.189193229310, 1e-11);

    auto config = scf::rhf_config(integrals, /*tol=*/1e-10);
    ATV::cpu::scf::SCFSolver result_solver{config};
    result_solver.build(integrals);
    const auto& result = result_solver.run(integrals);
    ASSERT_TRUE(result.converged);

    // PySCF 2.12.1 reference (cart=True to match our Cartesian d):
    //   E_RHF = -76.010481563475 Ha
    EXPECT_NEAR(result.energy, -76.010481563475, 1e-7);
}

TEST(BasisVsPySCF, H2O_def2SVP_cart_RHF) {
    const auto frame = make_h2o();
    auto integrals = scf::build_integrals(frame, BasisSet::Def2SVP);

    // O: 3s + 2p + 1d (cart) = 3 + 6 + 6 = 15 AOs; H: 2s + 1p = 2 + 3 = 5.
    // Total: 15 + 2*5 = 25.
    EXPECT_EQ(integrals.n_basis, 25);

    auto config = scf::rhf_config(integrals, /*tol=*/1e-10);
    ATV::cpu::scf::SCFSolver result_solver{config};
    result_solver.build(integrals);
    const auto& result = result_solver.run(integrals);
    ASSERT_TRUE(result.converged);

    // PySCF 2.12.1 reference (cart=True):
    //   E_RHF = -75.962206378613 Ha
    EXPECT_NEAR(result.energy, -75.962206378613, 1e-7);
}

TEST(BasisVsPySCF, H2O_def2TZVP_f_shell_build_integrals) {
    // def2-TZVP on O includes an f-shell (l=3). Since the f-shell Rys
    // support landed (kRysMaxL=3, kRysMaxRoots=7, iterative 6/7-root
    // finder), build_integrals must run to completion on a real f-shell
    // basis and yield a fully finite, non-zero integral pack. Energy
    // parity vs PySCF for def2-TZVP is deferred to a later sub-spec; this
    // test guards the integral-engine end-to-end path only.
    const auto frame = make_h2o();
    auto integrals = scf::build_integrals(frame, BasisSet::Def2TZVP,
                                          /*charge=*/0,
                                          /*compute_eri=*/true);
    ASSERT_GT(integrals.n_basis, 0);

    const int n = integrals.n_basis;
    bool eri_nonzero = false;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l) {
                    const double v = integrals.ERI(i, j, k, l);
                    ASSERT_TRUE(std::isfinite(v))
                        << "non-finite def2-TZVP ERI at (" << i << "," << j
                        << "," << k << "," << l << ")";
                    if (std::abs(v) > 1e-12)
                        eri_nonzero = true;
                }
    EXPECT_TRUE(eri_nonzero) << "def2-TZVP ERI tensor is all-zero";

    // Self-repulsion integrals (ii|ii) are strictly positive.
    for (int i = 0; i < n; ++i)
        EXPECT_GT(integrals.ERI(i, i, i, i), 0.0)
            << "(ii|ii) must be positive at i=" << i;
}

#endif // ATV_HAS_RHF
