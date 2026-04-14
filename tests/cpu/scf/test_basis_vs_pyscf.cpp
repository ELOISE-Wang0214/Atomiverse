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

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;
using namespace ATV::cpu;
using ATV::cpu::scf::BasisSet;

namespace {

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    // Geometry copied verbatim from the PySCF reference script.
    f.atoms.type[0] = 8;  // O
    f.atoms.type[1] = 1;  // H
    f.atoms.type[2] = 1;  // H
    f.atoms.x[0] = 0.000000; f.atoms.y[0] = 0.000000; f.atoms.z[0] = 0.117790;
    f.atoms.x[1] = 0.000000; f.atoms.y[1] = 0.755453; f.atoms.z[1] = -0.471161;
    f.atoms.x[2] = 0.000000; f.atoms.y[2] = -0.755453; f.atoms.z[2] = -0.471161;
    f.box = Mat3::cubic(100.0);
    return f;
}

} // namespace

TEST(BasisVsPySCF, H2O_631G_RHF) {
    const auto frame = make_h2o();
    auto integrals = scf::build_integrals(frame, BasisSet::Pople631G);

    // Structural sanity: basis size is deterministic from geometry + basis.
    EXPECT_EQ(integrals.n_basis, 13);

    // With kAngstromToBohr pinned to PySCF's 1/0.52917721092, V_nuc
    // should match the reference to the precision of the checked-in
    // literal (~1e-12).
    EXPECT_NEAR(integrals.V_nuc, 9.189193229310, 1e-11);

    auto config = scf::rhf_config(integrals, /*tol=*/1e-10);
    auto result = scf::scf_solve(integrals, config);
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
    auto result = scf::scf_solve(integrals, config);
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
    auto result = scf::scf_solve(integrals, config);
    ASSERT_TRUE(result.converged);

    // PySCF 2.12.1 reference (cart=True):
    //   E_RHF = -75.962206378613 Ha
    EXPECT_NEAR(result.energy, -75.962206378613, 1e-7);
}

TEST(BasisVsPySCF, H2O_def2TZVP_raises_for_f_shell) {
    // def2-TZVP on O includes an f-shell (l=3). The Rys engine caps at
    // kRysMaxL=2, so build_shell_list is expected to reject the basis
    // cleanly. f-shell support is tracked in phase1-productization.md
    // Stage B2c.
    const auto frame = make_h2o();
    EXPECT_THROW(
        scf::build_integrals(frame, BasisSet::Def2TZVP),
        std::runtime_error);
}

#endif // ATV_HAS_RHF
