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

#endif // ATV_HAS_RHF
