/// @file tests/cpu/scf/test_koopmans.cpp
/// @brief Default-gate anchors for the Koopmans-theorem tutorial (ch02-rhf).
///
/// ms-scale RHF −ε_HOMO checks against pinned PySCF 2.13.1 (cart=True,
/// conv_tol=1e-12) values: He/H2/H2O @ STO-3G plus He @ 6-31G. The last
/// anchor also end-to-end validates the newly added He (Z=2) row in the
/// 6-31G basis-data table. The full 5×2 −ε_HOMO + E_RHF parity sweep lives
/// in benchmarks/qc/rhf/validation/test_koopmans_parity.cpp (the single
/// source of truth for the reference numbers transcribed here).

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;
using ATV::cpu::scf::BasisSet;

namespace {

cpu::Frame make_he() {
    cpu::Frame f(1);
    f.atoms.type[0] = 2; // He
    f.atoms.x[0] = 0.0;
    f.atoms.y[0] = 0.0;
    f.atoms.z[0] = 0.0;
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.atoms.type = {1, 1};
    f.atoms.x = {0.0, 0.0};
    f.atoms.y = {0.0, 0.0};
    f.atoms.z = {0.0, 0.7414};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.type = {8, 1, 1};
    f.atoms.x = {0.0, 0.0, 0.0};
    f.atoms.y = {0.0, 0.7572, -0.7572};
    f.atoms.z = {0.1173, -0.4692, -0.4692};
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief −ε_HOMO [Ha] from a TIGHTLY converged RHF solve.
///
/// Comparing an orbital energy to a tightly converged (PySCF conv_tol=1e-12)
/// reference demands that the engine's own orbitals be reproducibly converged.
/// The default gradient tolerance `sqrt(tol)` = 1e-5 (for tol=1e-10) lets the
/// dual convergence check exit on |ΔE| alone — energy converges quadratically
/// to ~1e-13 while the orbital gradient (DIIS error) is still ~1e-7 — leaving
/// orbital energies good to only ~1e-5 and, worse, making that early exit fire
/// thread-order-dependently (non-reproducible ε_HOMO). Pinning a tight
/// `conv_tol_grad` forces the unique tightly-converged fixed point, so ε_HOMO
/// is reproducible to ~1e-10 across thread counts.
double neg_eps_homo(const cpu::Frame& frame, BasisSet basis) {
    auto integrals = cpu::scf::build_integrals(frame, basis);
    auto cfg = cpu::scf::rhf_config(integrals, 1e-10);
    cfg.max_iter = 100;
    cfg.conv_tol_grad = 1e-9;
    auto result = cpu::scf::scf_solve(integrals, cfg);
    EXPECT_TRUE(result.converged);
    return -result.eigenvalues_alpha(integrals.n_occ() - 1);
}

} // namespace

// PySCF 2.13.1 cart=True conv_tol=1e-12 (see test_koopmans_parity.cpp).
TEST(Koopmans, He_STO3G_NegEpsHomo) {
    EXPECT_NEAR(neg_eps_homo(make_he(), BasisSet::STO3G), 0.87603551, 1e-8);
}

TEST(Koopmans, H2_STO3G_NegEpsHomo) {
    EXPECT_NEAR(neg_eps_homo(make_h2(), BasisSet::STO3G), 0.57797481, 1e-8);
}

TEST(Koopmans, H2O_STO3G_NegEpsHomo) {
    EXPECT_NEAR(neg_eps_homo(make_h2o(), BasisSet::STO3G), 0.39123677, 1e-8);
}

// He @ 6-31G — also end-to-end validates the new He (Z=2) 6-31G basis row.
TEST(Koopmans, He_631G_NegEpsHomo) {
    EXPECT_NEAR(neg_eps_homo(make_he(), BasisSet::Pople631G), 0.91412663,
                1e-8);
}

// The Koopmans IP in eV uses the pinned Hartree→eV factor explicitly.
TEST(Koopmans, He_STO3G_IP_eV_Conversion) {
    const double neg_eps = neg_eps_homo(make_he(), BasisSet::STO3G);
    const double ip_eV = neg_eps * 27.211386245988; // 1 Ha = 27.211386245988 eV
    EXPECT_NEAR(ip_eV, 0.87603551 * 27.211386245988, 1e-6);
}

#endif // ATV_HAS_RHF
