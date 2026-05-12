/// @file test_uhf.cpp
/// @brief Integration tests for UHF: open-shell atoms, RHF/UHF comparison, spin diagnostics.
///
/// Required validations:
///   - H atom (doublet): open-shell n_alpha=1, n_beta=0
///   - Li atom (doublet): open-shell n_alpha=2, n_beta=1
///   - He atom: closed-shell, UHF must equal RHF
///   - H2 equilibrium: UHF must agree with RHF
///   - H2 stretched: UHF must give lower energy than RHF (symmetry breaking)
///   - Spin contamination: zero for closed-shell, nonzero for broken-symmetry solutions

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>
#include <cmath>

#include <xtensor/xnorm.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/uhf_solve.h"

using namespace ATV;

namespace {

// Convergence tolerance scaled to Real precision.
constexpr Real kTol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-9);
// Energy comparison tolerance: RHF and UHF should agree within 1e-3 Hartree on easy cases.
constexpr double kEnergyTol = 1e-3;

cpu::Frame make_h_atom() {
    cpu::Frame f(1);
    f.x[0] = Real(0); f.y[0] = Real(0); f.z[0] = Real(0);
    f.type_id[0] = 1;  // H: Z=1
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

cpu::Frame make_he_atom() {
    cpu::Frame f(1);
    f.x[0] = Real(0); f.y[0] = Real(0); f.z[0] = Real(0);
    f.type_id[0] = 2;  // He: Z=2
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

cpu::Frame make_li_atom() {
    cpu::Frame f(1);
    f.x[0] = Real(0); f.y[0] = Real(0); f.z[0] = Real(0);
    f.type_id[0] = 3;  // Li: Z=3
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

cpu::Frame make_h2(double bond_length_angstrom) {
    cpu::Frame f(2);
    f.x[0] = Real(0);                       f.y[0] = Real(0); f.z[0] = Real(0);
    f.x[1] = Real(bond_length_angstrom);    f.y[1] = Real(0); f.z[1] = Real(0);
    f.type_id[0] = 1; f.type_id[1] = 1;    // H: Z=1
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

} // namespace

// ── H atom: single alpha electron, no beta ─────────────────────────────────

TEST(UHF, HAtomConverges) {
    auto frame = make_h_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 0, kTol, 50);

    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 1);
    EXPECT_EQ(result.n_beta,  0);
}

TEST(UHF, HAtomEnergy) {
    // H atom STO-3G: energy equals the 1s core Hamiltonian element.
    // Reference: ~-0.4666 Hartree.
    auto frame = make_h_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 0, kTol, 50);

    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(static_cast<double>(result.energy), -0.4666, 5e-3);
}

TEST(UHF, HAtomElectronCount) {
    // Tr(P_alpha * S) = n_alpha = 1, Tr(P_beta * S) = n_beta = 0.
    auto frame = make_h_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 0, kTol, 50);

    Real tr_alpha = xt::sum(result.P_alpha * ints.S)();
    Real tr_beta  = xt::sum(result.P_beta  * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-5);
    EXPECT_NEAR(static_cast<double>(tr_beta),  0.0, 1e-5);
}

TEST(UHF, HAtomSpinContamination) {
    // H atom is a pure doublet: spin_contamination must be 0 (no beta electrons).
    auto frame = make_h_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 0, kTol, 50);

    EXPECT_NEAR(static_cast<double>(result.spin_contamination), 0.0, 1e-6);
}

// ── He atom: closed-shell, UHF must equal RHF ──────────────────────────────

TEST(UHF, HeAtomEqualsRHF) {
    auto frame = make_he_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto rhf   = cpu::scf::rhf_solve(ints, kTol, 50);
    auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, kTol, 50);

    EXPECT_TRUE(uhf.converged);
    EXPECT_NEAR(static_cast<double>(uhf.energy),
                static_cast<double>(rhf.energy), kEnergyTol);
}

TEST(UHF, HeAtomSpinContamination) {
    // He is closed-shell singlet: UHF spin contamination must be near 0.
    auto frame = make_he_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, kTol, 50);

    EXPECT_NEAR(static_cast<double>(uhf.spin_contamination), 0.0, 1e-4);
}

// ── Li atom: open-shell doublet ────────────────────────────────────────────

TEST(UHF, LiAtomConverges) {
    auto frame = make_li_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);

    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 2);
    EXPECT_EQ(result.n_beta,  1);
}

TEST(UHF, LiAtomElectronCount) {
    // Tr(P_alpha * S) = 2, Tr(P_beta * S) = 1.
    auto frame = make_li_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);

    Real tr_alpha = xt::sum(result.P_alpha * ints.S)();
    Real tr_beta  = xt::sum(result.P_beta  * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 2.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta),  1.0, 1e-4);
}

TEST(UHF, LiAtomEnergyRange) {
    // Li STO-3G UHF energy should be in the physically reasonable range.
    auto frame = make_li_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);

    EXPECT_TRUE(result.converged);
    EXPECT_LT(static_cast<double>(result.energy), -7.0);
    EXPECT_GT(static_cast<double>(result.energy), -8.0);
}

// ── H2 at equilibrium: UHF ≈ RHF ─────────────────────────────────────────

TEST(UHF, H2EquilibriumAgreesWithRHF) {
    // At equilibrium geometry, RHF and UHF should give essentially the same energy.
    auto frame = make_h2(0.74);
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto rhf   = cpu::scf::rhf_solve(ints, kTol, 200);
    auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, kTol, 200);

    EXPECT_TRUE(uhf.converged);
    EXPECT_NEAR(static_cast<double>(uhf.energy),
                static_cast<double>(rhf.energy), kEnergyTol);
}

TEST(UHF, H2EquilibriumSpinContamination) {
    // H2 at equilibrium is a closed-shell singlet: spin contamination near 0.
    auto frame = make_h2(0.74);
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, kTol, 200);

    EXPECT_NEAR(static_cast<double>(uhf.spin_contamination), 0.0, 1e-3);
}

// ── H2 stretched: UHF breaks spin symmetry, giving lower energy than RHF ──

TEST(UHF, H2StretchedUHFLowerThanRHF) {
    // At R=2.0 Angstrom, the UHF solution should break symmetry.
    // E_UHF <= E_RHF (variational principle guarantees this).
    auto frame = make_h2(2.0);
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto rhf   = cpu::scf::rhf_solve(ints, kTol, 50);
    auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, kTol, 100);

    EXPECT_TRUE(uhf.converged);
    // UHF energy must be at most as high as RHF (variational).
    EXPECT_LE(static_cast<double>(uhf.energy),
              static_cast<double>(rhf.energy) + 1e-6);
}

TEST(UHF, H2StretchedSpinContaminationNonzero) {
    // At large R, the UHF broken-symmetry solution has spin contamination > 0.
    auto frame = make_h2(2.0);
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, kTol, 100);

    EXPECT_TRUE(uhf.converged);
    // For strongly stretched H2, spin contamination should be positive.
    // A pure singlet would have 0, a pure triplet ~1. Expect somewhere between.
    EXPECT_GT(static_cast<double>(uhf.spin_contamination), 0.0);
}

// ── PES scan: UHF energy <= RHF energy across all bond lengths ─────────────

TEST(UHF, H2PESScanUHFVariational) {
    // Variational principle: E_UHF <= E_RHF at every geometry.
    // Near the Coulson-Fischer point (~1.2 Å for STO-3G), damped UHF needs more
    // iterations because the broken-symmetry and symmetric solutions are close
    // in energy and the HOMO-LUMO-mixed guess must slowly relax.
    for (double r : {0.6, 0.74, 1.0, 1.5, 2.0, 2.5, 3.0}) {
        auto frame = make_h2(r);
        auto ints  = cpu::scf::build_integrals_sto3g(frame);
        auto rhf   = cpu::scf::rhf_solve(ints, kTol, 200);
        auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, kTol, 400);

        EXPECT_TRUE(uhf.converged) << "UHF failed to converge at R=" << r;
        EXPECT_LE(static_cast<double>(uhf.energy),
                  static_cast<double>(rhf.energy) + 1e-5)
            << "UHF > RHF at R=" << r << " (violates variational principle)";
    }
}

// ── Result structure correctness ───────────────────────────────────────────

TEST(UHF, ResultContainsAlphaBetaEigenvalues) {
    auto frame = make_li_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);

    EXPECT_EQ(result.eigenvalues_alpha.size(), (size_t)ints.n_basis);
    EXPECT_EQ(result.eigenvalues_beta.size(),  (size_t)ints.n_basis);
}

TEST(UHF, ResultHasConvergenceTrace) {
    auto frame = make_h2(0.74);
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 1, kTol, 50);

    EXPECT_GT(result.trace.size(), 0u);
    EXPECT_EQ((int)result.trace.back().iter, result.iterations);
}

TEST(UHF, DensityMatricesHaveCorrectShape) {
    auto frame = make_li_atom();
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);

    Int n = ints.n_basis;
    EXPECT_EQ(result.P_alpha.shape(0), (size_t)n);
    EXPECT_EQ(result.P_alpha.shape(1), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(0),  (size_t)n);
    EXPECT_EQ(result.P_beta.shape(1),  (size_t)n);
}

#else

#include <gtest/gtest.h>
TEST(UHF, SkippedNoBuildRHF) {
    GTEST_SKIP() << "ATV_BUILD_RHF=OFF — UHF tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
