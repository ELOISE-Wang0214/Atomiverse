/// @file test_uks.cpp
/// @brief Integration tests for UKS: LSDA/PBE on H, He, Li atoms and H2 molecule.
///
/// Validates:
///   - SCF convergence for closed-shell and open-shell systems
///   - Energy decomposition: E_1e + E_J + E_xc + V_nuc = E_total
///   - Integrated electron counts from grid
///   - UKS vs UHF energy difference (different methods)
///   - LSDA vs PBE energy difference (functional dependence)
///   - Alpha/beta symmetry for closed-shell
///   - Spin density for open-shell

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>
#include <cmath>

#include <xtensor/xnorm.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/uhf_solve.h"
#include "atomiverse/cpu/scf/uks_solve.h"

using namespace ATV;

namespace {

constexpr Real kTol = (sizeof(Real) == sizeof(float)) ? Real(1e-4) : Real(1e-7);
constexpr double kEnergyTol = 5e-3;  // 5 mHa tolerance for grid-based method

cpu::Frame make_h_atom() {
    cpu::Frame f(1);
    f.x[0] = Real(0); f.y[0] = Real(0); f.z[0] = Real(0);
    f.type_id[0] = 1;
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

cpu::Frame make_he_atom() {
    cpu::Frame f(1);
    f.x[0] = Real(0); f.y[0] = Real(0); f.z[0] = Real(0);
    f.type_id[0] = 2;
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

cpu::Frame make_li_atom() {
    cpu::Frame f(1);
    f.x[0] = Real(0); f.y[0] = Real(0); f.z[0] = Real(0);
    f.type_id[0] = 3;
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

cpu::Frame make_h2(double bond_length_angstrom) {
    cpu::Frame f(2);
    f.x[0] = Real(0);                    f.y[0] = Real(0); f.z[0] = Real(0);
    f.x[1] = Real(bond_length_angstrom); f.y[1] = Real(0); f.z[1] = Real(0);
    f.type_id[0] = 1; f.type_id[1] = 1;
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

} // namespace

// ── Result structure correctness ────────────────────────────────────────

TEST(UKS, ResultContainsXCDecomposition) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    // All energy components should be populated
    EXPECT_NE(static_cast<double>(result.one_electron_energy), 0.0);
    EXPECT_NE(static_cast<double>(result.coulomb_energy), 0.0);
    EXPECT_NE(static_cast<double>(result.xc_energy), 0.0);
}

TEST(UKS, ResultHasIntegratedElectronCounts) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    EXPECT_GT(static_cast<double>(result.integrated_n_alpha), 0.0);
    EXPECT_GT(static_cast<double>(result.integrated_n_beta), 0.0);
}

TEST(UKS, DensityMatricesHaveCorrectShape) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    Int n = ints.n_basis;
    EXPECT_EQ(result.P_alpha.shape(0), (size_t)n);
    EXPECT_EQ(result.P_alpha.shape(1), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(0), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(1), (size_t)n);
}

TEST(UKS, EnergyDecompositionSumsCorrectly) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    double E_sum = static_cast<double>(result.one_electron_energy)
                 + static_cast<double>(result.coulomb_energy)
                 + static_cast<double>(result.xc_energy)
                 + static_cast<double>(ints.V_nuc);
    EXPECT_NEAR(E_sum, static_cast<double>(result.energy), 1e-6)
        << "Energy decomposition inconsistent";
}

// ── H atom (doublet: n_α=1, n_β=0) ─────────────────────────────────────

TEST(UKS, HAtomLSDAConverges) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 0, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 1);
    EXPECT_EQ(result.n_beta, 0);
}

TEST(UKS, HAtomLSDAEnergyReasonable) {
    // H atom UKS-LSDA/STO-3G should give energy around -0.47 Hartree.
    // PySCF reference: UKS(lda,vwn5)/STO-3G → ~-0.4696 Ha
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 0, cpu::scf::XCType::LSDA, 50, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -0.40) << "H atom LSDA energy too high";
    EXPECT_GT(E, -0.55) << "H atom LSDA energy too low";
}

TEST(UKS, HAtomSpinDensity) {
    // H atom: P_beta should be zero, P_alpha should be nonzero.
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 0, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    Real tr_alpha = xt::sum(result.P_alpha * ints.S)();
    Real tr_beta  = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta), 0.0, 1e-10);
}

// ── He atom (singlet: n_α=1, n_β=1) ────────────────────────────────────

TEST(UKS, HeAtomLSDAConverges) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
}

TEST(UKS, HeAtomLSDAEnergyReasonable) {
    // He atom UKS-LSDA/STO-3G: expect around -2.81 to -2.84 Hartree.
    // PySCF reference: UKS(lda,vwn5)/STO-3G → ~-2.8267 Ha
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 50, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -2.70) << "He atom LSDA energy too high";
    EXPECT_GT(E, -2.95) << "He atom LSDA energy too low";
}

TEST(UKS, HeClosedShellAlphaBetaSymmetry) {
    // He is closed-shell: P_alpha should equal P_beta.
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    Real diff = xt::norm_l2(result.P_alpha - result.P_beta)();
    EXPECT_NEAR(static_cast<double>(diff), 0.0, 1e-4)
        << "Alpha/beta density matrices differ for closed-shell He";
}

TEST(UKS, HeAtomElectronCount) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(result.converged);
    Real tr_alpha = xt::sum(result.P_alpha * ints.S)();
    Real tr_beta  = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta), 1.0, 1e-4);
}

// ── Li atom (doublet: n_α=2, n_β=1) ────────────────────────────────────

TEST(UKS, LiAtomLSDAConverges) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 2, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 100);

    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 2);
    EXPECT_EQ(result.n_beta, 1);
}

TEST(UKS, LiAtomLSDAEnergyReasonable) {
    // Li atom UKS-LSDA/STO-3G: expect around -7.2 to -7.4 Hartree.
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 2, 1, cpu::scf::XCType::LSDA, 50, 110, kTol, 100);

    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -7.0) << "Li atom LSDA energy too high";
    EXPECT_GT(E, -7.8) << "Li atom LSDA energy too low";
}

// ── H2 at equilibrium ───────────────────────────────────────────────────

TEST(UKS, H2EquilibriumLSDAConverges) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 100);

    EXPECT_TRUE(result.converged);
}

TEST(UKS, H2EquilibriumLSDAEnergyReasonable) {
    // H2 UKS-LSDA/STO-3G at 0.74 Å: expect around -1.13 to -1.15 Hartree.
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 50, 110, kTol, 100);

    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -1.05) << "H2 LSDA energy too high";
    EXPECT_GT(E, -1.25) << "H2 LSDA energy too low";
}

// ── UKS vs UHF comparison ──────────────────────────────────────────────

TEST(UKS, UKSvsUHFDifferentEnergy) {
    // UKS and UHF are different methods — they should give different energies.
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 50);
    auto uks = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_TRUE(uhf.converged);
    EXPECT_TRUE(uks.converged);
    // Energies should differ (UKS has XC correlation, UHF has exact exchange)
    EXPECT_GT(std::abs(static_cast<double>(uks.energy) - static_cast<double>(uhf.energy)), 1e-4);
}

// ── LSDA vs PBE comparison ─────────────────────────────────────────────

TEST(UKS, HeAtomLSDAvsPBE) {
    // PBE should give a different energy than LSDA.
    // Use single atom (He) where PBE converges reliably with damped SCF.
    // Multi-atom PBE requires DIIS for stable convergence (future work).
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto lsda = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    auto pbe  = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::PBE, 35, 110, kTol, 80);

    EXPECT_TRUE(lsda.converged);
    EXPECT_TRUE(pbe.converged);
    EXPECT_GT(std::abs(static_cast<double>(lsda.energy) - static_cast<double>(pbe.energy)), 1e-4);
}

// ── Convergence trace ──────────────────────────────────────────────────

TEST(UKS, ResultHasConvergenceTrace) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_GT(result.trace.size(), 0u);
    EXPECT_EQ(static_cast<int>(result.trace.back().iter), result.iterations);
}

TEST(UKS, EigenvaluesHaveCorrectShape) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);

    EXPECT_EQ(result.eigenvalues_alpha.size(), (size_t)ints.n_basis);
    EXPECT_EQ(result.eigenvalues_beta.size(), (size_t)ints.n_basis);
}

#else

#include <gtest/gtest.h>
TEST(UKS, SkippedNoBuildRHF) {
    GTEST_SKIP() << "ATV_BUILD_RHF=OFF — UKS tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
