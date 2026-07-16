/// @file tests/cpu/scf/test_uks_solver.cpp
/// @brief Tests for UKS configurations through the unified SCFSolver.

#ifdef ATV_HAS_RHF

#include <cmath>
#include <gtest/gtest.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/reducers/xnorm.hpp>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

class UKSSolverrTest : public ::testing::Test {
  protected:
    static constexpr double kTol =
        (sizeof(double) == sizeof(float)) ? double(1e-4) : double(1e-7);
    static constexpr double kEnergyTol = 5e-3;

    static cpu::Frame make_h_atom() {
        cpu::Frame f(1);
        f.span(ATV::kAtom_x)[0] = double(0);
        f.span(ATV::kAtom_y)[0] = double(0);
        f.span(ATV::kAtom_z)[0] = double(0);
        f.span(ATV::kAtom_type)[0] = 1;
        f.simbox = ATV::cpu::Box::cubic(100.0);
        return f;
    }

    static cpu::Frame make_he_atom() {
        cpu::Frame f(1);
        f.span(ATV::kAtom_x)[0] = double(0);
        f.span(ATV::kAtom_y)[0] = double(0);
        f.span(ATV::kAtom_z)[0] = double(0);
        f.span(ATV::kAtom_type)[0] = 2;
        f.simbox = ATV::cpu::Box::cubic(100.0);
        return f;
    }

    static cpu::Frame make_li_atom() {
        cpu::Frame f(1);
        f.span(ATV::kAtom_x)[0] = double(0);
        f.span(ATV::kAtom_y)[0] = double(0);
        f.span(ATV::kAtom_z)[0] = double(0);
        f.span(ATV::kAtom_type)[0] = 3;
        f.simbox = ATV::cpu::Box::cubic(100.0);
        return f;
    }

    static cpu::Frame make_h2(double bond_length_angstrom) {
        cpu::Frame f(2);
        f.span(ATV::kAtom_x)[0] = double(0);
        f.span(ATV::kAtom_y)[0] = double(0);
        f.span(ATV::kAtom_z)[0] = double(0);
        f.span(ATV::kAtom_x)[1] = double(bond_length_angstrom);
        f.span(ATV::kAtom_y)[1] = double(0);
        f.span(ATV::kAtom_z)[1] = double(0);
        f.span(ATV::kAtom_type)[0] = 1;
        f.span(ATV::kAtom_type)[1] = 1;
        f.simbox = ATV::cpu::Box::cubic(100.0);
        return f;
    }
};

// ── Result structure ────────────────────────────────────────────────

TEST_F(UKSSolverrTest, ResultContainsXCDecomposition) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    EXPECT_NE(static_cast<double>(result.one_electron_energy), 0.0);
    EXPECT_NE(static_cast<double>(result.coulomb_energy), 0.0);
    EXPECT_NE(static_cast<double>(result.xc_energy), 0.0);
}

TEST_F(UKSSolverrTest, ResultHasIntegratedElectronCounts) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    EXPECT_GT(static_cast<double>(result.integrated_n_alpha), 0.0);
    EXPECT_GT(static_cast<double>(result.integrated_n_beta), 0.0);
}

TEST_F(UKSSolverrTest, DensityMatricesHaveCorrectShape) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    int n = ints.n_basis;
    EXPECT_EQ(result.P_alpha.shape(0), (size_t)n);
    EXPECT_EQ(result.P_alpha.shape(1), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(0), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(1), (size_t)n);
}

TEST_F(UKSSolverrTest, EnergyDecompositionSumsCorrectly) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    double E_sum = static_cast<double>(result.one_electron_energy)
                   + static_cast<double>(result.coulomb_energy)
                   + static_cast<double>(result.xc_energy)
                   + static_cast<double>(ints.V_nuc);
    EXPECT_NEAR(E_sum, static_cast<double>(result.energy), 1e-6);
}

// ── H atom (doublet) ────────────────────────────────────────────────

TEST_F(UKSSolverrTest, HAtomLSDAConverges) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 0, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 1);
    EXPECT_EQ(result.n_beta, 0);
}

TEST_F(UKSSolverrTest, HAtomLSDAEnergyReasonable) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 0, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 50;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -0.40);
    EXPECT_GT(E, -0.55);
}

// ── Reference energies vs PySCF (libxc-xc-fix-05-uks-validation) ────
//
// All references generated by PySCF 2.13 with:
//   dft.UKS(mol), xc="lda,vwn", atom_grid=(50, 110), prune=None,
//   conv_tol=1e-12, basis="sto-3g"
// Values pinned from a one-off PySCF run; the generator script is gone by
// design — no third-party package is invoked at build or test time.
//
// **Tolerance budget**: 1e-6 .. 1e-4 Ha per system, matching the
// libxc-xc-fix-05 spec. Native uses Treutler-Ahlrichs radial quadrature
// (PySCF default per dft.gen_grid.Grids.radi_method) + tabulated
// Lebedev angular + PySCF BRAGG_RADII; XC primitives match libxc to
// 1e-12 (sub-spec 03). With the same grid scheme and the same XC
// primitives, single-atom energies match PySCF to ~1e-10 Ha
// (numerical-precision level) and molecular energies to similar
// precision modulo basis-product Coulomb integral roundoff.
namespace pyscf_ref {
constexpr double kH_atom = -4.3567023297957097e-01;  // doublet
constexpr double kHe_atom = -2.7718860444385971e+00; // singlet
constexpr double kLi_atom = -7.2213016497744018e+00; // doublet
constexpr double kH2_eq = -1.1212072895769538e+00;   // R = 0.74 A, singlet
// H2 PES scan, R (Angstrom) and singlet UKS-LSDA energy (Hartree):
struct H2Point {
    double R;
    double E;
};
constexpr H2Point kH2_pes[] = {
    {0.50, -1.0384688858111655e+00}, {0.60, -1.0997761369439225e+00},
    {0.74, -1.1212072895769538e+00}, {0.85, -1.1125034884575498e+00},
    {1.00, -1.0848978819597010e+00}, {1.25, -1.0251171926833260e+00},
    {1.50, -9.6579633051609504e-01}, {1.75, -9.1503830322964252e-01},
    {2.00, -8.7483133700381699e-01}, {2.50, -8.2254121662080815e-01},
};
} // namespace pyscf_ref

TEST_F(UKSSolverrTest, HAtomMatchesPyscfReference) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 0, cpu::scf::XCType::LSDA, 1e-10);
    config.n_radial = 50;
    config.n_angular = 110;
    config.max_iter = 200;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    ASSERT_TRUE(result.converged);
    EXPECT_LT(std::abs(result.energy - pyscf_ref::kH_atom), 1e-6)
        << "H atom energy: native=" << result.energy
        << " pyscf=" << pyscf_ref::kH_atom;
}

TEST_F(UKSSolverrTest, HeAtomMatchesPyscfReference) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, 1e-10);
    config.n_radial = 50;
    config.n_angular = 110;
    config.max_iter = 200;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    ASSERT_TRUE(result.converged);
    EXPECT_LT(std::abs(result.energy - pyscf_ref::kHe_atom), 1e-5)
        << "He atom energy: native=" << result.energy
        << " pyscf=" << pyscf_ref::kHe_atom;
}

TEST_F(UKSSolverrTest, LiAtomMatchesPyscfReference) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 2, 1, cpu::scf::XCType::LSDA, 1e-10);
    config.n_radial = 50;
    config.n_angular = 110;
    config.max_iter = 200;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    ASSERT_TRUE(result.converged);
    EXPECT_LT(std::abs(result.energy - pyscf_ref::kLi_atom), 1e-4)
        << "Li atom energy: native=" << result.energy
        << " pyscf=" << pyscf_ref::kLi_atom;
}

TEST_F(UKSSolverrTest, H2EquilibriumMatchesPyscfReference) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, 1e-10);
    config.n_radial = 50;
    config.n_angular = 110;
    config.max_iter = 200;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    ASSERT_TRUE(result.converged);
    EXPECT_LT(std::abs(result.energy - pyscf_ref::kH2_eq), 1e-5)
        << "H2 R=0.74 energy: native=" << result.energy
        << " pyscf=" << pyscf_ref::kH2_eq;
}

TEST_F(UKSSolverrTest, H2PESScanMatchesPyscfReference) {
    for (const auto& pt : pyscf_ref::kH2_pes) {
        auto frame = make_h2(pt.R);
        auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
        auto config =
            cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, 1e-10);
        config.n_radial = 50;
        config.n_angular = 110;
        config.max_iter = 200;
        cpu::scf::SCFSolver solver{config};
        solver.build(ints);
        const auto& result = solver.run(ints);
        ASSERT_TRUE(result.converged)
            << "PES point R=" << pt.R << " A failed to converge";
        EXPECT_LT(std::abs(result.energy - pt.E), 1e-4)
            << "H2 R=" << pt.R << " A: native=" << result.energy
            << " pyscf=" << pt.E << " diff=" << (result.energy - pt.E);
    }
}

TEST_F(UKSSolverrTest, HAtomSpinDensity) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 0, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    double tr_alpha = xt::sum(result.P_alpha * ints.S)();
    double tr_beta = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta), 0.0, 1e-10);
}

// ── He atom (singlet) ───────────────────────────────────────────────

TEST_F(UKSSolverrTest, HeAtomLSDAConverges) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
}

TEST_F(UKSSolverrTest, HeAtomLSDAEnergyReasonable) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 50;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -2.70);
    EXPECT_GT(E, -2.95);
}

TEST_F(UKSSolverrTest, HeClosedShellAlphaBetaSymmetry) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    double diff = xt::norm_l2(result.P_alpha - result.P_beta)();
    EXPECT_NEAR(static_cast<double>(diff), 0.0, 1e-4);
}

TEST_F(UKSSolverrTest, HeAtomElectronCount) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    double tr_alpha = xt::sum(result.P_alpha * ints.S)();
    double tr_beta = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta), 1.0, 1e-4);
}

// ── Li atom (doublet) ───────────────────────────────────────────────

TEST_F(UKSSolverrTest, LiAtomLSDAConverges) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 2, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 2);
    EXPECT_EQ(result.n_beta, 1);
}

TEST_F(UKSSolverrTest, LiAtomLSDAEnergyReasonable) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 2, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 50;
    config.n_angular = 110;
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -7.0);
    EXPECT_GT(E, -7.8);
}

// ── H2 at equilibrium ───────────────────────────────────────────────

TEST_F(UKSSolverrTest, H2EquilibriumLSDAConverges) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
}

TEST_F(UKSSolverrTest, H2EquilibriumLSDAEnergyReasonable) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 50;
    config.n_angular = 110;
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -1.05);
    EXPECT_GT(E, -1.25);
}

// ── Cross-method comparisons ────────────────────────────────────────

TEST_F(UKSSolverrTest, UKSvsUHFDifferentEnergy) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto uhf_config = cpu::scf::uhf_config(ints, 1, 1, kTol);
    uhf_config.max_iter = 50;
    cpu::scf::SCFSolver uhf_solver{uhf_config};
    uhf_solver.build(ints);
    const auto& uhf = uhf_solver.run(ints);
    auto uks_config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    uks_config.n_radial = 35;
    uks_config.n_angular = 110;
    uks_config.max_iter = 80;
    cpu::scf::SCFSolver uks_solver{uks_config};
    uks_solver.build(ints);
    const auto& uks = uks_solver.run(ints);
    EXPECT_TRUE(uhf.converged);
    EXPECT_TRUE(uks.converged);
    EXPECT_GT(std::abs(static_cast<double>(uks.energy)
                       - static_cast<double>(uhf.energy)),
              1e-4);
}

TEST_F(UKSSolverrTest, HeAtomLSDAvsPBE) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto lsda_config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    lsda_config.n_radial = 35;
    lsda_config.n_angular = 110;
    lsda_config.max_iter = 80;
    cpu::scf::SCFSolver lsda_solver{lsda_config};
    lsda_solver.build(ints);
    const auto& lsda = lsda_solver.run(ints);
    auto pbe_config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::PBE, kTol);
    pbe_config.n_radial = 35;
    pbe_config.n_angular = 110;
    pbe_config.max_iter = 80;
    cpu::scf::SCFSolver pbe_solver{pbe_config};
    pbe_solver.build(ints);
    const auto& pbe = pbe_solver.run(ints);
    EXPECT_TRUE(lsda.converged);
    EXPECT_TRUE(pbe.converged);
    EXPECT_GT(std::abs(static_cast<double>(lsda.energy)
                       - static_cast<double>(pbe.energy)),
              1e-4);
}

// ── Convergence and eigenvalues ─────────────────────────────────────

TEST_F(UKSSolverrTest, ResultHasConvergenceTrace) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_GT(result.trace.size(), 0u);
    EXPECT_EQ(static_cast<int>(result.trace.back().iter), result.iterations);
}

TEST_F(UKSSolverrTest, EigenvaluesHaveCorrectShape) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config =
        cpu::scf::uks_config(ints, 1, 1, cpu::scf::XCType::LSDA, kTol);
    config.n_radial = 35;
    config.n_angular = 110;
    config.max_iter = 80;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_EQ(result.eigenvalues_alpha.size(), (size_t)ints.n_basis);
    EXPECT_EQ(result.eigenvalues_beta.size(), (size_t)ints.n_basis);
}

#else

#include <gtest/gtest.h>
TEST(UKSSolver, SkippedNoBuildRHF) {
    GTEST_SKIP() << "ATV_BUILD_RHF=OFF — UKS tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
