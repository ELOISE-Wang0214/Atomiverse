/// @file test_rhf.cu
/// @brief Integration tests for H2/STO-3G/RHF: integrals, CPU solver, CUDA solver, parity.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>
#include <cmath>
#include <type_traits>

#include <xtensor/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cuda/scf/scf.cuh"
#include "atomiverse/cuda/scf/rhf_fock.cuh"
#include "atomiverse/cuda/scf/rhf_density.cuh"
#include "atomiverse/cuda/recorder/thermo.cuh"

using namespace ATV;

namespace {

cpu::Frame make_h2() {
    cpu::Frame frame(2);
    frame.atoms.x[0] = double(0);     frame.atoms.y[0] = double(0); frame.atoms.z[0] = double(0);
    frame.atoms.x[1] = double(0.74);  frame.atoms.y[1] = double(0); frame.atoms.z[1] = double(0);
    frame.atoms.mass[0] = double(1.008); frame.atoms.mass[1] = double(1.008);
    frame.atoms.type[0] = 1; frame.atoms.type[1] = 1;  // H: Z=1
    frame.box = Mat3::cubic(100.0);
    return frame;
}

// Textbook H2/STO-3G/RHF reference values
constexpr double kH2Energy = -1.1167;   // approx total energy [Hartree]
constexpr double kH2Vnuc   =  0.7151;   // approx nuclear repulsion [Hartree]

} // anonymous namespace

// ── Integral tests ────────────────────────────────────────────────────

TEST(RHF, IntegralSymmetry) {
    auto frame = make_h2();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);

    int n = integrals.n_basis;
    EXPECT_EQ(n, 2);
    EXPECT_EQ(integrals.n_electrons, 2);
    EXPECT_EQ(integrals.n_occ(), 1);

    // S, Hcore should be symmetric
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            EXPECT_NEAR(integrals.S(i,j), integrals.S(j,i), 1e-12);
            EXPECT_NEAR(integrals.Hcore(i,j), integrals.Hcore(j,i), 1e-12);
        }

    // Overlap diagonal should be ~1 for normalized basis
    for (int i = 0; i < n; ++i)
        EXPECT_NEAR(static_cast<double>(integrals.S(i,i)), 1.0, 0.01);
}

TEST(RHF, NuclearRepulsion) {
    auto frame = make_h2();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    // V_nuc = Z_A * Z_B / R_AB (in Bohr) = 1*1 / (0.74 * 1.8897) ~ 0.7151
    EXPECT_NEAR(static_cast<double>(integrals.V_nuc), kH2Vnuc, 0.01);
}

// ── Orthogonalization test ────────────────────────────────────────────

TEST(RHF, Orthogonalization) {
    auto frame = make_h2();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    int n = integrals.n_basis;

    // Compute X = S^{-1/2} the same way rhf_solve does
    using RealMatrix = xt::xtensor<double, 2, xt::layout_type::column_major>;
    RealMatrix S_cm = integrals.S;
    auto [s_eig, U] = xt::linalg::eigh(S_cm);
    RealMatrix s_inv_sqrt = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        s_inv_sqrt(i, i) = double(1) / std::sqrt(s_eig(i));
    RealMatrix X = xt::linalg::dot(xt::linalg::dot(U, s_inv_sqrt), xt::transpose(U));

    // Verify X^T S X = I
    double ortho_tol = (sizeof(double) == sizeof(float)) ? 1e-6 : 1e-10;
    RealMatrix XtSX = xt::linalg::dot(xt::transpose(X), xt::linalg::dot(S_cm, X));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(static_cast<double>(XtSX(i, j)), expected, ortho_tol)
                << "X^T S X not identity at (" << i << "," << j << ")";
        }
}

// ── Fock build test ──────────────────────────────────────────────────

TEST(RHF, FockBuild) {
    auto frame = make_h2();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0, true);
    int n = integrals.n_basis;

    // Use zero density -> G(0) = 0 -> F = Hcore
    using RealMatrix = xt::xtensor<double, 2, xt::layout_type::column_major>;
    RealMatrix P_zero = xt::zeros<double>({n, n});

    // G_ij = sum_kl P_kl * [(ij|kl) - 0.5*(ik|jl)]
    // With P=0, G=0, so F = Hcore
    RealMatrix G = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double val = 0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    val += P_zero(k, l) *
                           (integrals.ERI(i, j, k, l) -
                            double(0.5) * integrals.ERI(i, k, j, l));
            G(i, j) = val;
        }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            EXPECT_NEAR(static_cast<double>(G(i, j)), 0.0, 1e-15)
                << "G(P=0) should be zero at (" << i << "," << j << ")";

    // Now use converged density and verify G is nonzero and symmetric
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    auto result = cpu::scf::rhf_solve(integrals, tol, 100);
    RealMatrix P_conv = result.P;

    RealMatrix G_conv = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double val = 0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    val += P_conv(k, l) *
                           (integrals.ERI(i, j, k, l) -
                            double(0.5) * integrals.ERI(i, k, j, l));
            G_conv(i, j) = val;
        }

    // G must be symmetric
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            EXPECT_NEAR(static_cast<double>(G_conv(i, j)),
                        static_cast<double>(G_conv(j, i)), 1e-10)
                << "G not symmetric at (" << i << "," << j << ")";

    // F = Hcore + G should give self-consistent energy
    RealMatrix F_conv = integrals.Hcore + G_conv;
    double E_elec = double(0.5) * xt::sum((integrals.Hcore + F_conv) * P_conv)();
    double E_total = E_elec + integrals.V_nuc;
    EXPECT_NEAR(static_cast<double>(E_total), kH2Energy, 0.01)
        << "Reconstructed energy from Fock build does not match";
}

// ── CPU RHF tests ─────────────────────────────────────────────────────

TEST(RHF, CPUConvergence) {
    auto frame = make_h2();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);

    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    auto result = cpu::scf::rhf_solve(integrals, tol, 100);

    EXPECT_TRUE(result.converged);
    EXPECT_TRUE(std::isfinite(static_cast<double>(result.energy)));

    // Energy should match textbook within 0.01 Hartree
    EXPECT_NEAR(static_cast<double>(result.energy), kH2Energy, 0.01)
        << "CPU RHF energy " << result.energy << " does not match textbook " << kH2Energy;
}

TEST(RHF, CPUDensityTrace) {
    auto frame = make_h2();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);

    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    auto result = cpu::scf::rhf_solve(integrals, tol, 100);

    // Tr(P * S) should equal n_electrons
    double trace_PS = xt::sum(result.P * integrals.S)();
    EXPECT_NEAR(static_cast<double>(trace_PS),
                static_cast<double>(integrals.n_electrons), 0.01);
}

// ── CUDA RHF tests ────────────────────────────────────────────────────

TEST(RHF, CUDAConvergence) {
    auto frame = make_h2();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0, true);

    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    scf::SCF solver;
    solver.add(scf::RHFFock{integrals});
    solver.add(scf::EigenSolver{});
    solver.add(scf::RHFDensity{integrals.n_occ()});
    solver.add(scf::ConvergenceCheck{tol, false});
    solver.build();
    solver.run(frame, 100);

    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
    EXPECT_NEAR(static_cast<double>(frame.energy), kH2Energy, 0.01)
        << "CUDA RHF energy " << frame.energy << " does not match textbook";
}

// ── CPU/CUDA parity ───────────────────────────────────────────────────

TEST(RHF, CPUCUDAParity) {
    auto frame = make_h2();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0, true);

    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    // CPU
    auto cpu_result = cpu::scf::rhf_solve(integrals, tol, 100);

    // CUDA
    scf::SCF solver;
    solver.add(scf::RHFFock{integrals});
    solver.add(scf::EigenSolver{});
    solver.add(scf::RHFDensity{integrals.n_occ()});
    solver.add(scf::ConvergenceCheck{tol, false});
    solver.build();
    solver.run(frame, 100);

    double cpu_e = static_cast<double>(cpu_result.energy);
    double cuda_e = static_cast<double>(frame.energy);
    double rel_tol = (sizeof(double) == sizeof(float)) ? 1e-3 : 1e-8;

    EXPECT_NEAR(cpu_e, cuda_e, std::abs(cpu_e) * rel_tol + 1e-6)
        << "CPU energy: " << cpu_e << ", CUDA energy: " << cuda_e;
}

// ── CUDA-only path (GPU integrals, zero CPU computation) ──────────────

TEST(RHF, CUDAOnlyGPUIntegrals) {
    auto frame = make_h2();

    // STO-3G parameters
    double exp[3] = {3.42525091, 0.62391373, 0.16885540};
    double coeff[3] = {0.15432897, 0.53532814, 0.44463454};

    scf::SCF solver;
    solver.add(scf::RHFFockDevice{3, exp, coeff, 1});  // n_prim=3, n_occ=1
    solver.add(scf::EigenSolver{});
    solver.add(scf::RHFDensity{1});
    solver.add(scf::ConvergenceCheck{double(1e-5), false});
    solver.build();
    solver.run(frame, 100);

    EXPECT_NEAR(static_cast<double>(frame.energy), kH2Energy, 0.01)
        << "CUDA-only RHF energy " << frame.energy << " does not match textbook";
}

// ── Multi-system tests ───────────────────────────────────────────────

static cpu::Frame make_he() {
    cpu::Frame f(1);
    f.atoms.x[0] = 0; f.atoms.y[0] = 0; f.atoms.z[0] = 0;
    f.atoms.mass[0] = double(4.003); f.atoms.type[0] = 2;
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame make_heh_plus() {
    cpu::Frame f(2);
    f.atoms.x[0] = double(0);     f.atoms.y[0] = 0; f.atoms.z[0] = 0;
    f.atoms.x[1] = double(0.772); f.atoms.y[1] = 0; f.atoms.z[1] = 0;
    f.atoms.mass[0] = double(4.003); f.atoms.mass[1] = double(1.008);
    f.atoms.type[0] = 2; f.atoms.type[1] = 1;  // HeH+ has 2 electrons
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame make_lih() {
    cpu::Frame f(2);
    f.atoms.x[0] = double(0);     f.atoms.y[0] = 0; f.atoms.z[0] = 0;
    f.atoms.x[1] = double(1.595); f.atoms.y[1] = 0; f.atoms.z[1] = 0;
    f.atoms.mass[0] = double(6.941); f.atoms.mass[1] = double(1.008);
    f.atoms.type[0] = 3; f.atoms.type[1] = 1;
    f.box = Mat3::cubic(100.0);
    return f;
}

constexpr double kHeEnergy  = -2.8078;   // He/STO-3G/RHF [Hartree]
constexpr double kHeHEnergy = -2.8609;   // HeH+/STO-3G/RHF [Hartree]
constexpr double kLiHEnergy = -7.8634;   // LiH/STO-3G/RHF [Hartree]

TEST(RHF, HeAtom) {
    auto frame = make_he();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    EXPECT_EQ(integrals.n_basis, 1);
    EXPECT_EQ(integrals.n_electrons, 2);
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    auto result = cpu::scf::rhf_solve(integrals, tol, 100);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(static_cast<double>(result.energy), kHeEnergy, 0.01)
        << "He energy: " << result.energy;
}

TEST(RHF, HeHPlus) {
    auto frame = make_heh_plus();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 1);  // charge = +1
    EXPECT_EQ(integrals.n_basis, 2);
    EXPECT_EQ(integrals.n_electrons, 2);  // He(2) + H(1) - 1 = 2
    EXPECT_EQ(integrals.n_occ(), 1);
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    auto result = cpu::scf::rhf_solve(integrals, tol, 100);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(static_cast<double>(result.energy), kHeHEnergy, 0.02)
        << "HeH+ energy: " << result.energy;
}

TEST(RHF, LiH) {
    auto frame = make_lih();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    EXPECT_EQ(integrals.n_basis, 3);  // Li: 1s+2s=2, H: 1s=1 -> total 3
    EXPECT_EQ(integrals.n_electrons, 4);
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    auto result = cpu::scf::rhf_solve(integrals, tol, 100);
    EXPECT_TRUE(result.converged);
    double lih_tol = (sizeof(double) == sizeof(float)) ? 0.1 : 0.01;
    EXPECT_NEAR(static_cast<double>(result.energy), kLiHEnergy, lih_tol)
        << "LiH energy: " << result.energy;
}

TEST(RHF, H2PESScan) {
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    double prev_energy = 0;
    double min_energy = 1e10;
    double min_r = 0;
    for (int i = 0; i < 8; ++i) {
        double r = 0.5 + i * 0.2;  // 0.5 to 1.9 A
        cpu::Frame f(2);
        f.atoms.x[0] = double(0); f.atoms.y[0] = 0; f.atoms.z[0] = 0;
        f.atoms.x[1] = double(r); f.atoms.y[1] = 0; f.atoms.z[1] = 0;
        f.atoms.mass[0] = double(1.008); f.atoms.mass[1] = double(1.008);
        f.atoms.type[0] = 1; f.atoms.type[1] = 1;
        f.box = Mat3::cubic(100.0);
        auto integrals = cpu::scf::build_integrals(f, cpu::scf::BasisSet::STO3G);
        auto result = cpu::scf::rhf_solve(integrals, tol, 100);
        EXPECT_TRUE(result.converged) << "H2 scan failed at R=" << r;
        double e = static_cast<double>(result.energy);
        EXPECT_TRUE(std::isfinite(e));
        if (e < min_energy) { min_energy = e; min_r = r; }
        prev_energy = e;
    }
    // Minimum should be near 0.7 A (textbook: 0.74 A)
    EXPECT_NEAR(min_r, 0.7, 0.15) << "PES minimum at " << min_r << " A";
}

// ── Multi-system CUDA parity ─────────────────────────────────────────

TEST(RHF, CUDAHeAtom) {
    auto frame = make_he();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0, true);
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    auto cpu_result = cpu::scf::rhf_solve(integrals, tol, 100);

    scf::SCF solver;
    solver.add(scf::RHFFock{integrals});
    solver.add(scf::EigenSolver{});
    solver.add(scf::RHFDensity{integrals.n_occ()});
    solver.add(scf::ConvergenceCheck{tol, false});
    solver.build();
    solver.run(frame, 100);

    double rel_tol = (sizeof(double) == sizeof(float)) ? 1e-3 : 1e-8;
    EXPECT_NEAR(static_cast<double>(frame.energy),
                static_cast<double>(cpu_result.energy),
                std::abs(static_cast<double>(cpu_result.energy)) * rel_tol + 1e-6)
        << "He CPU/CUDA parity";
}

TEST(RHF, CUDALiH) {
    auto frame = make_lih();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0, true);
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    auto cpu_result = cpu::scf::rhf_solve(integrals, tol, 100);

    scf::SCF solver;
    solver.add(scf::RHFFock{integrals});
    solver.add(scf::EigenSolver{});
    solver.add(scf::RHFDensity{integrals.n_occ()});
    solver.add(scf::ConvergenceCheck{tol, false});
    solver.build();
    solver.run(frame, 100);

    double rel_tol = (sizeof(double) == sizeof(float)) ? 1e-2 : 1e-6;
    EXPECT_NEAR(static_cast<double>(frame.energy),
                static_cast<double>(cpu_result.energy),
                std::abs(static_cast<double>(cpu_result.energy)) * rel_tol + 1e-4)
        << "LiH CPU=" << cpu_result.energy << " CUDA=" << frame.energy;
}

TEST(RHF, CUDAHeHPlus) {
    auto frame = make_heh_plus();
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 1, true);
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    auto cpu_result = cpu::scf::rhf_solve(integrals, tol, 100);

    scf::SCF solver;
    solver.add(scf::RHFFock{integrals});
    solver.add(scf::EigenSolver{});
    solver.add(scf::RHFDensity{integrals.n_occ()});
    solver.add(scf::ConvergenceCheck{tol, false});
    solver.build();
    solver.run(frame, 100);

    double rel_tol = (sizeof(double) == sizeof(float)) ? 1e-3 : 1e-8;
    EXPECT_NEAR(static_cast<double>(frame.energy),
                static_cast<double>(cpu_result.energy),
                std::abs(static_cast<double>(cpu_result.energy)) * rel_tol + 1e-6)
        << "HeH+ CPU/CUDA parity";
}

// ── Discipline tests ─────────────────────────────────────────────────

TEST(RHF, DriverAcceptsCPUFrameOnly) {
    // Driver::run() takes cpu::Frame& -- compile-time guarantee.
    // This test verifies the API shape at the type level.
    static_assert(
        std::is_invocable_v<decltype(&scf::SCF::run), scf::SCF*, cpu::Frame&, int>,
        "Driver::run() must accept cpu::Frame&");
}

TEST(RHF, IntegralPackHasNoCUDA) {
    // IntegralPack fields are xtensor arrays, not device pointers.
    // Verify at compile time that S is an xtensor.
    static_assert(
        std::is_same_v<
            decltype(std::declval<cpu::scf::IntegralPack>().S),
            xt::xtensor<double, 2>>,
        "IntegralPack.S must be xtensor, not a device pointer");
}

TEST(RHF, RHFResultHasNoCUDA) {
    // RHFResult fields are xtensor arrays, not device pointers.
    static_assert(
        std::is_same_v<
            decltype(std::declval<cpu::scf::RHFResult>().P),
            xt::xtensor<double, 2>>,
        "RHFResult.P must be xtensor, not a device pointer");
    static_assert(
        std::is_same_v<
            decltype(std::declval<cpu::scf::RHFResult>().eigenvalues),
            xt::xtensor<double, 1>>,
        "RHFResult.eigenvalues must be xtensor, not a device pointer");
}

#else // !ATV_HAS_RHF

#include <gtest/gtest.h>
TEST(RHF, Disabled) {
    GTEST_SKIP() << "RHF tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
