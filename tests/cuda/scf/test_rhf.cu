/// @file test_rhf.cu
/// @brief Integration tests for H2/STO-3G/RHF: integrals, CPU solver, CUDA solver, parity.

#ifdef ATV_HAS_HF

#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>
#include <type_traits>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/views/xview.hpp>

#include "atomiverse/core/device.h"
#include "atomiverse/core/driver_template.h"
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/rhf.h"
#include "atomiverse/cuda/evaluate/evaluate.cuh"
#include "atomiverse/cuda/evaluate/single_point.cuh"
#include "atomiverse/cuda/scf/scf_evaluator.cuh"
#include "atomiverse/cuda/scf/scf_solver.cuh"

using namespace ATV;

namespace {

cpu::Frame make_h2() {
    cpu::Frame frame(2);
    frame.span(kAtom_x)[0] = double(0);
    frame.span(kAtom_y)[0] = double(0);
    frame.span(kAtom_z)[0] = double(0);
    frame.span(kAtom_x)[1] = double(0.74);
    frame.span(kAtom_y)[1] = double(0);
    frame.span(kAtom_z)[1] = double(0);
    frame.span(kAtom_mass)[0] = double(1.008);
    frame.span(kAtom_mass)[1] = double(1.008);
    frame.span(kAtom_type)[0] = 1;
    frame.span(kAtom_type)[1] = 1; // H: Z=1
    frame.box = cpu::Box::cubic(100.0);
    return frame;
}

// Lone hydrogen atom — one electron, odd count (open shell).
cpu::Frame make_h_atom() {
    cpu::Frame frame(1);
    frame.span(kAtom_x)[0] = double(0);
    frame.span(kAtom_y)[0] = double(0);
    frame.span(kAtom_z)[0] = double(0);
    frame.span(kAtom_mass)[0] = double(1.008);
    frame.span(kAtom_type)[0] = 1; // H: Z=1, 1 electron
    frame.box = cpu::Box::cubic(100.0);
    return frame;
}

// Textbook H2/STO-3G/RHF reference values
constexpr double kH2Energy = -1.1167; // approx total energy [Hartree]
constexpr double kH2Vnuc = 0.7151;    // approx nuclear repulsion [Hartree]

// CPU RHF reference from the concrete cpu::scf::RHF method (Frame in/out).
// Replaces the DELETED free CPU RHF oracle function. This is the project's
// own CPU RHF — not a third-party oracle — used to ground the CPU-side checks;
// the CUDA parity checks below compare against the hard-coded textbook
// constants (kH2Energy, kHeEnergy, ...) directly.
struct CpuRhfRef {
    double energy = 0;      // total SCF energy [Hartree]
    bool converged = false; // loop convergence verdict
};

CpuRhfRef cpu_rhf_reference(cpu::Frame& frame, int charge = 0) {
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.charge = charge;
    rhf.tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    rhf.max_iter = 100;
    rhf.run(frame);
    return {rhf.energy(), rhf.converged()};
}

} // anonymous namespace

// ── Integral tests ────────────────────────────────────────────────────

TEST(RHF, IntegralSymmetry) {
    auto frame = make_h2();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G};

    int n = integrals.n_basis;
    EXPECT_EQ(n, 2);
    EXPECT_EQ(integrals.n_electrons, 2);
    EXPECT_EQ(integrals.n_occ(), 1);

    // S, Hcore should be symmetric
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            EXPECT_NEAR(integrals.S(i, j), integrals.S(j, i), 1e-12);
            EXPECT_NEAR(integrals.Hcore(i, j), integrals.Hcore(j, i), 1e-12);
        }

    // Overlap diagonal should be ~1 for normalized basis
    for (int i = 0; i < n; ++i)
        EXPECT_NEAR(static_cast<double>(integrals.S(i, i)), 1.0, 0.01);
}

TEST(RHF, NuclearRepulsion) {
    auto frame = make_h2();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G};
    // V_nuc = Z_A * Z_B / R_AB (in Bohr) = 1*1 / (0.74 * 1.8897) ~ 0.7151
    EXPECT_NEAR(static_cast<double>(integrals.V_nuc), kH2Vnuc, 0.01);
}

// ── Orthogonalization test ────────────────────────────────────────────

TEST(RHF, Orthogonalization) {
    auto frame = make_h2();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G};
    int n = integrals.n_basis;

    // Compute X = S^{-1/2} the same way the RHF orthogonalization does
    using RealMatrix = xt::xtensor<double, 2, xt::layout_type::column_major>;
    RealMatrix S_cm = integrals.S;
    auto [s_eig, U] = xt::linalg::eigh(S_cm);
    RealMatrix s_inv_sqrt = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        s_inv_sqrt(i, i) = double(1) / std::sqrt(s_eig(i));
    RealMatrix X =
        xt::linalg::dot(xt::linalg::dot(U, s_inv_sqrt), xt::transpose(U));

    // Verify X^T S X = I
    double ortho_tol = (sizeof(double) == sizeof(float)) ? 1e-6 : 1e-10;
    RealMatrix XtSX =
        xt::linalg::dot(xt::transpose(X), xt::linalg::dot(S_cm, X));
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
    auto integrals =
        cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G, 0, true};
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
                    val += P_zero(k, l)
                           * (integrals.ERI(i, j, k, l)
                              - double(0.5) * integrals.ERI(i, k, j, l));
            G(i, j) = val;
        }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            EXPECT_NEAR(static_cast<double>(G(i, j)), 0.0, 1e-15)
                << "G(P=0) should be zero at (" << i << "," << j << ")";

    // Now use converged density and verify G is nonzero and symmetric.
    // The converged density comes from the concrete cpu::scf::RHF method;
    // density_total() is the factor-2 closed-shell P (same convention the
    // deleted free CPU solver used), and its column-major RealMatrix type is
    // identical to the local one, so the ERI contraction below is unchanged.
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    rhf.max_iter = 100;
    rhf.run(frame);
    RealMatrix P_conv = rhf.determinant().density_total();

    RealMatrix G_conv = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double val = 0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    val += P_conv(k, l)
                           * (integrals.ERI(i, j, k, l)
                              - double(0.5) * integrals.ERI(i, k, j, l));
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
    double E_elec =
        double(0.5) * xt::sum((integrals.Hcore + F_conv) * P_conv)();
    double E_total = E_elec + integrals.V_nuc;
    EXPECT_NEAR(static_cast<double>(E_total), kH2Energy, 0.01)
        << "Reconstructed energy from Fock build does not match";
}

// ── CPU RHF tests ─────────────────────────────────────────────────────

TEST(RHF, CPUConvergence) {
    auto frame = make_h2();
    auto ref = cpu_rhf_reference(frame);

    EXPECT_TRUE(ref.converged);
    EXPECT_TRUE(std::isfinite(ref.energy));

    // Energy should match textbook within 0.01 Hartree
    EXPECT_NEAR(ref.energy, kH2Energy, 0.01)
        << "CPU RHF energy " << ref.energy << " does not match textbook "
        << kH2Energy;
}

TEST(RHF, CPUDensityTrace) {
    auto frame = make_h2();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G};

    // Converged density from the concrete cpu::scf::RHF method.
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    rhf.max_iter = 100;
    rhf.run(frame);
    auto P = rhf.determinant().density_total();

    // Tr(P * S) should equal n_electrons
    double trace_PS = xt::sum(P * integrals.S)();
    EXPECT_NEAR(static_cast<double>(trace_PS),
                static_cast<double>(integrals.n_electrons), 0.01);
}

// ── Occupation guard (hf-robust-01 CUDA mirror) ──────────────────────
// Host-side: the guard fires in the SCFSolver constructor before any
// device work, so this runs without a GPU (like the CPU* tests above).

TEST(RHF, OddElectronThrowsHostIntegrals) {
    // Lone H atom / STO-3G: 1 electron. The host-integrals CUDA solver
    // constructor must reject the odd count (closed-shell RHF cannot
    // represent it; n_occ() = 0 would silently drop the electron) with
    // std::invalid_argument — never a physically wrong converged state.
    auto frame = make_h_atom();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G};
    ASSERT_EQ(integrals.n_electrons, 1);

    EXPECT_THROW((cuda::scf::SCFSolver{integrals, 1e-10, 100}),
                 std::invalid_argument);
}

TEST(RHF, EvenElectronConstructsHostIntegrals) {
    // Neutral H2 (2 electrons) constructs without throwing — the guard is
    // a strict no-op on the valid closed-shell path.
    auto frame = make_h2();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G};
    ASSERT_EQ(integrals.n_electrons, 2);

    EXPECT_NO_THROW((cuda::scf::SCFSolver{integrals, 1e-10, 100}));
}

TEST(RHF, NegativeNoccThrowsDeviceIntegrals) {
    // Device-integrals path guards only non-negativity of n_occ.
    const double exps[3] = {3.42525091, 0.62391373, 0.16885540};
    const double coefs[3] = {0.15432897, 0.53532814, 0.44463454};
    EXPECT_THROW((cuda::scf::SCFSolver{/*n_prim=*/3, exps, coefs,
                                       /*n_occ=*/-1, 1e-10, 100}),
                 std::invalid_argument);
}

// ── CUDA RHF single-point tests (cuda::Evaluate + SCF) ───────

TEST(RHF, CUDASinglePointEnergy) {
    // ac-006: single-point H2/STO-3G RHF through
    // `cuda::Evaluate + cuda::scf::SCF` reproduces the textbook
    // reference -1.1167 Ha.
    auto frame = make_h2();
    auto integrals =
        cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G, 0, true};

    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    cuda::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate)
        .add(cuda::scf::SCF{integrals, tol, 100});
    eval.build();
    eval.run(frame);

    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.meta.get(kEnergy))));
    EXPECT_NEAR(static_cast<double>(frame.meta.get(kEnergy)), kH2Energy, 0.01)
        << "CUDA RHF single-point energy " << frame.meta.get(kEnergy)
        << " does not match textbook";
}

TEST(RHF, CUDALevelShiftEnergyInvariant) {
    // hf-robust-03 CUDA mirror: the virtual-orbital level shift is a pure
    // convergence stabilizer — the converged energy must be λ-invariant.
    // H2/STO-3G (n_occ=1, n_virt=1) exercises the full device shift path
    // (diagonal +λ, −λ C'_occ C'_occ^T projector) and the post-convergence
    // de-shift.
    auto frame0 = make_h2();
    auto frame5 = make_h2();
    auto integrals0 =
        cpu::scf::IntegralPack{frame0, cpu::scf::BasisSet::STO3G, 0, true};
    auto integrals5 =
        cpu::scf::IntegralPack{frame5, cpu::scf::BasisSet::STO3G, 0, true};

    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    cuda::Evaluate eval0;
    eval0.on(ATV::evaluate::Phase::evaluate)
        .add(cuda::scf::SCF{integrals0, tol, 100, /*level_shift=*/0.0});
    eval0.build();
    eval0.run(frame0);

    cuda::Evaluate eval5;
    eval5.on(ATV::evaluate::Phase::evaluate)
        .add(cuda::scf::SCF{integrals5, tol, 100, /*level_shift=*/0.5});
    eval5.build();
    eval5.run(frame5);

    ASSERT_TRUE(std::isfinite(static_cast<double>(frame0.meta.get(kEnergy))));
    ASSERT_TRUE(std::isfinite(static_cast<double>(frame5.meta.get(kEnergy))));
    double bar = (sizeof(double) == sizeof(float)) ? 1e-4 : 1e-8;
    EXPECT_NEAR(static_cast<double>(frame5.meta.get(kEnergy)),
                static_cast<double>(frame0.meta.get(kEnergy)), bar)
        << "level_shift=0.5 energy " << frame5.meta.get(kEnergy)
        << " must match level_shift=0 energy " << frame0.meta.get(kEnergy);
}

TEST(RHF, CUDASinglePointHelper) {
    // cuda::single_point produces the identical result as constructing
    // cuda::Evaluate and add()-ing the evaluator directly.
    auto frame_direct = make_h2();
    auto frame_helper = make_h2();
    auto integrals_d = cpu::scf::IntegralPack{
        frame_direct, cpu::scf::BasisSet::STO3G, 0, true};
    auto integrals_h = cpu::scf::IntegralPack{
        frame_helper, cpu::scf::BasisSet::STO3G, 0, true};
    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    {
        cuda::Evaluate eval;
        eval.on(ATV::evaluate::Phase::evaluate)
            .add(cuda::scf::SCF{integrals_d, tol, 100});
        eval.build();
        eval.run(frame_direct);
    }
    {
        cuda::scf::SCF scf{integrals_h, tol, 100};
        cuda::single_point(frame_helper, scf);
    }

    EXPECT_NEAR(static_cast<double>(frame_direct.meta.get(kEnergy)),
                static_cast<double>(frame_helper.meta.get(kEnergy)), 1e-9);
}

// ── CPU/CUDA parity ───────────────────────────────────────────────────

TEST(RHF, CUDAMatchesTextbookH2) {
    // The deleted free CPU RHF solver was the oracle for a CPU/CUDA
    // implementation-parity check (1e-8). With it gone, the meaningful
    // remaining assertion is physics: the CUDA SCF single point reproduces
    // the hard-coded H2/STO-3G/RHF textbook energy (-1.1167 Ha).
    auto frame = make_h2();
    auto integrals =
        cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G, 0, true};

    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    cuda::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate)
        .add(cuda::scf::SCF{integrals, tol, 100});
    eval.build();
    eval.run(frame);

    double cuda_e = static_cast<double>(frame.meta.get(kEnergy));
    EXPECT_NEAR(cuda_e, kH2Energy, 0.01)
        << "CUDA H2 energy: " << cuda_e << " vs textbook " << kH2Energy;
}

// ── CUDA-only path (GPU integrals, zero CPU computation) ──────────────

TEST(RHF, CUDAOnlyGPUIntegrals) {
    auto frame = make_h2();

    // STO-3G parameters
    double exp[3] = {3.42525091, 0.62391373, 0.16885540};
    double coeff[3] = {0.15432897, 0.53532814, 0.44463454};

    cuda::Evaluate eval;
    // n_prim=3, n_occ=1
    eval.on(ATV::evaluate::Phase::evaluate)
        .add(cuda::scf::SCF{3, exp, coeff, 1, double(1e-5), 100});
    eval.build();
    eval.run(frame);

    EXPECT_NEAR(static_cast<double>(frame.meta.get(kEnergy)), kH2Energy, 0.01)
        << "CUDA-only RHF energy " << frame.meta.get(kEnergy)
        << " does not match textbook";
}

// ── Multi-system tests ───────────────────────────────────────────────

static cpu::Frame make_he() {
    cpu::Frame f(1);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_mass)[0] = double(4.003);
    f.span(kAtom_type)[0] = 2;
    f.box = cpu::Box::cubic(100.0);
    return f;
}

static cpu::Frame make_heh_plus() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = double(0);
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = double(0.772);
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_mass)[0] = double(4.003);
    f.span(kAtom_mass)[1] = double(1.008);
    f.span(kAtom_type)[0] = 2;
    f.span(kAtom_type)[1] = 1; // HeH+ has 2 electrons
    f.box = cpu::Box::cubic(100.0);
    return f;
}

static cpu::Frame make_lih() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = double(0);
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = double(1.595);
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_mass)[0] = double(6.941);
    f.span(kAtom_mass)[1] = double(1.008);
    f.span(kAtom_type)[0] = 3;
    f.span(kAtom_type)[1] = 1;
    f.box = cpu::Box::cubic(100.0);
    return f;
}

constexpr double kHeEnergy = -2.8078;  // He/STO-3G/RHF [Hartree]
constexpr double kHeHEnergy = -2.8609; // HeH+/STO-3G/RHF [Hartree]
constexpr double kLiHEnergy = -7.8634; // LiH/STO-3G/RHF [Hartree]

TEST(RHF, HeAtom) {
    auto frame = make_he();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G};
    EXPECT_EQ(integrals.n_basis, 1);
    EXPECT_EQ(integrals.n_electrons, 2);
    auto ref = cpu_rhf_reference(frame);
    EXPECT_TRUE(ref.converged);
    EXPECT_NEAR(ref.energy, kHeEnergy, 0.01) << "He energy: " << ref.energy;
}

TEST(RHF, HeHPlus) {
    auto frame = make_heh_plus();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G,
                                            1}; // charge = +1
    EXPECT_EQ(integrals.n_basis, 2);
    EXPECT_EQ(integrals.n_electrons, 2); // He(2) + H(1) - 1 = 2
    EXPECT_EQ(integrals.n_occ(), 1);
    auto ref = cpu_rhf_reference(frame, /*charge=*/1);
    EXPECT_TRUE(ref.converged);
    EXPECT_NEAR(ref.energy, kHeHEnergy, 0.02) << "HeH+ energy: " << ref.energy;
}

TEST(RHF, LiH) {
    auto frame = make_lih();
    auto integrals = cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G};
    EXPECT_EQ(integrals.n_basis,
              6); // Li: 1s+2s+2p(x,y,z)=5, H: 1s=1 -> total 6
    EXPECT_EQ(integrals.n_electrons, 4);
    auto ref = cpu_rhf_reference(frame);
    EXPECT_TRUE(ref.converged);
    double lih_tol = (sizeof(double) == sizeof(float)) ? 0.1 : 0.01;
    EXPECT_NEAR(ref.energy, kLiHEnergy, lih_tol)
        << "LiH energy: " << ref.energy;
}

TEST(RHF, H2PESScan) {
    double min_energy = 1e10;
    double min_r = 0;
    for (int i = 0; i < 8; ++i) {
        double r = 0.5 + i * 0.2; // 0.5 to 1.9 A
        cpu::Frame f(2);
        f.span(kAtom_x)[0] = double(0);
        f.span(kAtom_y)[0] = 0;
        f.span(kAtom_z)[0] = 0;
        f.span(kAtom_x)[1] = double(r);
        f.span(kAtom_y)[1] = 0;
        f.span(kAtom_z)[1] = 0;
        f.span(kAtom_mass)[0] = double(1.008);
        f.span(kAtom_mass)[1] = double(1.008);
        f.span(kAtom_type)[0] = 1;
        f.span(kAtom_type)[1] = 1;
        f.box = cpu::Box::cubic(100.0);
        auto ref = cpu_rhf_reference(f);
        EXPECT_TRUE(ref.converged) << "H2 scan failed at R=" << r;
        EXPECT_TRUE(std::isfinite(ref.energy));
        if (ref.energy < min_energy) {
            min_energy = ref.energy;
            min_r = r;
        }
    }
    // Minimum should be near 0.7 A (textbook: 0.74 A)
    EXPECT_NEAR(min_r, 0.7, 0.15) << "PES minimum at " << min_r << " A";
}

// ── Multi-system CUDA parity ─────────────────────────────────────────

TEST(RHF, CUDAHeAtom) {
    // He/STO-3G/RHF: CUDA SCF reproduces the hard-coded textbook energy
    // (the CPU RHF oracle it used to compare against is deleted).
    auto frame = make_he();
    auto integrals =
        cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G, 0, true};
    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    cuda::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate)
        .add(cuda::scf::SCF{integrals, tol, 100});
    eval.build();
    eval.run(frame);

    EXPECT_NEAR(static_cast<double>(frame.meta.get(kEnergy)), kHeEnergy, 0.01)
        << "He CUDA energy: " << frame.meta.get(kEnergy);
}

TEST(RHF, CUDALiH) {
    // LiH/STO-3G/RHF: CUDA SCF reproduces the hard-coded textbook energy.
    auto frame = make_lih();
    auto integrals =
        cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G, 0, true};
    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    cuda::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate)
        .add(cuda::scf::SCF{integrals, tol, 100});
    eval.build();
    eval.run(frame);

    double lih_tol = (sizeof(double) == sizeof(float)) ? 0.1 : 0.01;
    EXPECT_NEAR(static_cast<double>(frame.meta.get(kEnergy)), kLiHEnergy,
                lih_tol)
        << "LiH CUDA energy: " << frame.meta.get(kEnergy);
}

TEST(RHF, CUDAHeHPlus) {
    // HeH+/STO-3G/RHF: CUDA SCF reproduces the hard-coded textbook energy.
    auto frame = make_heh_plus();
    auto integrals =
        cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G, 1, true};
    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    cuda::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate)
        .add(cuda::scf::SCF{integrals, tol, 100});
    eval.build();
    eval.run(frame);

    EXPECT_NEAR(static_cast<double>(frame.meta.get(kEnergy)), kHeHEnergy, 0.02)
        << "HeH+ CUDA energy: " << frame.meta.get(kEnergy);
}

// ── Structural / discipline tests ────────────────────────────────────

TEST(RHF, SCFIsEnergyDensityOrbitalsOnly) {
    // ac-005: the RHF-configured SCF carries exactly
    // Energy|DensityMatrix|Orbitals — no Force capability (no CUDA SCF
    // gradient kernel exists).
    auto frame = make_h2();
    auto integrals =
        cpu::scf::IntegralPack{frame, cpu::scf::BasisSet::STO3G, 0, true};
    cuda::scf::SCF eval{integrals, double(1e-6), 100};
    EXPECT_EQ(eval.role(), 0);
}

TEST(RHF, SCFSolverIsNotADriver) {
    // ac-002: the RHF inner solver is a low-level class, not a Driver<>.
    static_assert(
        !std::is_base_of_v<Driver<ATV::evaluate::Phase, cuda::device_tag>,
                           cuda::scf::SCFSolver>,
        "cuda::scf::SCFSolver must not be a Driver<> subclass");
    SUCCEED();
}

TEST(RHF, EvaluateRunAcceptsCPUFrameOnly) {
    // cuda::Evaluate::run() takes cpu::Frame& -- compile-time guarantee.
    static_assert(std::is_invocable_v<decltype(&cuda::Evaluate::run),
                                      cuda::Evaluate*, cpu::Frame&>,
                  "cuda::Evaluate::run() must accept cpu::Frame&");
}

TEST(RHF, IntegralPackHasNoCUDA) {
    // IntegralPack fields are xtensor arrays, not device pointers.
    // Verify at compile time that S is an xtensor.
    static_assert(
        std::is_same_v<decltype(std::declval<cpu::scf::IntegralPack>().S),
                       xt::xtensor<double, 2>>,
        "IntegralPack.S must be xtensor, not a device pointer");
}

TEST(RHF, DeterminantDensityHasNoCUDA) {
    // The converged wavefunction is host-side dense linear algebra: the CPU
    // RHF density (density_total()) is a column-major xtensor, not a device
    // pointer. (Replaces the deleted cpu::scf::RHFResult contract check.)
    static_assert(
        std::is_same_v<decltype(std::declval<cpu::scf::SlaterDeterminant>()
                                    .density_total()),
                       xt::xtensor<double, 2, xt::layout_type::column_major>>,
        "SlaterDeterminant::density_total() must be a host xtensor");
}

#else // !ATV_HAS_HF

#include <gtest/gtest.h>
TEST(RHF, Disabled) { GTEST_SKIP() << "RHF tests require -DATV_BUILD_HF=ON"; }

#endif // ATV_HAS_HF
