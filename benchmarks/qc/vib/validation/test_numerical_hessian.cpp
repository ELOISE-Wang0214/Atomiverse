/// @file tests/cpu/vib/test_numerical_hessian.cpp
/// @brief Tests for the Response driver category and the finite-difference
///        numerical Hessian (spec phase4-vibrations-01-numerical-hessian).
///
/// Coverage map to acceptance criteria:
///   ac-003 BuildsAndProducesHessian      — driver builds, [3N,3N] shape
///   ac-004 MatchesGradientDifferenceOracle — driver vs in-test scf_gradient
///                                            central-difference oracle
///   ac-005 ReturnedHessianIsSymmetric
///   ac-006 TranslationalSumRule
///   ac-007 ValidatesRequiredComponents    — missing Force / missing Hessian
///   ac-008 SingleAtomFrame

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/driver/response_driver.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/vib/hessian.h"
#include "atomiverse/cpu/vib/hessian_result.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;

namespace {

/// @brief Binding step size — must match cpu::vib::DisplaceStep (1.0e-3 Angstrom).
constexpr double kHessianStepAngstrom = 1.0e-3;

/// @brief Near-equilibrium water geometry (O at origin, two H, ~104.5 deg).
cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.x = {0.0, 0.7572, -0.7572};
    f.atoms.y = {0.0, 0.5865, 0.5865};
    f.atoms.z = {0.0, 0.0, 0.0};
    f.atoms.mass = {15.999, 1.008, 1.008};
    f.atoms.type = {8, 1, 1};
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief A single hydrogen atom (degenerate N=1 case for ac-008).
cpu::Frame make_single_atom() {
    cpu::Frame f(1);
    f.atoms.x = {0.0};
    f.atoms.y = {0.0};
    f.atoms.z = {0.0};
    f.atoms.mass = {1.008};
    f.atoms.type = {1};
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief RHF config with a tight gradient-convergence tolerance.
///
/// The numerical Hessian amplifies gradient noise `delta_g` by the
/// `1/(2h)` factor, so the spec requires every displaced geometry to be
/// converged at least as tightly as the equilibrium one (gradient
/// tolerance <= 1e-7 Ha/Bohr). `conv_tol_grad` is set explicitly so the
/// converged density is geometry-correct regardless of SCF warm-start —
/// the driver's `SCFEvaluator` warm-starts across displacements while
/// the oracle solves fresh, and only a tight tolerance makes the two
/// agree to the acceptance bar. Used identically by the driver and the
/// oracle (the "same convergence" the acceptance contract demands).
///
/// @param ints Integral pack supplying the electron count.
/// @return RHF SCFConfig with `tol = 1e-12`, `conv_tol_grad = 1e-9`.
cpu::scf::SCFConfig tight_rhf_config(const cpu::scf::IntegralPack& ints) {
    auto cfg = cpu::scf::rhf_config(ints, /*tol=*/1e-12);
    cfg.conv_tol_grad = 1e-9; // density-stationarity tolerance.
    return cfg;
}

/// @brief In-test oracle: central-difference of `scf_gradient` directly.
///
/// Reproduces, independently of the driver, the finite-difference Hessian
/// using exactly the same step `h`, the same RHF/STO-3G method, and the
/// same SCF convergence. For Cartesian coordinate `c` (atom c/3, axis c%3)
/// it displaces by +/-h, rebuilds integrals, runs SCF, computes the
/// analytic gradient, and forms the central-difference row
/// `H[c,:] = -(f_plus - f_minus) / (2h)` in eV/Angstrom^2.
///
/// @param frame   Equilibrium geometry (positions in Angstrom).
/// @param h       Displacement step in Angstrom.
/// @return Symmetrized [3N,3N] Hessian in eV/Angstrom^2.
xt::xtensor<double, 2> gradient_difference_oracle(const cpu::Frame& frame,
                                                  double h) {
    const int n = frame.atoms.size();
    const int dim = 3 * n;
    xt::xtensor<double, 2> H = xt::zeros<double>({dim, dim});

    auto force_at = [&](const cpu::Frame& geom) {
        auto ints = cpu::scf::build_integrals(geom, cpu::scf::BasisSet::STO3G);
        auto cfg = tight_rhf_config(ints);
        auto res = cpu::scf::scf_solve(ints, cfg);
        return cpu::scf::scf_gradient(geom, ints, res, cfg)
            .forces; // [n,3] eV/A
    };

    for (int c = 0; c < dim; ++c) {
        const int atom = c / 3;
        const int axis = c % 3;

        cpu::Frame fwd = frame;
        cpu::Frame bwd = frame;
        auto disp = [&](cpu::Frame& g, double s) {
            if (axis == 0)
                g.atoms.x(atom) += s;
            else if (axis == 1)
                g.atoms.y(atom) += s;
            else
                g.atoms.z(atom) += s;
        };
        disp(fwd, +h);
        disp(bwd, -h);

        auto f_plus = force_at(fwd);
        auto f_minus = force_at(bwd);

        for (int j = 0; j < n; ++j)
            for (int b = 0; b < 3; ++b)
                H(c, 3 * j + b) = -(f_plus(j, b) - f_minus(j, b)) / (2.0 * h);
    }

    // Symmetrize (same projection the driver's Symmetrize component applies).
    xt::xtensor<double, 2> Hs = xt::zeros<double>({dim, dim});
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            Hs(i, j) = 0.5 * (H(i, j) + H(j, i));
    return Hs;
}

} // namespace

// ── ac-003: driver builds and produces a [3N,3N] Hessian ────────────────

TEST(NumericalHessian, BuildsAndProducesHessian) {
    auto frame = make_h2o();

    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = tight_rhf_config(ints);

    cpu::vib::Hessian driver;
    driver.on(ATV::vib::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
    driver.build();
    auto result = driver.run(frame);

    EXPECT_EQ(result.n_atoms, 3);
    ASSERT_EQ(result.hessian.shape()[0], 9u);
    ASSERT_EQ(result.hessian.shape()[1], 9u);
    EXPECT_EQ(result.units, "eV/Angstrom^2");
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            EXPECT_TRUE(std::isfinite(result.hessian(i, j)));
}

// ── ac-004: driver Hessian matches the gradient-difference oracle ───────

TEST(NumericalHessian, MatchesGradientDifferenceOracle) {
    auto frame = make_h2o();

    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = tight_rhf_config(ints);

    cpu::vib::Hessian driver;
    driver.on(ATV::vib::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
    driver.build();
    auto result = driver.run(frame);

    auto oracle = gradient_difference_oracle(make_h2o(), kHessianStepAngstrom);

    ASSERT_EQ(result.hessian.shape()[0], oracle.shape()[0]);
    double max_err = 0.0;
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            max_err = std::max(max_err,
                               std::abs(result.hessian(i, j) - oracle(i, j)));
    EXPECT_LE(max_err, 1e-6)
        << "driver Hessian deviates from gradient-difference oracle by "
        << max_err << " eV/Angstrom^2";
}

// ── ac-005: returned Hessian is symmetric ───────────────────────────────

TEST(NumericalHessian, ReturnedHessianIsSymmetric) {
    auto frame = make_h2o();

    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = tight_rhf_config(ints);

    cpu::vib::Hessian driver;
    driver.on(ATV::vib::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
    driver.build();
    auto result = driver.run(frame);

    double max_asym = 0.0;
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            max_asym = std::max(max_asym, std::abs(result.hessian(i, j)
                                                   - result.hessian(j, i)));
    EXPECT_LE(max_asym, 1e-12);
}

// ── ac-006: translational sum rule near a stationary geometry ───────────

TEST(NumericalHessian, TranslationalSumRule) {
    auto frame = make_h2o();

    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = tight_rhf_config(ints);

    cpu::vib::Hessian driver;
    driver.on(ATV::vib::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
    driver.build();
    auto result = driver.run(frame);

    const int n = 3;
    // For each atom i and axis alpha, sum_j H[i.alpha, j.alpha] ~ 0.
    for (int i = 0; i < n; ++i)
        for (int axis = 0; axis < 3; ++axis) {
            double row_sum = 0.0;
            for (int j = 0; j < n; ++j)
                row_sum += result.hessian(3 * i + axis, 3 * j + axis);
            EXPECT_LE(std::abs(row_sum), 1e-4)
                << "translational sum rule violated for atom " << i << " axis "
                << axis << ": " << row_sum;
        }
}

// ── ac-007: ResponseDriver validates required components ─────────────────

TEST(NumericalHessian, ValidateThrowsWithoutForceEvaluator) {
    // No Capability::Force component registered → build() must throw.
    cpu::vib::Hessian driver;
    EXPECT_THROW(driver.build(), std::runtime_error);
}

TEST(NumericalHessian, ValidateThrowsWithoutHessianAccumulator) {
    // Hessian accumulator removed → no Capability::Hessian component → throw.
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = tight_rhf_config(ints);

    cpu::vib::Hessian driver{cpu::vib::Hessian::WithoutAccumulator{}};
    driver.on(ATV::vib::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
    EXPECT_THROW(driver.build(), std::runtime_error);
}

TEST(NumericalHessian, ValidateAcceptsCorrectlyConfiguredDriver) {
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = tight_rhf_config(ints);

    cpu::vib::Hessian driver;
    driver.on(ATV::vib::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
    EXPECT_NO_THROW(driver.build());
}

// ── ac-008: single-atom (N=1) frame runs without crash ──────────────────

TEST(NumericalHessian, SingleAtomFrame) {
    auto frame = make_single_atom();

    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::uhf_config(ints, 1, 0); // single H: one alpha electron

    cpu::vib::Hessian driver;
    driver.on(ATV::vib::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
    driver.build();
    auto result = driver.run(frame);

    EXPECT_EQ(result.n_atoms, 1);
    ASSERT_EQ(result.hessian.shape()[0], 3u);
    ASSERT_EQ(result.hessian.shape()[1], 3u);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_TRUE(std::isfinite(result.hessian(i, j)));
}

// ── Capability bit static presence (ac-001 support) ─────────────────────

TEST(NumericalHessian, HessianCapabilityBitIsDistinct) {
    EXPECT_EQ(static_cast<uint32_t>(Capability::Hessian), 1u << 7);
}
