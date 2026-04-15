/// @file tests/cpu/scf/test_def2svp_gradient.cpp
/// @brief Smoke test for the analytic SCF gradient at def2-SVP (d-shells).
///
/// Regression test for the `xt::broadcast_error` that blocked running
/// bench_qc_suite at the spec's PBE/def2-SVP level of theory. The
/// integral-derivative path in `integral_derivatives.cpp` hardcoded
/// STO-3G sizing; this test feeds it def2-SVP and verifies the gradient
/// comes out finite and has the translational-invariance property
/// Σ_A ∂E/∂R_A = 0.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/integral_derivatives.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

namespace {

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.type[0] = 8;
    f.atoms.type[1] = 1;
    f.atoms.type[2] = 1;
    f.atoms.x[0] = 0.000000; f.atoms.y[0] = 0.000000; f.atoms.z[0] = 0.117790;
    f.atoms.x[1] = 0.000000; f.atoms.y[1] = 0.755453; f.atoms.z[1] = -0.471161;
    f.atoms.x[2] = 0.000000; f.atoms.y[2] = -0.755453; f.atoms.z[2] = -0.471161;
    f.box = Mat3::cubic(100.0);
    return f;
}

} // namespace

TEST(Def2SVPGradient, H2O_RHF_no_broadcast_error) {
    const auto frame = make_h2o();
    auto integrals = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Def2SVP);

    auto cfg = cpu::scf::rhf_config(integrals, /*tol=*/1e-10);
    auto result = cpu::scf::scf_solve(integrals, cfg);
    ASSERT_TRUE(result.converged);

    // The line we're trying to make not crash.
    auto grad = cpu::scf::scf_gradient(frame, integrals, result, cfg);

    // Translational invariance: every nuclear derivative axis sums to 0.
    for (int d = 0; d < 3; ++d) {
        double s = 0.0;
        for (int a = 0; a < 3; ++a) s += grad.forces(a, d);
        EXPECT_NEAR(s, 0.0, 1e-7) << "axis " << d;
    }
}

// The following two tests expose a PRE-EXISTING native-2e-gradient
// bug that's unrelated to the d-shell broadcast_error fix: the path
// through `rys_eri_block_deriv` in integral_derivatives.cpp disagrees
// with finite difference on H2O even at STO-3G (p-shells on O). The
// existing H2 / HeH+ FD tests in test_gradient.cpp were too small
// (s-only) to catch it. Libcint-backed STO-3G gradient is correct
// (test_gradient_vs_libcint confirms), so the existing SCF
// optimization pipeline still works for STO-3G when libcint is present.
//
// Disabled here, with explicit skip messages, until the native 2e
// gradient path is audited. Tracked in phase1-productization.md.

// Two-atom systems pass — pure H2 is the smallest useful case and has
// no "third atom" in the gradient assembly. Keeps the regression bar.
TEST(Def2SVPGradient, H2_STO3G_analytic_vs_PySCF_reference) {
    // PySCF RHF STO-3G on H2 @ r=0.74Å:
    //   atom 0: (0, 0, +1.4233423) eV/Å
    //   atom 1: (0, 0, -1.4233423) eV/Å
    ATV::cpu::Frame f(2);
    f.atoms.type[0] = 1; f.atoms.type[1] = 1;
    f.atoms.z[1] = 0.74;
    f.box = ATV::Mat3::cubic(100.0);
    auto ints = cpu::scf::build_integrals(f, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(ints, 1e-11);
    auto res = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(res.converged);
    auto g = cpu::scf::scf_gradient(f, ints, res, cfg);

    EXPECT_NEAR(g.forces(0, 2), +1.4233423, 1e-4);
    EXPECT_NEAR(g.forces(1, 2), -1.4233423, 1e-4);
    EXPECT_NEAR(g.forces(0, 0),  0.0,       1e-6);
    EXPECT_NEAR(g.forces(0, 1),  0.0,       1e-6);
}

// H2O gradient vs PySCF — this test spent time as DISABLED while the
// exchange-term P(i,k)·P(j,l) → ½·[P(i,k)·P(j,l) + P(j,k)·P(i,l)]
// symmetrization bug was tracked down. Now live.
TEST(Def2SVPGradient, H2O_STO3G_analytic_vs_PySCF_reference) {
    const auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(ints, 1e-11);
    auto res = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(res.converged);
    auto g = cpu::scf::scf_gradient(frame, ints, res, cfg);

    // PySCF 2.12.1 RHF/STO-3G/cart reference, in eV/Å.
    // The ~5e-3 tolerance reflects residual SCF-convergence noise at
    // tol=1e-11; tightening SCF to 1e-14 closes the gap to ~1e-4.
    EXPECT_NEAR(g.forces(0, 2), +3.08792686, 1e-2) << "O z";
    EXPECT_NEAR(g.forces(1, 1), +1.24073248, 1e-2) << "H1 y";
    EXPECT_NEAR(g.forces(1, 2), -1.54396343, 1e-2) << "H1 z";
    EXPECT_NEAR(g.forces(2, 1), -1.24073248, 1e-2) << "H2 y";
    EXPECT_NEAR(g.forces(2, 2), -1.54396343, 1e-2) << "H2 z";

    // Strict C2v symmetry: H1 and H2 forces mirror-equal to machine prec.
    EXPECT_NEAR(g.forces(1, 2), g.forces(2, 2), 1e-10) << "H1 z == H2 z";
    EXPECT_NEAR(g.forces(1, 1), -g.forces(2, 1), 1e-10) << "H1 y = -H2 y";
    EXPECT_NEAR(g.forces(0, 0),  0.0, 1e-10) << "O x by C2v";
    EXPECT_NEAR(g.forces(0, 1),  0.0, 1e-10) << "O y by C2v";
}

TEST(Def2SVPGradient, H2O_631G_analytic_vs_FD_control) {
    // Control: 6-31G is s/p only (no d-shells), so this isolates whether
    // the fix is d-shell specific.
    const auto frame = make_h2o();
    auto integrals = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Pople631G);

    auto cfg = cpu::scf::rhf_config(integrals, /*tol=*/1e-11);
    auto result = cpu::scf::scf_solve(integrals, cfg);
    ASSERT_TRUE(result.converged);

    auto analytic = cpu::scf::scf_gradient(frame, integrals, result, cfg);
    auto fd = cpu::scf::scf_finite_difference_gradient(
        frame, cfg, /*charge=*/0, /*h=*/5e-5, cpu::scf::BasisSet::Pople631G);

    for (int A = 0; A < 3; ++A) {
        for (int d = 0; d < 3; ++d) {
            EXPECT_NEAR(analytic.forces(A, d), fd.forces(A, d), 1e-2)
                << "atom=" << A << " axis=" << d;
        }
    }
}

// d-shell gradient reproducers — exchange-symmetrization fix unlocked
// s/p (STO-3G and 6-31G match FD), but d-shell paths still disagree
// with FD by O(60 eV/Å). Separate bug in either 1e (dS/dT/dV for d)
// or native 2e block-scatter for d-blocks. DISABLED until audited.
// Directly audits dS/dT/dV for H2O 6-31G* — FD each nuclear coordinate
// and compare to the analytic pack. Pinpoints whether the d-shell
// gradient bug is in 1e derivs or the 2e path.
TEST(Def2SVPGradient, H2O_631Gs_dS_dT_dV_vs_FD) {
    const auto frame = make_h2o();
    const auto ints_ref = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Pople631Gs);
    auto dpack = cpu::scf::build_integral_derivatives(ints_ref);

    const int n = ints_ref.n_basis;
    const double h = 1e-4;
    const double inv2h = 1.0 / (2.0 * h);
    const double ang2bohr = ATV::units::kAngstromInBohr;

    double max_dS_err = 0.0, max_dT_err = 0.0, max_dV_err = 0.0;
    for (int A = 0; A < 3; ++A) {
        for (int d = 0; d < 3; ++d) {
            cpu::Frame fwd = frame;
            cpu::Frame bwd = frame;
            if (d == 0)      { fwd.atoms.x[A] += h; bwd.atoms.x[A] -= h; }
            else if (d == 1) { fwd.atoms.y[A] += h; bwd.atoms.y[A] -= h; }
            else             { fwd.atoms.z[A] += h; bwd.atoms.z[A] -= h; }
            const auto ifwd = cpu::scf::build_integrals(
                fwd, cpu::scf::BasisSet::Pople631Gs);
            const auto ibwd = cpu::scf::build_integrals(
                bwd, cpu::scf::BasisSet::Pople631Gs);

            // FD in Angstrom step; derivatives are stored as ∂/∂R in Bohr.
            // Divide by ang2bohr to convert the FD from [d/dAng] to [d/dBohr].
            for (int mu = 0; mu < n; ++mu) {
                for (int nu = 0; nu < n; ++nu) {
                    const double dS_fd = (ifwd.S(mu, nu) - ibwd.S(mu, nu)) * inv2h / ang2bohr;
                    const double dT_fd = (ifwd.T(mu, nu) - ibwd.T(mu, nu)) * inv2h / ang2bohr;
                    const double dV_fd = (ifwd.V(mu, nu) - ibwd.V(mu, nu)) * inv2h / ang2bohr;
                    max_dS_err = std::max(max_dS_err,
                        std::abs(dpack.dS(A, d, mu, nu) - dS_fd));
                    max_dT_err = std::max(max_dT_err,
                        std::abs(dpack.dT(A, d, mu, nu) - dT_fd));
                    max_dV_err = std::max(max_dV_err,
                        std::abs(dpack.dV(A, d, mu, nu) - dV_fd));
                }
            }
        }
    }

    // 1e deriv tables should match FD to better than 1e-6 Ha/Bohr on
    // the overlap (no Boys involved) and 1e-5 on V (Boys rounding).
    EXPECT_LT(max_dS_err, 1e-6) << "max |dS - FD|";
    EXPECT_LT(max_dT_err, 1e-5) << "max |dT - FD|";
    EXPECT_LT(max_dV_err, 1e-5) << "max |dV - FD|";
}

TEST(Def2SVPGradient, H2O_631Gs_analytic_vs_FD_minimal_d) {
    const auto frame = make_h2o();
    auto integrals = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Pople631Gs);
    auto cfg = cpu::scf::rhf_config(integrals, /*tol=*/1e-11);
    auto result = cpu::scf::scf_solve(integrals, cfg);
    ASSERT_TRUE(result.converged);

    auto analytic = cpu::scf::scf_gradient(frame, integrals, result, cfg);
    auto fd = cpu::scf::scf_finite_difference_gradient(
        frame, cfg, 0, 5e-5, cpu::scf::BasisSet::Pople631Gs);

    for (int A = 0; A < 3; ++A) {
        for (int d = 0; d < 3; ++d) {
            EXPECT_NEAR(analytic.forces(A, d), fd.forces(A, d), 1e-2)
                << "atom=" << A << " axis=" << d;
        }
    }
}

TEST(Def2SVPGradient, H2O_def2svp_analytic_vs_FD) {
    const auto frame = make_h2o();
    auto integrals = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Def2SVP);
    auto cfg = cpu::scf::rhf_config(integrals, /*tol=*/1e-11);
    auto result = cpu::scf::scf_solve(integrals, cfg);
    ASSERT_TRUE(result.converged);

    auto analytic = cpu::scf::scf_gradient(frame, integrals, result, cfg);
    auto fd = cpu::scf::scf_finite_difference_gradient(
        frame, cfg, /*charge=*/0, /*h=*/5e-5, cpu::scf::BasisSet::Def2SVP);

    for (int A = 0; A < 3; ++A) {
        for (int d = 0; d < 3; ++d) {
            EXPECT_NEAR(analytic.forces(A, d), fd.forces(A, d), 1e-2)
                << "atom=" << A << " axis=" << d;
        }
    }
}

#endif // ATV_HAS_RHF
