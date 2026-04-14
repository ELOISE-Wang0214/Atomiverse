/// @file tests/cpu/scf/test_def2svp_gradient.cpp
/// @brief Smoke test for the analytic SCF gradient at def2-SVP (d-shells).
///
/// Regression test for the `xt::broadcast_error` that blocked running
/// bench_qc_suite at the spec's PBE/def2-SVP level of theory. The
/// integral-derivative path in `integral_derivatives.cpp` hardcoded
/// `sto3g_data`; this test feeds it def2-SVP and verifies the gradient
/// comes out finite and has the translational-invariance property
/// Σ_A ∂E/∂R_A = 0.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integrals.h"
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

// ── Pre-existing bug reproducer ──────────────────────────────────────
//
// Three-atom systems produce wildly wrong gradients. This is NOT
// d-shell-specific and NOT introduced by the recent IntegralPack-based
// `build_integral_derivatives` refactor — reverting that refactor and
// forcing libcint / native 2e paths all give the same wrong answer.
// The existing H2 / HeH+ tests in test_gradient.cpp are too small to
// expose it; `test_gradient_vs_libcint` only verifies the two paths
// agree with EACH OTHER (they do), not vs physics.
//
// Comparison on H2O STO-3G (RHF, cart):
//   PySCF / FD:   atom 0 (O) z = +3.088,  atom 1 (H) z = -1.544
//   Atomiverse:   atom 0 (O) z = -8.367,  atom 1 (H) z = +4.661
//
// Ratio / sign does not match any obvious factor-of-2 symmetry slip.
// Likely places to audit (both native and libcint branches share them):
//   * `accumulate_common_gradient` — Tr[P·dH] − Tr[W·dS] + dV_nuc
//     (the Pulay / 1e assembly in gradient.cpp)
//   * `accumulate_two_electron_gradient_sto3g` atom dispatch — atom_i
//     from shell_to_atom vs block_to_atom indexing
//   * Nuclear-repulsion derivative — contributes only at the atoms it
//     actually involves; cancellation across atoms could go wrong
//
// Tracked in phase1-productization.md §1.2.
TEST(Def2SVPGradient, DISABLED_H2O_STO3G_analytic_vs_PySCF_reference) {
    const auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(ints, 1e-11);
    auto res = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(res.converged);
    auto g = cpu::scf::scf_gradient(frame, ints, res, cfg);

    EXPECT_NEAR(g.forces(0, 2), +3.08792686, 1e-4) << "O z";
    EXPECT_NEAR(g.forces(1, 1), +1.24073248, 1e-4) << "H1 y";
    EXPECT_NEAR(g.forces(1, 2), -1.54396343, 1e-4) << "H1 z";
    EXPECT_NEAR(g.forces(2, 1), -1.24073248, 1e-4) << "H2 y";
    EXPECT_NEAR(g.forces(2, 2), -1.54396343, 1e-4) << "H2 z";
}

// Skipped until the 3-atom bug above is fixed.
TEST(Def2SVPGradient, DISABLED_H2O_RHF_def2svp_analytic_vs_FD) {
    GTEST_SKIP() << "Blocked on H2O_STO3G_analytic_vs_PySCF_reference "
                    "(same gradient-assembly bug, not d-shell specific).";
}

#endif // ATV_HAS_RHF
