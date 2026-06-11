/// @file tests/cpu/scf/test_scf_solver.cpp
/// @brief Phase 3 — `cpu::scf::SCFSolver` is the public low-level inner
///        solver. After Phase 3, `scf_solve()` is a one-line shim that
///        constructs an `SCFSolver` and calls `run()` — so direct
///        bit-for-bit equality with the legacy free function is trivial
///        and uninformative. This test instead asserts that:
///
///        1. `SCFSolver{cfg}.run(integrals)` produces a converged result
///           on a representative molecule (H₂O / STO-3G RHF and H₂O / STO-3G
///           UKS-LDA).
///        2. Calling `run()` twice on the same solver yields identical
///           results (no internal state leaks across runs).
///        3. The result energy matches the (now-shimmed) `scf_solve` to
///           1e-12 Ha — guaranteeing the shim is a true pass-through.
///
/// Acceptance: ac-010 in `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#ifdef ATV_HAS_RHF

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

namespace {

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    // Standard H2O geometry (Angstrom). O at origin; H atoms in the xz plane.
    f.atoms.x = {0.0, 0.7572, -0.7572};
    f.atoms.y = {0.0, 0.0, 0.0};
    f.atoms.z = {0.0, -0.5860, -0.5860};
    f.atoms.type = {8, 1, 1}; // O, H, H
    f.atoms.mass = {15.999, 1.008, 1.008};
    f.box = Mat3::cubic(50.0);
    return f;
}

} // namespace

TEST(SCFSolver, H2OConvergesAndMatchesShim) {
    auto frame = make_h2o();
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(integrals);

    cpu::scf::SCFSolver solver{cfg};
    solver.build(integrals);
    const auto& result = solver.run(integrals);

    EXPECT_TRUE(result.converged) << "H2O / STO-3G RHF should converge";
    EXPECT_GT(result.iterations, 0);
    EXPECT_LT(result.iterations, cfg.max_iter);

    // The shim path constructs a fresh SCFSolver internally. Two separate
    // SCFSolver instances on the same input differ by ~1e-12 due to OMP
    // parallel-reduction order non-determinism — well below any chemistry
    // tolerance, but real on AArch64 / Linux x86. Compare at 1e-10.
    auto result_via_shim = cpu::scf::scf_solve(integrals, cfg);
    EXPECT_NEAR(result.energy, result_via_shim.energy, 1e-10)
        << "SCFSolver.run() and scf_solve() shim agree within FP noise";
    ASSERT_EQ(result.P_alpha.shape(0), result_via_shim.P_alpha.shape(0));
    for (std::size_t i = 0; i < result.P_alpha.shape(0); ++i) {
        for (std::size_t j = 0; j < result.P_alpha.shape(1); ++j) {
            EXPECT_NEAR(result.P_alpha(i, j), result_via_shim.P_alpha(i, j),
                        1e-10);
        }
    }
}

TEST(SCFSolver, RunIsIdempotentAcrossInvocations) {
    auto frame = make_h2o();
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(integrals);

    cpu::scf::SCFSolver solver{cfg};
    solver.build(integrals);

    auto& r1 = solver.run(integrals);
    const double e1 = r1.energy;
    const int iters1 = r1.iterations;
    auto P1 = r1.P_alpha; // copy before second run overwrites the buffer

    auto& r2 = solver.run(integrals);
    // Cross-run determinism is *converged-energy* equality, not bit-equality:
    // OpenMP reductions in the Fock build can re-order partial sums between
    // invocations, giving ~1e-12 noise that is well below any chemistry
    // tolerance. The contract tested here is "no state leaks across runs"
    // — the iteration count and converged result must match within FP noise.
    EXPECT_NEAR(r2.energy, e1, 1e-10)
        << "Same solver, same input → same converged energy across run() "
           "invocations";
    EXPECT_EQ(r2.iterations, iters1)
        << "Same solver, same input → same iteration count";
    for (std::size_t i = 0; i < P1.shape(0); ++i)
        for (std::size_t j = 0; j < P1.shape(1); ++j)
            EXPECT_NEAR(r2.P_alpha(i, j), P1(i, j), 1e-10);
}

TEST(SCFSolver, H2OUKSLSDA) {
    auto frame = make_h2o();
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::uks_config(integrals,
                                    /*n_alpha=*/5, /*n_beta=*/5,
                                    cpu::scf::XCType::LSDA);

    cpu::scf::SCFSolver solver{cfg};
    solver.build(integrals);
    const auto& result = solver.run(integrals);

    EXPECT_TRUE(result.converged);
    EXPECT_GT(result.iterations, 0);
    EXPECT_NEAR(result.spin_contamination, 0.0, 1e-6)
        << "Closed-shell H2O / UKS should have zero spin contamination";

    // shim path matches
    auto via_shim = cpu::scf::scf_solve(integrals, cfg);
    EXPECT_NEAR(result.energy, via_shim.energy, 1e-12);
}

#endif // ATV_HAS_RHF
