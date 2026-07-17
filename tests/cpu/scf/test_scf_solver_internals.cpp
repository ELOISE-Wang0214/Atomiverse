/// @file tests/cpu/scf/test_scf_solver_internals.cpp
/// @brief Phase 3 — direct unit tests of `cpu::scf::internal::*` helpers.
///
/// These verify each internal helper in isolation without going through
/// the full SCF cycle. Tests are allowed to include `cpu/scf/internal/*.h`
/// directly — facade-public *production* headers must not (see
/// `test_scf_solver_internal_isolation.cpp`).
///
/// Acceptance: ac-012 (first half) in
/// `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#ifdef ATV_HAS_RHF

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/internal/density_builder.h"
#include "atomiverse/cpu/scf/internal/dual_convergence_check.h"
#include "atomiverse/cpu/scf/internal/eigen_solver.h"

using namespace ATV::cpu::scf;
using namespace ATV::cpu::scf::internal;

// ── DensityBuilder ──────────────────────────────────────────────────

TEST(SCFSolverInternals, DensityBuilderRestrictedFactor2) {
    // C with 1 occupied orbital: c0 = (1, 0)^T.
    RealMatrix C = xt::zeros<double>({2, 2});
    C(0, 0) = 1.0;
    C(0, 1) = 0.0;
    DensityBuilder db;
    db.prepare(2);
    auto P = db.step(C, /*n_occ=*/1, /*n_basis=*/2, /*factor=*/2.0);
    EXPECT_NEAR(P(0, 0), 2.0, 1e-12);
    EXPECT_NEAR(P(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(P(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(P(1, 1), 0.0, 1e-12);
}

TEST(SCFSolverInternals, DensityBuilderEmptyOccupation) {
    RealMatrix C = xt::eye<double>(3);
    DensityBuilder db;
    auto P = db.step(C, /*n_occ=*/0, /*n_basis=*/3);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_NEAR(P(i, j), 0.0, 1e-12);
}

// ── EigenSolver ─────────────────────────────────────────────────────

TEST(SCFSolverInternals, EigenSolverGeneralizedEig) {
    // F = [[2, 0], [0, 1]]; S = I → eigenvalues {1, 2}.
    RealMatrix F = xt::eye<double>(2);
    F(0, 0) = 2.0;
    F(1, 1) = 1.0;
    RealMatrix S = xt::eye<double>(2);
    EigenSolver es;
    auto [eig, C] = es.step(F, S);
    EXPECT_NEAR(eig(0), 1.0, 1e-12);
    EXPECT_NEAR(eig(1), 2.0, 1e-12);
}

// ── DualConvergenceCheck ────────────────────────────────────────────

TEST(SCFSolverInternals, DualConvergenceFirstIterAlwaysFalse) {
    DualConvergenceCheck cc{1e-9, 1e-4};
    RealMatrix err = xt::zeros<double>({2, 2}); // tiny error
    auto v = cc.step(/*E=*/-1.0, err);
    EXPECT_FALSE(v.converged) << "no prior energy on iter 0 → never converged";
    EXPECT_TRUE(std::isinf(v.dE));
}

TEST(SCFSolverInternals, DualConvergenceTriggersWhenBothCriteriaMet) {
    DualConvergenceCheck cc{1e-9, 1e-4};
    RealMatrix big_err = xt::ones<double>({2, 2});   // ‖err‖ = 2 ≥ 1e-4
    RealMatrix tiny_err = xt::zeros<double>({2, 2}); // ‖err‖ = 0 < 1e-4

    // Iter 0: prime prev_energy, no convergence reported.
    cc.step(-1.0, big_err);
    // Iter 1: dE ≈ 0 < 1e-9, but err norm 2 → not converged.
    auto v1 = cc.step(-1.0, big_err);
    EXPECT_FALSE(v1.converged);
    // Iter 2: dE ≈ 0 < 1e-9, err norm 0 → converged.
    auto v2 = cc.step(-1.0, tiny_err);
    EXPECT_TRUE(v2.converged);
    EXPECT_LT(std::abs(v2.dE), 1e-9);
}

TEST(SCFSolverInternals, DualConvergenceResetClearsHistory) {
    DualConvergenceCheck cc{1e-9, 1e-4};
    RealMatrix err = xt::zeros<double>({2, 2});
    cc.step(-1.0, err); // primes prev_energy
    cc.reset();
    // After reset, iter 0 must again return ∞ dE.
    auto v = cc.step(-1.0, err);
    EXPECT_TRUE(std::isinf(v.dE));
    EXPECT_FALSE(v.converged);
}

#endif // ATV_HAS_RHF
