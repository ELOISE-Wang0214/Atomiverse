/// @file tests/cpu/scf/test_level_shift.cpp
/// @brief Solver-level tests for the virtual-orbital level shift
///        (spec hf-robust-03-level-shift, ac-001/ac-003/ac-004/ac-005).
///
/// PySCF 2.13.1 references (cart=True, conv_tol=1e-12), regenerated
/// 2026-07-02 for the exact geometries below:
///   E(H2O/6-31G RHF)          = -75.983974472722 Ha
///   |E(λ=0.5) − E(λ=0)|       = 2.8e-14 Ha   (spec bar: 1e-9)
///   max|ε(λ=0.5) − ε(λ=0)|    = 2.5e-7 Ha    (spec bar: 1e-6)
///
/// Hard case (tutorial ch02 "SCF convergence is a numerical problem"):
/// stretched H2O (both O–H = 1.7 Å), 6-31G RHF, Hcore guess, no DIIS,
/// zero damping — the bare iteration oscillates past max_iter = 100
/// (converged == false), while λ = 0.5 converges in 57 iterations to the
/// same energy a DIIS reference run reaches. All numbers quoted in the
/// tutorial section are the constants asserted here (lockstep rule).

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

// PySCF 2.13.1 anchor for H2O/6-31G RHF at the geometry below (ac-001).
constexpr double kPyscfH2O631gAnchor = -75.9839744727;

// Rescue shift and its iteration count (pinned tutorial numbers; ac-005).
// Measured scan on the stretched system (single-thread probe):
//   λ=0.0 → no convergence in 100 iters (oscillates ~4.5 Ha too high)
//   λ=0.3 → 52 iters   λ=0.5 → 57 iters
//   λ=1.0 → 105 iters  λ=2.0 → 145 iters (over-damped)
// λ=0.5 is the tutorial choice: large enough to stabilize, small enough
// not to crawl.
constexpr double kRescueLambda = 0.5; // Hartree
constexpr int kRescueIterations = 57;

/// H2O / 6-31G anchor geometry (Å): O (0,0,0.1173); H (0,±0.7572,−0.4692).
/// `stretch` scales both O–H displacement vectors about O
/// (1.0 = equilibrium; 1.7 Å bonds at 1.7/|r_OH| ≈ 1.775).
cpu::Frame make_h2o(double stretch) {
    cpu::Frame f(3);
    const double oy = 0.0, oz = 0.1173;
    const double hy = 0.7572, hz = -0.4692;
    f.atoms.x = {0.0, 0.0, 0.0};
    f.atoms.y = {oy, oy + stretch * (hy - oy), oy - stretch * (hy - oy)};
    f.atoms.z = {oz, oz + stretch * (hz - oz), oz + stretch * (hz - oz)};
    f.atoms.type = {8, 1, 1};
    f.atoms.mass = {15.999, 1.008, 1.008};
    f.box = Mat3::cubic(50.0);
    return f;
}

/// Equilibrium O–H distance of the anchor geometry [Å].
double equilibrium_oh() {
    return std::sqrt(0.7572 * 0.7572 + (0.1173 + 0.4692) * (0.1173 + 0.4692));
}

cpu::scf::SCFConfig rhf_631g_config() {
    cpu::scf::SCFConfig c;
    c.spin = cpu::scf::SpinType::Restricted;
    c.method = cpu::scf::MeanFieldType::HF;
    c.tol = 1e-9;
    c.damping_factor = 0.0;
    c.init_guess = cpu::scf::InitGuess::MINAO;
    return c;
}

/// The tutorial hard-case config: bare fixed-point iteration — Hcore
/// guess, no DIIS, no damping, fixed iteration budget.
cpu::scf::SCFConfig bare_config() {
    cpu::scf::SCFConfig c;
    c.spin = cpu::scf::SpinType::Restricted;
    c.method = cpu::scf::MeanFieldType::HF;
    c.tol = 1e-9;
    c.use_diis = false;
    c.damping_factor = 0.0;
    c.init_guess = cpu::scf::InitGuess::Hcore;
    c.max_iter = 100;
    return c;
}

} // namespace

// ---------------------------------------------------------------------------
// ac-001 — λ=0 default anchors; ac-003/ac-004 — λ-invariance + de-shift.
// ---------------------------------------------------------------------------

// ac-001: default level_shift = 0.0; PySCF anchor to 1e-6 Ha.
TEST(LevelShiftSolver, H2OAnchorMatchesPyscfAtLambdaZero) {
    auto frame = make_h2o(1.0);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);
    auto cfg = rhf_631g_config();
    cfg.n_alpha = ints.n_occ();
    cfg.n_beta = ints.n_occ();

    const auto result = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(result.converged);
    EXPECT_NEAR(result.energy, kPyscfH2O631gAnchor, 1e-6);
}

// ac-003: converged energy invariant under λ = 0.5 (< 1e-9 Ha);
// ac-004: reported ε shift-free after the post-convergence λ=0
// diagonalization (max over ALL orbitals ≤ 1e-6 Ha).
// tutorial ch02 "SCF convergence is a numerical problem" worked numbers;
// regression anchor test_cpu_scf_level_shift.
TEST(LevelShiftSolver, EnergyAndReportedEpsilonInvariantUnderShift) {
    auto frame = make_h2o(1.0);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);

    auto cfg0 = rhf_631g_config();
    cfg0.n_alpha = ints.n_occ();
    cfg0.n_beta = ints.n_occ();
    auto cfg5 = cfg0;
    cfg5.level_shift = 0.5; // Hartree

    const auto r0 = cpu::scf::scf_solve(ints, cfg0);
    const auto r5 = cpu::scf::scf_solve(ints, cfg5);
    ASSERT_TRUE(r0.converged);
    ASSERT_TRUE(r5.converged);

    EXPECT_LT(std::abs(r5.energy - r0.energy), 1e-9)
        << "converged energy must be λ-invariant (PySCF measured 2.8e-14)";

    ASSERT_EQ(r5.eigenvalues_alpha.size(), r0.eigenvalues_alpha.size());
    double max_deps = 0.0;
    for (std::size_t i = 0; i < r0.eigenvalues_alpha.size(); ++i)
        max_deps = std::max(max_deps,
                            std::abs(r5.eigenvalues_alpha(i)
                                     - r0.eigenvalues_alpha(i)));
    EXPECT_LE(max_deps, 1e-6)
        << "virtual ε must be de-shifted after convergence "
           "(PySCF measured 2.5e-7)";
}

// ---------------------------------------------------------------------------
// ac-005 — hard case: λ rescues a divergent stretched-H2O SCF.
// ---------------------------------------------------------------------------

// tutorial ch02 "SCF convergence is a numerical problem" hard case;
// regression anchor test_cpu_scf_level_shift (ac-005, ac-007).
TEST(LevelShiftSolver, LevelShiftRescuesStretchedH2O) {
    const double stretch = 1.7 / equilibrium_oh(); // both O–H → 1.7 Å
    auto frame = make_h2o(stretch);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);

    // λ = 0: the bare iteration oscillates and never converges.
    auto cfg_bare = bare_config();
    cfg_bare.n_alpha = ints.n_occ();
    cfg_bare.n_beta = ints.n_occ();
    const auto r_bare = cpu::scf::scf_solve(ints, cfg_bare);
    EXPECT_FALSE(r_bare.converged)
        << "bare Hcore/no-DIIS/no-damping iteration must NOT converge "
           "within 100 iterations (that is what makes the case hard)";

    // Same settings + λ = 0.5 Ha: converges.
    auto cfg_shift = cfg_bare;
    cfg_shift.level_shift = kRescueLambda;
    const auto r_shift = cpu::scf::scf_solve(ints, cfg_shift);
    ASSERT_TRUE(r_shift.converged)
        << "λ = 0.5 must rescue the divergent iteration";
    EXPECT_EQ(r_shift.iterations, kRescueIterations)
        << "tutorial-pinned iteration count";

    // DIIS reference on the same geometry: energies agree < 1e-8 Ha.
    auto cfg_diis = rhf_631g_config();
    cfg_diis.n_alpha = ints.n_occ();
    cfg_diis.n_beta = ints.n_occ();
    const auto r_diis = cpu::scf::scf_solve(ints, cfg_diis);
    ASSERT_TRUE(r_diis.converged);
    EXPECT_LT(std::abs(r_shift.energy - r_diis.energy), 1e-8)
        << "level-shifted solution must reach the DIIS reference energy";
}

#endif // ATV_HAS_RHF
