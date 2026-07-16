/// @file tests/cpu/scf/test_spin_resolution.cpp
/// @brief Tests for spec `hf-robust-01-charge-multiplicity`: the
///        (N_electrons, multiplicity) -> (n_alpha, n_beta) resolver struct
///        `cpu::scf::SpinResolution` (validating constructor — direct
///        instantiation, no factory function per repo convention) and the
///        `SCFSolver::run` entry guard that rejects inconsistent occupation
///        configurations with `std::invalid_argument`.
///
/// Test-category map (mol harness):
///   - Happy path:        SpinResolution.WorkedTableMatchesTutorialCh03,
///                        SpinResolution.OhMinusChargeFoldsIntoElectronCount
///                        (ac-001)
///   - Edge cases:        SpinResolution.RejectsInvalidCombinations,
///                        SpinResolution.RejectsDegenerateInputs (ac-002),
///                        and each guard condition triggered independently
///                        (SCFSolverGuard.*)
///   - Domain validation: the accept/reject table below IS the domain
///                        validation — it is the PySCF 2.13.1-verified
///                        (charge, M) -> (n_alpha, n_beta) table from the
///                        spec's Domain basis (same accept, same reject;
///                        exception TYPE asserted, never message text).
///   - Immutability:      SpinResolution.IsPureAndDeterministic —
///                        construction is a pure computation over two
///                        by-value ints; two constructions from the same
///                        arguments give identical members.
///
/// Acceptance: ac-001..ac-005, ac-008 in
/// `.claude/specs/hf-robust-01-charge-multiplicity.acceptance.md`.

#ifdef ATV_HAS_RHF

#include <stdexcept>

#include <gtest/gtest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

namespace {

// Pre-change baseline: single-thread (omp_set_num_threads(1)) run of
// the hf-robust-01 probe at base 2c98906, H2O / STO-3G rhf_config, at
// 17 significant digits (ac-004c).
constexpr double kH2ORhfSto3gPreChange = -74.962964802040815;

// Baseline agreement bar [Ha]. Far above the observed cross-architecture drift
// (~1.4e-12 Ha here, x86_64 vs aarch64 BLAS reassociation) and far below
// anything a real regression in the guard would move.
constexpr double kEnergyTol = 1e-8;

/// Standard H2O geometry (Angstrom) — same as test_scf_solver.cpp.
cpu::Frame make_h2o() {
    cpu::Frame f(3);
    // O at origin; H atoms in the xz plane.
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = 0.7572;
    f.span(ATV::kAtom_x)[2] = -0.7572;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.0;
    f.span(ATV::kAtom_y)[2] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_z)[1] = -0.5860;
    f.span(ATV::kAtom_z)[2] = -0.5860;
    f.span(ATV::kAtom_type)[0] = 8;
    f.span(ATV::kAtom_type)[1] = 1;
    f.span(ATV::kAtom_type)[2] = 1; // O, H, H
    f.span(ATV::kAtom_mass)[0] = 15.999;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.span(ATV::kAtom_mass)[2] = 1.008;
    f.simbox = ATV::cpu::Box::cubic(50.0);
    return f;
}

/// Hydroxide geometry (Angstrom): O at origin, H at (0, 0, 0.97).
cpu::Frame make_oh() {
    cpu::Frame f(2);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = 0.0;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_z)[1] = 0.97;
    f.span(ATV::kAtom_type)[0] = 8;
    f.span(ATV::kAtom_type)[1] = 1; // O, H
    f.span(ATV::kAtom_mass)[0] = 15.999;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.simbox = ATV::cpu::Box::cubic(50.0);
    return f;
}

/// Single oxygen atom.
cpu::Frame make_o_atom() {
    cpu::Frame f(1);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_type)[0] = 8; // O
    f.span(ATV::kAtom_mass)[0] = 15.999;
    f.simbox = ATV::cpu::Box::cubic(50.0);
    return f;
}

} // namespace

// ---------------------------------------------------------------------------
// ac-001 — resolver worked table (happy path + domain validation).
// ---------------------------------------------------------------------------

// tutorial ch03 "Charge and multiplicity conventions" worked example;
// regression anchor test_cpu_scf_spin_resolution (ac-001, ac-008).
// Domain validation: this accept table is verified against PySCF 2.13.1
// (mol.spin = 2S = n_alpha - n_beta convention; see spec Domain basis).
TEST(SpinResolution, WorkedTableMatchesTutorialCh03) {
    // H2O, charge 0, M = 1  ->  (5, 5)
    {
        const cpu::scf::SpinResolution occ(10, 1);
        EXPECT_EQ(occ.n_alpha, 5);
        EXPECT_EQ(occ.n_beta, 5);
    }
    // H2O+, charge +1, M = 2  ->  (5, 4)
    {
        const cpu::scf::SpinResolution occ(9, 2);
        EXPECT_EQ(occ.n_alpha, 5);
        EXPECT_EQ(occ.n_beta, 4);
    }
    // O atom, charge 0, M = 3 (triplet ground state)  ->  (5, 3)
    {
        const cpu::scf::SpinResolution occ(8, 3);
        EXPECT_EQ(occ.n_alpha, 5);
        EXPECT_EQ(occ.n_beta, 3);
    }
}

// tutorial ch03 "Charge and multiplicity conventions" worked example;
// regression anchor test_cpu_scf_spin_resolution (ac-001, ac-008).
// OH-, charge -1, M = 1 -> (5, 5): the molecular charge is folded into
// IntegralPack::n_electrons by build_integrals, so the resolver only ever
// sees (N_e, M).
TEST(SpinResolution, OhMinusChargeFoldsIntoElectronCount) {
    auto frame = make_oh();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G,
                                          /*charge=*/-1);
    ASSERT_EQ(ints.n_electrons, 10)
        << "OH-: Z(O) + Z(H) - charge = 8 + 1 - (-1) = 10";

    const cpu::scf::SpinResolution occ(ints.n_electrons, 1);
    EXPECT_EQ(occ.n_alpha, 5);
    EXPECT_EQ(occ.n_beta, 5);
}

// ---------------------------------------------------------------------------
// ac-002 — resolver rejection table (edge cases + domain validation).
// ---------------------------------------------------------------------------

// tutorial ch03 "Charge and multiplicity conventions" worked example
// (rejected combinations); regression anchor test_cpu_scf_spin_resolution
// (ac-002, ac-008). Domain validation: PySCF 2.13.1 rejects the same three
// combinations (parity -> RuntimeError, n_beta<0 -> AssertionError); we
// assert only our exception TYPE, never any message text.
TEST(SpinResolution, RejectsInvalidCombinations) {
    // H2O(0, 2): N_e = 10 even, M = 2 -> parity violation (n_alpha/n_beta
    // would be non-integer).
    EXPECT_THROW(cpu::scf::SpinResolution(10, 2), std::invalid_argument);
    // H2O(+1, 1): N_e = 9 odd, M = 1 -> parity violation.
    EXPECT_THROW(cpu::scf::SpinResolution(9, 1), std::invalid_argument);
    // H(0, 4): N_e = 1, M = 4 -> n_beta = (1 - 4 + 1)/2 = -1 < 0.
    EXPECT_THROW(cpu::scf::SpinResolution(1, 4), std::invalid_argument);
}

// ac-002 (degenerate inputs): multiplicity < 1 and n_electrons < 0 are
// rejected with the same exception type.
TEST(SpinResolution, RejectsDegenerateInputs) {
    // multiplicity < 1 (M = 2S + 1 >= 1 always).
    EXPECT_THROW(cpu::scf::SpinResolution(10, 0), std::invalid_argument);
    // negative electron count.
    EXPECT_THROW(cpu::scf::SpinResolution(-1, 1), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Immutability — SpinResolution construction is pure.
// ---------------------------------------------------------------------------

// Constructing the resolver twice with identical arguments must produce
// identical members (no hidden state). Input-mutation assertions are N/A:
// both parameters are passed by value as plain ints, so the caller's values
// cannot be modified by construction.
TEST(SpinResolution, IsPureAndDeterministic) {
    const cpu::scf::SpinResolution first(9, 2);
    const cpu::scf::SpinResolution second(9, 2);
    EXPECT_EQ(first.n_alpha, second.n_alpha);
    EXPECT_EQ(first.n_beta, second.n_beta);
}

// ---------------------------------------------------------------------------
// ac-005 — SpinResolution composes with the existing raw-count config
// builders (no new factory overloads; direct instantiation per repo rule).
// ---------------------------------------------------------------------------

// tutorial ch03 "Charge and multiplicity conventions" worked example
// (SpinResolution + uhf_config usage); regression anchor
// test_cpu_scf_spin_resolution (ac-005, ac-008).
TEST(SpinResolution, FeedsUhfConfigRawCounts) {
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G,
                                          /*charge=*/+1);
    ASSERT_EQ(ints.n_electrons, 9);

    // H2O+ doublet: resolve, then feed the existing raw-count builder.
    const cpu::scf::SpinResolution occ(ints.n_electrons, /*multiplicity=*/2);
    const auto cfg = cpu::scf::uhf_config(ints, occ.n_alpha, occ.n_beta);
    EXPECT_EQ(cfg.n_alpha, 5);
    EXPECT_EQ(cfg.n_beta, 4);
    EXPECT_EQ(cfg.spin, cpu::scf::SpinType::Unrestricted);
    EXPECT_EQ(cfg.method, cpu::scf::MeanFieldType::HF);
}

// tutorial ch03 "Charge and multiplicity conventions" worked example
// (SpinResolution + uks_config usage, O-atom triplet); regression anchor
// test_cpu_scf_spin_resolution (ac-005, ac-008).
TEST(SpinResolution, FeedsUksConfigRawCounts) {
    auto frame = make_o_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G,
                                          /*charge=*/0);
    ASSERT_EQ(ints.n_electrons, 8);

    const cpu::scf::SpinResolution occ(ints.n_electrons, /*multiplicity=*/3);
    const auto cfg = cpu::scf::uks_config(ints, occ.n_alpha, occ.n_beta);
    EXPECT_EQ(cfg.n_alpha, 5);
    EXPECT_EQ(cfg.n_beta, 3);
    EXPECT_EQ(cfg.spin, cpu::scf::SpinType::Unrestricted);
    EXPECT_EQ(cfg.method, cpu::scf::MeanFieldType::KS);
}

// ---------------------------------------------------------------------------
// ac-003 / ac-004 — SCFSolver::run entry guard.
// ---------------------------------------------------------------------------

// ac-003: the core defect fix. rhf_config on an odd-electron pack silently
// truncates (n_occ() = n_electrons / 2 -> n_alpha + n_beta = N_e - 1); the
// solver entry guard must now reject it before any SCF iteration runs.
// tutorial ch02 "Closed-shell validity" worked example (H2O+ via
// rhf_config/SCFSolver::run throws); regression anchor
// test_cpu_scf_spin_resolution (ac-003, ac-008).
TEST(SCFSolverGuard, OddElectronRhfThrows) {
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G,
                                          /*charge=*/+1);
    ASSERT_EQ(ints.n_electrons, 9) << "H2O+: 10 - 1 = 9 electrons";

    auto cfg = cpu::scf::rhf_config(ints, 1e-10);
    cfg.max_iter = 100;
    cpu::scf::SCFSolver solver{cfg};
    solver.build(ints);
    EXPECT_THROW(solver.run(ints), std::invalid_argument);
}

// ac-004(a): guard condition 2 — n_alpha + n_beta != integrals.n_electrons.
// Neutral H2O has 10 electrons; uhf_config(ints, 5, 4) sums to 9.
TEST(SCFSolverGuard, OccupationSumMismatchThrows) {
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    ASSERT_EQ(ints.n_electrons, 10);

    const auto cfg = cpu::scf::uhf_config(ints, 5, 4); // sum = 9 != 10
    cpu::scf::SCFSolver solver{cfg};
    solver.build(ints);
    EXPECT_THROW(solver.run(ints), std::invalid_argument);
}

// Guard condition 1 triggered independently: negative occupation count with
// a *correct* sum (n_alpha + n_beta = -1 + 11 = 10 = n_electrons), so only
// the negativity check can fire.
TEST(SCFSolverGuard, NegativeOccupationThrows) {
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    ASSERT_EQ(ints.n_electrons, 10);

    auto cfg = cpu::scf::uhf_config(ints, 5, 5);
    cfg.n_alpha = -1;
    cfg.n_beta = 11;
    cpu::scf::SCFSolver solver{cfg};
    solver.build(ints);
    EXPECT_THROW(solver.run(ints), std::invalid_argument);
}

// ac-004(b): guard condition 3 triggered independently — Restricted spin
// with n_alpha != n_beta. Using the H2O+ pack (n_electrons = 9) with
// n_alpha = 5, n_beta = 4 keeps the sum correct (5 + 4 = 9), so ONLY the
// restricted-inequality condition can fire.
TEST(SCFSolverGuard, RestrictedOpenShellThrows) {
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G,
                                          /*charge=*/+1);
    ASSERT_EQ(ints.n_electrons, 9);

    auto cfg = cpu::scf::rhf_config(ints); // spin = Restricted, method = HF
    cfg.n_alpha = 5;
    cfg.n_beta = 4;
    cpu::scf::SCFSolver solver{cfg};
    solver.build(ints);
    EXPECT_THROW(solver.run(ints), std::invalid_argument);
}

// ac-004(c): the guard must have ZERO behavioural impact on valid runs —
// neutral H2O / STO-3G RHF converges to the pre-change baseline within
// kEnergyTol. Bit-identity is NOT asserted: it does not survive a change of
// architecture (x86_64 vs aarch64 reassociate the BLAS reductions and differ
// by ~1e-12 Ha here), while the SCF itself only converges to 1e-9. Single OMP
// thread still pins the reduction order within one machine.
TEST(SCFSolverGuard, NeutralH2ORhfEnergyMatchesBaseline) {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    ASSERT_EQ(ints.n_electrons, 10);

    cpu::scf::SCFSolver solver{cpu::scf::rhf_config(ints)};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged) << "H2O / STO-3G RHF should converge";
    EXPECT_NEAR(result.energy, kH2ORhfSto3gPreChange, kEnergyTol)
        << "Guard must not perturb the valid closed-shell path "
           "(pre-change single-thread baseline)";
}

#endif // ATV_HAS_RHF
