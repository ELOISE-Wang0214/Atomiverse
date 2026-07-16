/// @file tests/cpu/scf/test_fci.cpp
/// @brief Fast unit tests for the tiny full-CI reference solver (fci_energy,
///        fci_hamiltonian). Spec: hf-fci-tiny.
///
/// Structural, millisecond-scale invariants ONLY: the He single-determinant
/// boundary, CI-Hamiltonian Hermiticity, the determinant count, the no-fallback
/// guards, and the variational + correlation-sign inequalities. NO external-
/// engine (PySCF) parity numbers live here -- all such parity checks are
/// validation benchmarks in benchmarks/qc/fci/validation/test_fci_vs_pyscf.cpp.

#ifdef ATV_HAS_RHF

#include <cmath>
#include <cstdint>
#include <stdexcept>

#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/fci.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

namespace {

cpu::Frame make_he() {
    cpu::Frame f(1);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_type)[0] = 2;
    f.span(ATV::kAtom_mass)[0] = 4.0026;
    f.simbox = ATV::cpu::Box::cubic(50.0);
    return f;
}

cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = 0.0;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_z)[1] = 0.7414;
    f.span(ATV::kAtom_type)[0] = 1;
    f.span(ATV::kAtom_type)[1] = 1;
    f.span(ATV::kAtom_mass)[0] = 1.008;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.simbox = ATV::cpu::Box::cubic(50.0);
    return f;
}

// Solve a restricted reference and return the unified SCFResult.
cpu::scf::SCFResult solve_rhf(const cpu::scf::IntegralPack& ints) {
    cpu::scf::SCFSolver solver{cpu::scf::rhf_config(ints, 1e-12)};
    solver.build(ints);
    return solver.run(ints);
}

// Binomial coefficient for the determinant-count invariant.
long long choose(int n, int k) {
    if (k < 0 || k > n)
        return 0;
    long long r = 1;
    for (int i = 0; i < k; ++i)
        r = r * (n - i) / (i + 1);
    return r;
}

} // namespace

// He/STO-3G has a single determinant (n_mo = 1): full CI == RHF, correlation is
// exactly zero, and E_FCI == E_RHF.
TEST(FCI, He_STO3G_SingleDeterminant) {
    auto frame = make_he();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);
    ASSERT_EQ(scf.C_alpha.shape(1), 1u); // no virtuals

    auto res = cpu::scf::fci_energy(ints, scf);
    EXPECT_EQ(res.e_corr, 0.0);
    EXPECT_TRUE(std::isfinite(res.e_corr));
    EXPECT_TRUE(std::isfinite(res.e_fci));
    EXPECT_NEAR(res.e_fci, scf.energy, 1e-10);

    // N_det == C(1,1) * C(1,1) == 1.
    auto hci = cpu::scf::fci_hamiltonian(ints, scf);
    EXPECT_EQ(hci.shape(0), 1u);
    EXPECT_EQ(hci.shape(1), 1u);
}

// H2/STO-3G: N_det == C(n_mo,n_alpha) * C(n_mo,n_beta) == 4, and the CI
// Hamiltonian is symmetric to ~1e-12.
TEST(FCI, H2_STO3G_CountAndHermiticity) {
    auto frame = make_h2();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    const int n_mo = static_cast<int>(scf.C_alpha.shape(1));
    const long long n_det =
        choose(n_mo, scf.n_alpha) * choose(n_mo, scf.n_beta);
    EXPECT_EQ(n_det, 4);

    auto hci = cpu::scf::fci_hamiltonian(ints, scf);
    ASSERT_EQ(hci.shape(0), static_cast<std::size_t>(n_det));
    ASSERT_EQ(hci.shape(1), static_cast<std::size_t>(n_det));

    double max_asym = 0.0;
    for (std::size_t i = 0; i < hci.shape(0); ++i)
        for (std::size_t j = 0; j < hci.shape(1); ++j)
            max_asym = std::max(max_asym, std::abs(hci(i, j) - hci(j, i)));
    EXPECT_LE(max_asym, 1e-12);
}

// H2/STO-3G: variational bound (E_FCI <= E_RHF) and correlation sign (<= 0).
TEST(FCI, H2_STO3G_VariationalAndSign) {
    auto frame = make_h2();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    auto res = cpu::scf::fci_energy(ints, scf);
    EXPECT_TRUE(std::isfinite(res.e_fci));
    EXPECT_LE(res.e_fci, scf.energy + 1e-10); // variational
    EXPECT_LE(res.e_corr, 0.0);               // correlation lowers energy
    EXPECT_NEAR(res.e_corr, res.e_fci - scf.energy, 1e-12);
}

// MP2-vs-FCI two-level cross-check (H2/STO-3G): FCI is exact, so it recovers at
// least as much correlation as MP2 -- e_corr_fci <= e_corr_mp2 < 0 (MP2
// under-binds H2/STO-3G at equilibrium). NEVER assert |E_FCI| >= |E_MP2| via
// a reversed inequality on the (negative) e_corr values.
TEST(FCI, H2_STO3G_Mp2VsFci) {
    auto frame = make_h2();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    const auto mp2 = cpu::scf::mp2_energy(ints, scf);
    const auto fci = cpu::scf::fci_energy(ints, scf);
    EXPECT_LT(fci.e_corr, 0.0);
    // FCI is exact -> recovers at least as much correlation as MP2 (MP2
    // under-binds H2/STO-3G at equilibrium).
    EXPECT_LE(fci.e_corr, mp2.e_corr + 1e-12);
}

// No-fallback guards all throw std::invalid_argument.
TEST(FCI, GuardsThrow) {
    auto frame = make_h2();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    // (1) ERI not materialized (compute_eri=false).
    auto ints_no_eri =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                  /*compute_eri=*/false);
    EXPECT_THROW(cpu::scf::fci_energy(ints_no_eri, scf), std::invalid_argument);

    // (2) SCF result not converged.
    auto scf_unconverged = scf;
    scf_unconverged.converged = false;
    EXPECT_THROW(cpu::scf::fci_energy(ints, scf_unconverged),
                 std::invalid_argument);

    // (3) Open-shell reference (n_alpha != n_beta).
    auto scf_open = scf;
    scf_open.n_beta = scf.n_alpha - 1;
    EXPECT_THROW(cpu::scf::fci_energy(ints, scf_open), std::invalid_argument);

    // (4) Over the TINY cap (n_mo > 14).
    auto scf_big = scf;
    scf_big.C_alpha = xt::zeros<double>({std::size_t(2), std::size_t(20)});
    EXPECT_THROW(cpu::scf::fci_energy(ints, scf_big), std::invalid_argument);
}

#endif // ATV_HAS_RHF
