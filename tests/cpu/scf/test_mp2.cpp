/// @file tests/cpu/scf/test_mp2.cpp
/// @brief Fast unit tests for the AO->MO transform (ao2mo_ovov) and
///        closed-shell MP2 (mp2_energy). Spec: hf-ao2mo-mp2.
///
/// Structural, millisecond-scale invariants only: the He zero-virtual edge,
/// ovov permutation symmetry, and the no-fallback guards. All external-engine
/// (PySCF) parity checks are validation benchmarks, not unit tests, and live in
/// benchmarks/qc/mp2/validation/test_mp2_vs_pyscf.cpp.

#ifdef ATV_HAS_RHF

#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

namespace {

// H2O / spec geometry (Angstrom): O(0,0,0.1173), H(0,+-0.7572,-0.4692).
cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = 0.0;
    f.span(ATV::kAtom_x)[2] = 0.0;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.7572;
    f.span(ATV::kAtom_y)[2] = -0.7572;
    f.span(ATV::kAtom_z)[0] = 0.1173;
    f.span(ATV::kAtom_z)[1] = -0.4692;
    f.span(ATV::kAtom_z)[2] = -0.4692;
    f.span(ATV::kAtom_type)[0] = 8;
    f.span(ATV::kAtom_type)[1] = 1;
    f.span(ATV::kAtom_type)[2] = 1;
    f.span(ATV::kAtom_mass)[0] = 15.999;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.span(ATV::kAtom_mass)[2] = 1.008;
    f.simbox = ATV::cpu::Box::cubic(50.0);
    return f;
}

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

// Solve a restricted reference and return the unified SCFResult.
cpu::scf::SCFResult solve_rhf(const cpu::scf::IntegralPack& ints) {
    cpu::scf::SCFSolver solver{cpu::scf::rhf_config(ints, 1e-12)};
    solver.build(ints);
    return solver.run(ints);
}

} // namespace

// ac-002: He/STO-3G has no virtuals -> correlation is exactly 0.0, never NaN.
TEST(MP2, He_STO3G_ZeroVirtuals) {
    auto frame = make_he();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    auto res = cpu::scf::mp2_energy(ints, scf);
    EXPECT_EQ(res.e_corr, 0.0);
    EXPECT_TRUE(std::isfinite(res.e_corr));
    EXPECT_TRUE(std::isfinite(res.e_tot));
}

// ac-003: ao2mo_ovov satisfies the (ia|jb) = (jb|ia) charge-cloud symmetry.
TEST(MP2, Ao2moOvovSymmetry) {
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    const int n_occ = scf.n_alpha;
    const int n_mo = static_cast<int>(scf.C_alpha.shape(1));
    const int n_vir = n_mo - n_occ;
    ASSERT_GT(n_vir, 0);

    auto ovov = cpu::scf::ao2mo_ovov(ints.ERI, scf.C_alpha, n_occ, n_mo);
    double max_diff = 0.0;
    for (int i = 0; i < n_occ; ++i)
        for (int a = 0; a < n_vir; ++a)
            for (int j = 0; j < n_occ; ++j)
                for (int b = 0; b < n_vir; ++b)
                    max_diff = std::max(max_diff, std::abs(ovov(i, a, j, b)
                                                           - ovov(j, b, i, a)));
    EXPECT_LE(max_diff, 1e-12);
}

// ac-004: no-fallback guards throw std::invalid_argument.
TEST(MP2, GuardsThrow) {
    auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    // (1) ERI not materialized (compute_eri=false).
    auto ints_no_eri =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G, 0,
                                  /*compute_eri=*/false);
    EXPECT_THROW(cpu::scf::mp2_energy(ints_no_eri, scf), std::invalid_argument);

    // (2) SCF result not converged.
    auto scf_unconverged = scf;
    scf_unconverged.converged = false;
    EXPECT_THROW(cpu::scf::mp2_energy(ints, scf_unconverged),
                 std::invalid_argument);

    // (3) Unrestricted reference (n_alpha != n_beta).
    auto scf_open = scf;
    scf_open.n_beta = scf.n_alpha - 1;
    EXPECT_THROW(cpu::scf::mp2_energy(ints, scf_open), std::invalid_argument);
}

#endif // ATV_HAS_RHF
