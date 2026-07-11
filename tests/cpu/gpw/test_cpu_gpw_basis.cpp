/// @file tests/cpu/gpw/test_cpu_gpw_basis.cpp
/// @brief Acceptance ac-001 (gpw spec, D1): CP2K Gaussian basis parser parses
///        DZVP-GTH-PADE (Ar) and DZVP-MOLOPT-GTH (O, H) exactly.
///
/// Reference numbers are read by hand from the CP2K data files at
/// tests/QS/regtest-gpw-1/Ar.inp); O/H from BASIS_MOLOPT (the H2O tutorial).

#include <string>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/basis_data/gpw_data.h"

using ATV::cpu::gpw::GaussianBasis;
using ATV::cpu::gpw::builtin_gth_basis;

namespace {
} // namespace

// Ar DZVP-GTH-PADE: two sets. Set 0 is n=3, l=0..1, 4 exps, 2 s + 2 p shells
// (4 coef columns ordered s1,s2,p1,p2). Set 1 is a single d shell.
TEST(GthBasisParser, ArgonDzvpGthPade) {
    GaussianBasis b = builtin_gth_basis("Ar", "DZVP-GTH-PADE");
    ASSERT_EQ(b.sets.size(), 2u);

    const auto& s0 = b.sets[0];
    EXPECT_EQ(s0.n, 3);
    EXPECT_EQ(s0.lmin, 0);
    EXPECT_EQ(s0.lmax, 1);
    EXPECT_EQ(s0.nexp(), 4);
    ASSERT_EQ(s0.nshell.size(), 2u);
    EXPECT_EQ(s0.nshell_l(0), 2);
    EXPECT_EQ(s0.nshell_l(1), 2);
    EXPECT_EQ(s0.ncoef, 4);
    EXPECT_DOUBLE_EQ(s0.exps[0], 2.6436290990);
    EXPECT_DOUBLE_EQ(s0.exps[3], 0.1663306775);
    // row 0: s1, s2, p1, p2
    EXPECT_DOUBLE_EQ(s0.coef_at(0, 0), 0.3737616862);
    EXPECT_DOUBLE_EQ(s0.coef_at(0, 1), 0.0000000000);
    EXPECT_DOUBLE_EQ(s0.coef_at(0, 2), 0.1735201807);
    EXPECT_DOUBLE_EQ(s0.coef_at(0, 3), 0.0000000000);
    // row 3
    EXPECT_DOUBLE_EQ(s0.coef_at(3, 0), -0.3339930421);
    EXPECT_DOUBLE_EQ(s0.coef_at(3, 1), 1.0000000000);
    EXPECT_DOUBLE_EQ(s0.coef_at(3, 3), 1.0000000000);

    const auto& s1 = b.sets[1];
    EXPECT_EQ(s1.lmin, 2);
    EXPECT_EQ(s1.lmax, 2);
    EXPECT_EQ(s1.nexp(), 1);
    EXPECT_EQ(s1.ncoef, 1);
    EXPECT_DOUBLE_EQ(s1.exps[0], 0.8500000000);
    EXPECT_DOUBLE_EQ(s1.coef_at(0, 0), 1.0000000000);
}

// O DZVP-MOLOPT-GTH: single shared-exponent set, n=2, l=0..2, 7 exps,
// 2 s + 2 p + 1 d shells (5 coef columns).
TEST(GthBasisParser, OxygenDzvpMoloptGth) {
    GaussianBasis b = builtin_gth_basis("O", "DZVP-MOLOPT-GTH");
    ASSERT_EQ(b.sets.size(), 1u);
    const auto& s = b.sets[0];
    EXPECT_EQ(s.n, 2);
    EXPECT_EQ(s.lmin, 0);
    EXPECT_EQ(s.lmax, 2);
    EXPECT_EQ(s.nexp(), 7);
    ASSERT_EQ(s.nshell.size(), 3u);
    EXPECT_EQ(s.nshell_l(0), 2);
    EXPECT_EQ(s.nshell_l(1), 2);
    EXPECT_EQ(s.nshell_l(2), 1);
    EXPECT_EQ(s.ncoef, 5);
    EXPECT_DOUBLE_EQ(s.exps[0], 12.015954705512);
    EXPECT_DOUBLE_EQ(s.exps[6], 0.046760918300);
    EXPECT_DOUBLE_EQ(s.coef_at(0, 0), -0.060190841200);
    EXPECT_DOUBLE_EQ(s.coef_at(0, 4), 0.014807054400);
    EXPECT_DOUBLE_EQ(s.coef_at(6, 0), -0.000255945800);
    EXPECT_DOUBLE_EQ(s.coef_at(6, 4), -0.005771736600);
}

// H DZVP-MOLOPT-GTH: single set, n=2, l=0..1, 7 exps, 2 s + 1 p (3 coef cols).
TEST(GthBasisParser, HydrogenDzvpMoloptGth) {
    GaussianBasis b = builtin_gth_basis("H", "DZVP-MOLOPT-GTH");
    ASSERT_EQ(b.sets.size(), 1u);
    const auto& s = b.sets[0];
    EXPECT_EQ(s.lmin, 0);
    EXPECT_EQ(s.lmax, 1);
    EXPECT_EQ(s.nexp(), 7);
    EXPECT_EQ(s.ncoef, 3);
    EXPECT_EQ(s.nshell_l(0), 2);
    EXPECT_EQ(s.nshell_l(1), 1);
    EXPECT_DOUBLE_EQ(s.exps[0], 11.478000339908);
    EXPECT_DOUBLE_EQ(s.coef_at(0, 0), 0.024916243200);
    EXPECT_DOUBLE_EQ(s.coef_at(0, 2), 0.024510918200);
    EXPECT_DOUBLE_EQ(s.coef_at(6, 1), -0.059693171300);
}

TEST(GthBasisParser, MissingEntryThrows) {
    EXPECT_THROW(builtin_gth_basis("Ar", "NO-SUCH-BASIS"), std::runtime_error);
}
