/// @file test_pw_smearing.cpp
/// @brief Unit tests for metallic occupation smearing (FD / Gaussian / MP-1).
#include "atomiverse/cpu/planewave/internal/smearing.h"
#include "atomiverse/cpu/planewave/pw_solver.h"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using namespace ATV;
using namespace ATV::cpu::pw;
using namespace ATV::cpu::pw::internal;

TEST(PWSmearing, FermiDiracClosedForm) {
    // f(0)=1/2, f(±∞) limits, f(x)+f(-x)=1
    EXPECT_NEAR(static_cast<double>(smear_f_spinless(
                    AccReal(0), SmearKind::fermi_dirac)),
                0.5, 1e-12);
    EXPECT_NEAR(static_cast<double>(smear_f_spinless(
                    AccReal(40), SmearKind::fermi_dirac)),
                0.0, 1e-12);
    EXPECT_NEAR(static_cast<double>(smear_f_spinless(
                    AccReal(-40), SmearKind::fermi_dirac)),
                1.0, 1e-12);
    const AccReal x = AccReal(1.2);
    const AccReal fp =
        smear_f_spinless(x, SmearKind::fermi_dirac);
    const AccReal fm =
        smear_f_spinless(-x, SmearKind::fermi_dirac);
    EXPECT_NEAR(static_cast<double>(fp + fm), 1.0, 1e-12);
    // Band occupation is 2× spinless
    EXPECT_NEAR(static_cast<double>(smear_f_band(
                    AccReal(0), AccReal(0), AccReal(0.02),
                    SmearKind::fermi_dirac)),
                1.0, 1e-12);
}

TEST(PWSmearing, GaussianClosedForm) {
    // f(0)=1/2, f(x)+f(-x)=1 for Gaussian/erfc
    EXPECT_NEAR(static_cast<double>(smear_f_spinless(
                    AccReal(0), SmearKind::gaussian)),
                0.5, 1e-12);
    const AccReal x = AccReal(0.7);
    const AccReal fp = smear_f_spinless(x, SmearKind::gaussian);
    const AccReal fm = smear_f_spinless(-x, SmearKind::gaussian);
    EXPECT_NEAR(static_cast<double>(fp + fm), 1.0, 1e-12);
}

TEST(PWSmearing, MethfesselPaxtonOrder1) {
    // MP-1 at x=0: 1/2 erfc(0) + A1 H1(0) exp(0) = 1/2  (H1(0)=0)
    EXPECT_NEAR(static_cast<double>(smear_f_spinless(
                    AccReal(0), SmearKind::methfessel_paxton, 1)),
                0.5, 1e-12);
    // MP-1: f(x) = 1/2 erfc(x) - (x/(2√π)) e^{-x²}
    const double x = 0.5;
    const double sqrt_pi = std::sqrt(4.0 * std::atan(1.0));
    const double ref = 0.5 * std::erfc(x)
                       - (x / (2.0 * sqrt_pi)) * std::exp(-x * x);
    EXPECT_NEAR(static_cast<double>(smear_f_spinless(
                    AccReal(x), SmearKind::methfessel_paxton, 1)),
                ref, 1e-12);
}

TEST(PWSmearing, FindFermiRecoversNelec) {
    // Three levels at -1, 0, +1 Ha; want Ne=3 with FD degauss=0.1
    // With spin degeneracy, roughly 1.5 bands occupied.
    std::vector<AccReal> evals = {AccReal(-1), AccReal(0), AccReal(1),
                                  AccReal(2)};
    std::vector<AccReal> w(evals.size(), AccReal(1)); // single k, weight 1
    const AccReal ne = AccReal(3);
    const AccReal deg = AccReal(0.1);
    const AccReal mu = find_fermi_level(evals, w, ne, deg,
                                       SmearKind::fermi_dirac, 1);
    auto fill =
        fill_occupations(evals, w, mu, deg, SmearKind::fermi_dirac, 1);
    EXPECT_NEAR(static_cast<double>(fill.n_occ), 3.0, 1e-8);
    // Mid-spectrum levels should be fractional.
    EXPECT_GT(static_cast<double>(fill.occ[1]), 0.01);
    EXPECT_LT(static_cast<double>(fill.occ[1]), 1.99);
}

TEST(PWSmearing, PublicMembersNoConfigBag) {
    PW pw;
    pw.smearing = Smearing::fermi_dirac;
    pw.degauss = 0.02;
    pw.mp_order = 1;
    pw.extra_bands = 4;
    EXPECT_EQ(pw.smearing, Smearing::fermi_dirac);
    EXPECT_DOUBLE_EQ(pw.degauss, 0.02);
    EXPECT_EQ(pw.extra_bands, 4);
    // Result fields default-constructible
    PWResult r;
    EXPECT_DOUBLE_EQ(r.e_fermi, 0.0);
    EXPECT_DOUBLE_EQ(r.n_elec_occ, 0.0);
    EXPECT_DOUBLE_EQ(r.e_entropy, 0.0);
}
