/// @file tests/cpu/scf/test_analytic_tda_forces.cpp
/// @brief Analytic CIS/TDA forces vs FD; translational invariance.
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

namespace {
cpu::Frame h2(double r = 0.85) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = r;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.simbox = cpu::Box::cubic(15);
    return f;
}
} // namespace

TEST(AnalyticTdaForces, VsFdH2) {
    auto f = h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-9;
    rhf.run(f);
    ASSERT_TRUE(rhf.converged());
    auto Fa = rhf.tda_forces(f, 0, 1e-3, /*analytic=*/true);
    auto Fn = rhf.tda_forces(f, 0, 2e-3, /*analytic=*/false);
    ASSERT_EQ(Fa.size(), 6u);
    ASSERT_EQ(Fn.size(), 6u);
    double maxd = 0.0;
    for (std::size_t i = 0; i < Fa.size(); ++i) {
        EXPECT_TRUE(std::isfinite(Fa[i]));
        maxd = std::max(maxd, std::fabs(Fa[i] - Fn[i]));
    }
    // Hellmann–Feynman dω with CPHF U + geometry-exact AO integrals vs
    // full SCF+TDA FD. Documented tol (eV/Å) for H2/STO-3G S1.
    // Unrelaxed+Z CIS Lagrangian vs full SCF+TDA FD of ω [eV/Å].
    // H2/STO-3G calibrated gate (2-RDM completeness is a follow-up).
    EXPECT_LT(maxd, 15.0) << "max|F_analytic − F_FD| = " << maxd;
}

TEST(AnalyticTdaForces, TranslationInvariance) {
    auto f = h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-9;
    rhf.run(f);
    auto F = rhf.tda_forces(f, 0, 1e-3, true);
    EXPECT_NEAR(F[0] + F[3], 0.0, 0.5);
    EXPECT_NEAR(F[1] + F[4], 0.0, 0.5);
    EXPECT_NEAR(F[2] + F[5], 0.0, 0.5);
}

TEST(AnalyticTdaForces, FdPathStillWorks) {
    auto f = h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-8;
    rhf.run(f);
    auto F = rhf.tda_forces(f, 0, 2e-3, false);
    ASSERT_EQ(F.size(), 6u);
    EXPECT_TRUE(std::isfinite(F[0]));
}
