/// @file tests/cpu/gpw/test_cpu_gpw_dftu.cpp
/// @brief Unit tests for Dudarev DFT+U helpers and Ni/O PADE table lookup.
///        No full NiO SCF (that is the regression).

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/dft_plus_u.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/molrs.hpp"

using ATV::cpu::gpw::GaussianBasis;
using ATV::cpu::gpw::GthPotential;
using ATV::cpu::gpw::PlusU;
using ATV::cpu::gpw::internal::dudarev_energy_potential;
using ATV::cpu::gpw::internal::mulliken_occupation;
using ATV::cpu::linalg::RealMatrix;

TEST(GpwDftU, NickelPadeTableLookup) {
    GthPotential pot(molrs::Element::Ni, "GTH-PADE-q18");
    EXPECT_EQ(pot.zion, 18);
    EXPECT_DOUBLE_EQ(pot.rloc, 0.35);
    ASSERT_EQ(pot.channels.size(), 3u);
    EXPECT_EQ(pot.channels[2].l, 2); // d projector
    EXPECT_EQ(pot.channels[2].nproj, 1);

    GaussianBasis bas(molrs::Element::Ni, "SZV-GTH-PADE-q18");
    ASSERT_EQ(bas.sets.size(), 1u);
    EXPECT_EQ(bas.sets[0].lmin, 0);
    EXPECT_EQ(bas.sets[0].lmax, 2);
    EXPECT_EQ(bas.sets[0].nshell_l(2), 1); // one d shell
}

TEST(GpwDftU, OxygenPadeTableLookup) {
    GthPotential pot(molrs::Element::O, "GTH-PADE-q6");
    EXPECT_EQ(pot.zion, 6);
    GaussianBasis bas(molrs::Element::O, "SZV-GTH-PADE");
    ASSERT_EQ(bas.sets.size(), 1u);
    EXPECT_EQ(bas.sets[0].lmax, 1);
}

TEST(GpwDftU, DudarevEmptyShellZero) {
    RealMatrix occ = xt::zeros<double>({0, 0});
    RealMatrix V;
    const double e = dudarev_energy_potential(0.2, occ, V);
    EXPECT_DOUBLE_EQ(e, 0.0);
    EXPECT_EQ(V.shape(0), 0u);
}

TEST(GpwDftU, DudarevHalfFilledSingleOrbital) {
    // n = 0.5 on a 1×1 shell → E = (U/2)(0.5 - 0.25) = U/8; V = U(0.5 - 0.5)=0
    RealMatrix occ = xt::zeros<double>({1, 1});
    occ(0, 0) = 0.5;
    RealMatrix V;
    const double U = 0.2; // Ha
    const double e = dudarev_energy_potential(U, occ, V);
    EXPECT_NEAR(e, U / 8.0, 1e-14);
    ASSERT_EQ(V.shape(0), 1u);
    EXPECT_NEAR(V(0, 0), 0.0, 1e-14);
}

TEST(GpwDftU, DudarevEmptyAndFull) {
    // empty: Tr=0 → E=0, V=U/2 on diagonal
    {
        RealMatrix occ = xt::zeros<double>({2, 2});
        RealMatrix V;
        const double e = dudarev_energy_potential(0.4, occ, V);
        EXPECT_NEAR(e, 0.0, 1e-14);
        EXPECT_NEAR(V(0, 0), 0.2, 1e-14);
        EXPECT_NEAR(V(1, 1), 0.2, 1e-14);
    }
    // full identity: Tr=2, Tr(n²)=2 → E=0; V = U(0.5 - 1) = -U/2
    {
        RealMatrix occ = xt::zeros<double>({2, 2});
        occ(0, 0) = 1.0;
        occ(1, 1) = 1.0;
        RealMatrix V;
        const double e = dudarev_energy_potential(0.4, occ, V);
        EXPECT_NEAR(e, 0.0, 1e-14);
        EXPECT_NEAR(V(0, 0), -0.2, 1e-14);
    }
}

TEST(GpwDftU, PlusUPublicMembers) {
    PlusU u(molrs::Element::Ni, 2, 6.0);
    EXPECT_EQ(static_cast<int>(u.element), 28);
    EXPECT_EQ(u.l, 2);
    EXPECT_DOUBLE_EQ(u.u_minus_j_ev, 6.0);
    // eV → Ha conversion used by the solver
    const double U_ha = u.u_minus_j_ev / ATV::units::kHartreeInEv;
    EXPECT_NEAR(U_ha * ATV::units::kHartreeInEv, 6.0, 1e-12);
}

TEST(GpwDftU, MullikenOccupationIdentityS) {
    // P = diag(1,0), S = I → n = diag(1,0) on full space as subspace
    RealMatrix P = xt::zeros<double>({2, 2});
    P(0, 0) = 1.0;
    RealMatrix S = xt::zeros<double>({2, 2});
    S(0, 0) = 1.0;
    S(1, 1) = 1.0;
    const std::vector<int> idx = {0, 1};
    RealMatrix n = mulliken_occupation(P, S, idx);
    EXPECT_NEAR(n(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(n(1, 1), 0.0, 1e-14);
    EXPECT_NEAR(n(0, 1), 0.0, 1e-14);
}
