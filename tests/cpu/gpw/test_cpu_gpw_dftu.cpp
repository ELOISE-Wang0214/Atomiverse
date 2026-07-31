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

TEST(GpwDftU, MullikenOccupationKMatchesReal) {
    // Complex Γ with real P,S must recover the real Mulliken block.
    using Cplx = std::complex<double>;
    using ATV::cpu::linalg::CplxMatrix;
    using ATV::cpu::gpw::internal::accumulate_mulliken_occupation_k;
    RealMatrix P = xt::zeros<double>({2, 2});
    P(0, 0) = 0.7;
    P(1, 1) = 0.3;
    P(0, 1) = P(1, 0) = 0.1;
    RealMatrix S = xt::zeros<double>({2, 2});
    S(0, 0) = 1.0;
    S(1, 1) = 1.0;
    S(0, 1) = S(1, 0) = 0.2;
    const std::vector<int> idx = {0, 1};
    RealMatrix n_real = mulliken_occupation(P, S, idx);

    CplxMatrix Pk = xt::zeros<Cplx>({2, 2});
    CplxMatrix Sk = xt::zeros<Cplx>({2, 2});
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            Pk(i, j) = P(i, j);
            Sk(i, j) = S(i, j);
        }
    RealMatrix n_k = xt::zeros<double>({2, 2});
    accumulate_mulliken_occupation_k(n_k, Pk, Sk, idx, /*w=*/1.0);
    EXPECT_NEAR(n_k(0, 0), n_real(0, 0), 1e-12);
    EXPECT_NEAR(n_k(1, 1), n_real(1, 1), 1e-12);
    EXPECT_NEAR(n_k(0, 1), n_real(0, 1), 1e-12);
}

TEST(GpwDftU, SubtractDftuPulayHandChecked) {
    // nb = 2, one correlated AO (atom 0, l = 0, idx = {0}), U = 1 Ha.
    // n_00 = ½[(PS)+(SP)]_00 = p00 + p01·s, V̂ = diag(U(½ − n00), 0),
    // M = ½(V̂P + PV̂):  M00 = v·p00, M01 = M10 = ½·v·p01, M11 = 0.
    using ATV::cpu::gpw::internal::AoShellTag;
    using ATV::cpu::gpw::internal::PlusUHa;
    using ATV::cpu::gpw::internal::subtract_dftu_pulay;
    const double p00 = 0.7, p01 = 0.2, p11 = 0.4, s = 0.3;
    RealMatrix P = xt::zeros<double>({2, 2});
    P(0, 0) = p00;
    P(0, 1) = P(1, 0) = p01;
    P(1, 1) = p11;
    RealMatrix S = xt::zeros<double>({2, 2});
    S(0, 0) = S(1, 1) = 1.0;
    S(0, 1) = S(1, 0) = s;
    const std::vector<AoShellTag> tags = {{0, 0}, {1, 0}};
    const std::vector<int> z_of_atom = {28, 8};
    const std::vector<PlusUHa> slots = {{/*l=*/0, /*U_ha=*/1.0}};
    const std::vector<int> z_keys = {28};

    RealMatrix W = xt::zeros<double>({2, 2});
    subtract_dftu_pulay(W, P, S, /*channel_scale=*/1.0, tags, z_of_atom,
                        slots, z_keys);
    const double n00 = p00 + p01 * s;
    const double v = 1.0 * (0.5 - n00);
    EXPECT_NEAR(W(0, 0), -v * p00, 1e-14);
    EXPECT_NEAR(W(0, 1), -0.5 * v * p01, 1e-14);
    EXPECT_NEAR(W(1, 0), -0.5 * v * p01, 1e-14);
    EXPECT_NEAR(W(1, 1), 0.0, 1e-14);
}

TEST(GpwDftU, SubtractDftuPulayKGammaMatchesReal) {
    // One real k = Γ (weight 1) must reproduce the real fold exactly.
    using Cplx = std::complex<double>;
    using ATV::cpu::gpw::internal::AoShellTag;
    using ATV::cpu::gpw::internal::PlusUHa;
    using ATV::cpu::gpw::internal::subtract_dftu_pulay;
    using ATV::cpu::gpw::internal::subtract_dftu_pulay_k;
    using ATV::cpu::linalg::CplxMatrix;
    RealMatrix P = xt::zeros<double>({2, 2});
    P(0, 0) = 0.7;
    P(0, 1) = P(1, 0) = 0.2;
    P(1, 1) = 0.4;
    RealMatrix S = xt::zeros<double>({2, 2});
    S(0, 0) = S(1, 1) = 1.0;
    S(0, 1) = S(1, 0) = 0.3;
    const std::vector<AoShellTag> tags = {{0, 0}, {1, 0}};
    const std::vector<int> z_of_atom = {28, 8};
    const std::vector<PlusUHa> slots = {{/*l=*/0, /*U_ha=*/1.0}};
    const std::vector<int> z_keys = {28};

    RealMatrix W_real = xt::zeros<double>({2, 2});
    subtract_dftu_pulay(W_real, P, S, /*channel_scale=*/2.0, tags, z_of_atom,
                        slots, z_keys);

    std::vector<CplxMatrix> P_k(1), S_k(1),
        W_k(1, CplxMatrix(xt::zeros<Cplx>({2, 2})));
    P_k[0] = xt::zeros<Cplx>({2, 2});
    S_k[0] = xt::zeros<Cplx>({2, 2});
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            P_k[0](i, j) = P(i, j);
            S_k[0](i, j) = S(i, j);
        }
    subtract_dftu_pulay_k(W_k, P_k, S_k, {1.0}, /*channel_scale=*/2.0, tags,
                          z_of_atom, slots, z_keys);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            EXPECT_NEAR(W_k[0](i, j).real(), W_real(i, j), 1e-14);
            EXPECT_NEAR(W_k[0](i, j).imag(), 0.0, 1e-14);
        }
}

TEST(GpwDftU, BuildRejectsGammaRksPlusU) {
    // Active plus_u on the Γ closed-shell RKS path must throw (the Γ RKS SCF
    // has no +U fold — silent ignore is forbidden).
    using namespace ATV::cpu::gpw;
    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 8.0;
    cfg.n_occ = 1;
    cfg.nonperiodic = false;
    cfg.plus_u.push_back(PlusU{molrs::Element::H, /*l=*/0, /*u_ev=*/4.0});
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double r0[3] = {3.3, 4.0, 4.0};
    const double r1[3] = {4.7, 4.0, 4.0};
    std::vector<GpwAtom> atoms{GpwAtom(r0, basis, pot, 1),
                               GpwAtom(r1, basis, pot, 1)};
    EXPECT_THROW(cfg.build(atoms), std::invalid_argument);
}

TEST(GpwDftU, BuildAcceptsKmeshPlusU) {
    using namespace ATV::cpu::gpw;
    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 8.0;
    cfg.n_alpha = 2;
    cfg.n_beta = 0;
    cfg.k_mesh[0] = cfg.k_mesh[1] = cfg.k_mesh[2] = 2;
    cfg.plus_u.push_back(PlusU{molrs::Element::H, /*l=*/0, /*u_ev=*/4.0});
    cfg.local_pp_on_grid = true;
    cfg.nonperiodic = false;
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double r0[3] = {3.3, 4.0, 4.0};
    const double r1[3] = {4.7, 4.0, 4.0};
    std::vector<GpwAtom> atoms{GpwAtom(r0, basis, pot, 1),
                               GpwAtom(r1, basis, pot, 1)};
    EXPECT_NO_THROW(cfg.build(atoms));
}

TEST(GpwDftU, KpointGammaLimitMatchesRealGammaUks) {
    // Explicit k={Γ} complex path + PLUS_U vs real Γ-UKS + PLUS_U.
    using namespace ATV::cpu::gpw;
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double L = 8.0;
    const double r0[3] = {0.5 * L - 0.7, 0.5 * L, 0.5 * L};
    const double r1[3] = {0.5 * L + 0.7, 0.5 * L, 0.5 * L};
    std::vector<GpwAtom> atoms{GpwAtom(r0, basis, pot, 1),
                               GpwAtom(r1, basis, pot, 1)};

    auto base = [&]() {
        GpwDft cfg;
        cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = L;
        cfg.cutoff_ry = 100.0;
        cfg.n_alpha = 2;
        cfg.n_beta = 0;
        // Variational ½(V̂S+SV̂) scatter (2026-07-31): ~90+ DIIS iters here.
        cfg.max_iter = 300;
        cfg.tol = 1e-7;
        cfg.xc = XcKind::PADE;
        cfg.plus_u.push_back(PlusU{molrs::Element::H, /*l=*/0, /*u_ev=*/4.0});
        return cfg;
    };

    GpwDft g = base();
    g.build(atoms);
    const GpwResult rg = g.run(atoms);

    GpwDft k = base();
    k.k_points = {{0.0, 0.0, 0.0}};
    k.k_point_weights = {1.0};
    k.build(atoms);
    const GpwResult rk = k.run(atoms);

    ASSERT_TRUE(rg.converged);
    ASSERT_TRUE(rk.converged);
    EXPECT_NEAR(rk.e_total, rg.e_total, 1e-6);
    // Total energy is the decisive identity; e_dftu decomposition can
    // differ slightly under Hermitian projection / grid path.
    EXPECT_NEAR(rk.e_dftu, rg.e_dftu, 1e-4);
    EXPECT_GT(std::fabs(rg.e_dftu), 1e-8);
}
