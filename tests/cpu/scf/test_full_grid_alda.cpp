/// @file tests/cpu/scf/test_full_grid_alda.cpp
/// @brief Unit checks for full MO-product ALDA kernel (symmetry, diagonal limit).
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/internal/alda_kernel.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;
using namespace ATV::cpu::scf;

namespace {

cpu::Frame h2_frame(double r = 0.74) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = r;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.simbox = cpu::Box::cubic(20.0);
    return f;
}

} // namespace

TEST(FullGridAlda, KernelSymmetric) {
    auto f = h2_frame();
    RKS rks(XCType::LSDA);
    rks.basis = BasisSet::STO3G;
    rks.tol = 1e-9;
    rks.run(f);
    ASSERT_TRUE(rks.converged());

    IntegralPack pack{f, BasisSet::STO3G};
    const MolecularGrid grid = build_molecular_grid(
        pack.atoms.data(), static_cast<int>(pack.atoms.size()),
        rks.grid_level, rks.n_radial, rks.n_angular);
    const RealMatrix K = internal::alda_kernel_matrix(
        rks.determinant(), pack, grid, ExcitationChannel::singlet);
    ASSERT_EQ(K.shape()[0], K.shape()[1]);
    const std::size_t n = K.shape()[0];
    double max_asym = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            max_asym = std::max(max_asym, std::fabs(K(i, j) - K(j, i)));
    EXPECT_LT(max_asym, 1e-12);
}

TEST(FullGridAlda, DiracFxAnalytic) {
    // f_αα(ρ) = −(4/9) 2^{1/3} C_x ρ^{-2/3}
    constexpr double Cx = 0.7385587663820224;
    const double rho = 0.2;
    const double expect =
        -(4.0 / 9.0) * std::pow(2.0, 1.0 / 3.0) * Cx * std::pow(rho, -2.0 / 3.0);
    EXPECT_NEAR(internal::alda_fx_dirac(rho), expect, 1e-14);
}

TEST(FullGridAlda, FullGridDiffersFromDiagonalOnly) {
    // Full-grid ALDA must move the S1 root relative to pure CIS@RKS (no fxc).
    auto f = h2_frame();
    RKS lsda(XCType::LSDA);
    lsda.basis = BasisSet::STO3G;
    lsda.tol = 1e-9;
    lsda.run(f);
    ASSERT_TRUE(lsda.converged());
    auto e_full = lsda.tda(1);
    EXPECT_TRUE(e_full.alda);
    ASSERT_FALSE(e_full.energy_ev.empty());

    RKS pbe(XCType::PBE);
    pbe.basis = BasisSet::STO3G;
    pbe.tol = 1e-9;
    pbe.run(f);
    auto e_pbe = pbe.tda(1);
    EXPECT_FALSE(e_pbe.alda);
    EXPECT_GT(std::fabs(e_full.energy_ev[0] - e_pbe.energy_ev[0]), 1e-3);

    // H2/STO-3G: n_ov = 1; K(0,0) finite and nonzero.
    IntegralPack pack{f, BasisSet::STO3G};
    const MolecularGrid grid = build_molecular_grid(
        pack.atoms.data(), static_cast<int>(pack.atoms.size()),
        lsda.grid_level, lsda.n_radial, lsda.n_angular);
    const RealMatrix K = internal::alda_kernel_matrix(
        lsda.determinant(), pack, grid, ExcitationChannel::singlet);
    EXPECT_EQ(K.shape()[0], 1u);
    EXPECT_TRUE(std::isfinite(K(0, 0)));
    EXPECT_NE(K(0, 0), 0.0);
}

TEST(FullGridAlda, LsdaTdaPositiveFinite) {
    auto f = h2_frame();
    RKS rks(XCType::LSDA);
    rks.basis = BasisSet::STO3G;
    rks.tol = 1e-9;
    rks.run(f);
    ASSERT_TRUE(rks.converged());
    auto es = rks.tda(1);
    EXPECT_TRUE(es.alda);
    EXPECT_GT(es.energy_ev[0], 0.5);
    EXPECT_LT(es.energy_ev[0], 40.0);
}
