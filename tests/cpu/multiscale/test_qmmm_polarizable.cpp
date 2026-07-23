/// @file tests/cpu/multiscale/test_qmmm_polarizable.cpp
/// @brief Unit tests for induced-dipole SCC and QmmmPolarizable.
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/qmmm_polarizable.h"

using namespace ATV;

TEST(InducedDipole, SingleSiteUniformField) {
    // μ = α E for one site (no mutual terms).
    const double alpha = 1.5; // Å³
    const double e_perm[3] = {0.2, -0.1, 0.05};
    const double pos[3] = {0, 0, 0};
    std::vector<double> mu;
    const double e_pol =
        cpu::solve_induced_dipoles(1, pos, &alpha, e_perm, /*thole=*/0.0,
                                   /*tol=*/1e-12, /*max_iter=*/10, mu);
    ASSERT_EQ(mu.size(), 3u);
    EXPECT_NEAR(mu[0], alpha * e_perm[0], 1e-14);
    EXPECT_NEAR(mu[1], alpha * e_perm[1], 1e-14);
    EXPECT_NEAR(mu[2], alpha * e_perm[2], 1e-14);
    const double e_expect =
        units::coulomb * (-0.5 * (mu[0] * e_perm[0] + mu[1] * e_perm[1] +
                                  mu[2] * e_perm[2]));
    EXPECT_NEAR(e_pol, e_expect, 1e-12);
    EXPECT_LT(e_pol, 0.0);
}

TEST(InducedDipole, SccThrowsWhenUnconverged) {
    // Two close polarizable sites with huge α: force max_iter=1 then throw.
    const double pos[6] = {0, 0, 0, 0.5, 0, 0};
    const double alpha[2] = {50.0, 50.0};
    const double e_perm[6] = {1, 0, 0, -1, 0, 0};
    std::vector<double> mu;
    EXPECT_THROW(cpu::solve_induced_dipoles(2, pos, alpha, e_perm, 0.39, 1e-16,
                                            /*max_iter=*/1, mu),
                 std::runtime_error);
}

TEST(QmmmPolarizable, PolarizationLowersEnergy) {
    // H2 QM + one MM site with charge and α; E(α>0) < E(α=0).
    auto energy = [](double alpha) {
        cpu::Frame f(3);
        f.span(kAtom_x)[0] = 0;
        f.span(kAtom_y)[0] = 0;
        f.span(kAtom_z)[0] = 0;
        f.span(kAtom_x)[1] = 0.74;
        f.span(kAtom_y)[1] = 0;
        f.span(kAtom_z)[1] = 0;
        f.span(kAtom_x)[2] = 3.0;
        f.span(kAtom_y)[2] = 0;
        f.span(kAtom_z)[2] = 0;
        f.span(kAtom_type)[0] = 1;
        f.span(kAtom_type)[1] = 1;
        f.span(kAtom_type)[2] = 1;
        f.span(kAtom_tag)[0] = 1;
        f.span(kAtom_tag)[1] = 2;
        f.span(kAtom_tag)[2] = 3;
        f.simbox = cpu::Box::cubic(20);
        f.meta.get(kEnergy) = 0;
        cpu::Context ctx;
        cpu::QmmmPolarizable qm({1, 2}, /*mm_q=*/0.5, alpha, /*thole=*/0.39,
                                0.01, 3.4, 10.0);
        qm.build({ctx, f});
        qm.launch({ctx, f, 0, 0});
        return std::make_pair(f.meta.get(kEnergy), qm.last_pol_energy());
    };
    const auto [e0, p0] = energy(0.0);
    const auto [e1, p1] = energy(2.0);
    EXPECT_TRUE(std::isfinite(e0) && std::isfinite(e1));
    EXPECT_NEAR(p0, 0.0, 1e-14);
    EXPECT_LT(p1, 0.0);
    EXPECT_LT(e1, e0);
}
