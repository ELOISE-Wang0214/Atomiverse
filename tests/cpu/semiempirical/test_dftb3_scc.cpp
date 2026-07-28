/// @file tests/cpu/semiempirical/test_dftb3_scc.cpp
/// Unit gates on the DFTB3 second-/third-order SCC kernel: the potential must
/// be the exact derivative of the energy it is paired with.
#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "atomiverse/cpu/semiempirical/internal/dftb3_scc.h"

using ATV::cpu::semiempirical::internal::Dftb3Scc;

namespace {

// O/H/H mio-style Hubbard data at the water geometry used by the DFTB3
// regression (O at origin, H at (0.96,0,0), H at (-0.24,0.93,0)) [Å].
Dftb3Scc water_kernel() {
    const double r_oh1 = 0.96;
    const double r_oh2 = std::sqrt(0.24 * 0.24 + 0.93 * 0.93);
    const double r_hh = std::sqrt(1.2 * 1.2 + 0.93 * 0.93);
    return Dftb3Scc({13.4, 10.591, 10.591}, {-3.81, -1.36, -1.36},
                    {0.0, r_oh1, r_oh2, r_oh1, 0.0, r_hh, r_oh2, r_hh, 0.0});
}

/// Polar, charge-neutral test fluctuation [e].
const std::vector<double> kDq = {-0.6, 0.3, 0.3};

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        s += a[i] * b[i];
    return s;
}

} // namespace

TEST(Dftb3Scc, OnsiteKernelIsHubbardData) {
    const Dftb3Scc scc = water_kernel();
    ASSERT_EQ(scc.size(), 3);
    // γ_AA = U_A and Γ_AA = ∂γ_AA/∂q_A = U^d_A — the convention E³ is built on.
    EXPECT_DOUBLE_EQ(scc.gamma(0, 0), 13.4);
    EXPECT_DOUBLE_EQ(scc.gamma(1, 1), 10.591);
    EXPECT_DOUBLE_EQ(scc.gamma_charge_derivative(0, 0), -3.81);
    EXPECT_DOUBLE_EQ(scc.gamma_charge_derivative(2, 2), -1.36);
}

TEST(Dftb3Scc, GammaChargeDerivativeIsNotSymmetric) {
    const Dftb3Scc scc = water_kernel();
    // Γ_AB ≠ Γ_BA is exactly why ∂E³/∂Δq_C needs two distinct contractions.
    EXPECT_GT(std::fabs(scc.gamma_charge_derivative(0, 1)
                        - scc.gamma_charge_derivative(1, 0)),
              1e-6);
    EXPECT_DOUBLE_EQ(scc.gamma(0, 1), scc.gamma(1, 0));
}

TEST(Dftb3Scc, ThirdOrderPotentialSatisfiesEulerIdentity) {
    const Dftb3Scc scc = water_kernel();
    // E³ is homogeneous of degree 3 in Δq ⇒ Σ_A Δq_A V³_A = 3 E³ exactly.
    const double lhs = dot(kDq, scc.third_order_potential(kDq));
    const double rhs = 3.0 * scc.third_order_energy(kDq);
    ASSERT_GT(std::fabs(rhs), 1e-6);
    EXPECT_LE(std::fabs(lhs - rhs) / std::fabs(rhs), 1e-13);
}

TEST(Dftb3Scc, SecondOrderPotentialSatisfiesEulerIdentity) {
    const Dftb3Scc scc = water_kernel();
    // E² is homogeneous of degree 2 in Δq ⇒ Σ_A Δq_A V²_A = 2 E² exactly.
    const double lhs = dot(kDq, scc.second_order_potential(kDq));
    const double rhs = 2.0 * scc.second_order_energy(kDq);
    ASSERT_GT(std::fabs(rhs), 1e-6);
    EXPECT_LE(std::fabs(lhs - rhs) / std::fabs(rhs), 1e-13);
}

TEST(Dftb3Scc, SecondOrderPotentialMatchesFiniteDifference) {
    const Dftb3Scc scc = water_kernel();
    const std::vector<double> v = scc.second_order_potential(kDq);
    constexpr double h = 1e-6; // e
    for (int a = 0; a < scc.size(); ++a) {
        std::vector<double> up = kDq, dn = kDq;
        up[static_cast<std::size_t>(a)] += h;
        dn[static_cast<std::size_t>(a)] -= h;
        const double num =
            (scc.second_order_energy(up) - scc.second_order_energy(dn))
            / (2.0 * h);
        EXPECT_NEAR(v[static_cast<std::size_t>(a)], num, 1e-8);
    }
}

TEST(Dftb3Scc, PotentialMatchesFiniteDifferenceOfEnergy) {
    const Dftb3Scc scc = water_kernel();
    const std::vector<double> v = scc.potential(kDq);
    constexpr double h = 1e-6; // e
    auto energy = [&](const std::vector<double>& q) {
        return scc.second_order_energy(q) + scc.third_order_energy(q);
    };
    for (int a = 0; a < scc.size(); ++a) {
        std::vector<double> up = kDq, dn = kDq;
        up[static_cast<std::size_t>(a)] += h;
        dn[static_cast<std::size_t>(a)] -= h;
        const double num = (energy(up) - energy(dn)) / (2.0 * h);
        EXPECT_NEAR(v[static_cast<std::size_t>(a)], num, 1e-8);
    }
}

TEST(Dftb3Scc, VanishingChargeFluctuationGivesZeroEnergyAndPotential) {
    const Dftb3Scc scc = water_kernel();
    const std::vector<double> zero(3, 0.0);
    // Homonuclear symmetry ⇒ Δq ≡ 0 ⇒ DFTB3 reduces exactly to DFTB2/DFTB0.
    EXPECT_DOUBLE_EQ(scc.second_order_energy(zero), 0.0);
    EXPECT_DOUBLE_EQ(scc.third_order_energy(zero), 0.0);
    for (const double vi : scc.potential(zero))
        EXPECT_DOUBLE_EQ(vi, 0.0);
}

TEST(Dftb3Scc, RejectsMismatchedInputSizes) {
    EXPECT_THROW(Dftb3Scc({10.0, 10.0}, {-1.0}, std::vector<double>(4, 0.0)),
                 std::invalid_argument);
    EXPECT_THROW(Dftb3Scc({10.0, 10.0}, {-1.0, -1.0}, {0.0, 1.0, 1.0}),
                 std::invalid_argument);
    EXPECT_THROW(Dftb3Scc({}, {}, {}), std::invalid_argument);
    const Dftb3Scc scc = water_kernel();
    EXPECT_THROW(scc.potential({0.1, -0.1}), std::invalid_argument);
}
