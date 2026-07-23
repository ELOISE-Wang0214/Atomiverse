/// @file tests/cpu/potentials/test_pair_born.cpp
/// @brief Unit checks for cpu::PairBorn energy/force vs analytic shifted form.
#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/potentials/pair_born.h"

using namespace ATV;

namespace {

cpu::Frame dimer(double r) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = r;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1.0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 0;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.simbox = cpu::Box::cubic(50.0);
    return f;
}

double U_un(double A, double rho, double sigma, double C, double D, double r) {
    return A * std::exp((sigma - r) / rho) - C / std::pow(r, 6) +
           D / std::pow(r, 8);
}

double dUdr(double A, double rho, double sigma, double C, double D, double r) {
    return -A / rho * std::exp((sigma - r) / rho) + 6.0 * C / std::pow(r, 7) -
           8.0 * D / std::pow(r, 9);
}

} // namespace

TEST(PairBorn, AnalyticEnergyAndForce) {
    const double A = 1000.0, rho = 3.0, sigma = 1.0, C = 50.0, D = 20.0;
    const double rc = 8.0, r = 2.5;
    const double U = U_un(A, rho, sigma, C, D, r) - U_un(A, rho, sigma, C, D, rc);
    const double du = dUdr(A, rho, sigma, C, D, r);

    auto f = dimer(r);
    f.meta.get(kEnergy) = 0;
    for (int i = 0; i < 2; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;

    cpu::Context ctx;
    cpu::PairBorn born(A, rho, sigma, C, D, rc);
    born.build({ctx, f});
    born.launch({ctx, f, 0, 0});

    EXPECT_NEAR(f.meta.get(kEnergy), U, 1e-12);
    EXPECT_NEAR(f.span(kAtom_fx)[0], du, 1e-12);
    EXPECT_NEAR(f.span(kAtom_fx)[1], -du, 1e-12);
    EXPECT_NEAR(f.span(kAtom_fy)[0], 0.0, 1e-15);
    EXPECT_NEAR(f.span(kAtom_fz)[0], 0.0, 1e-15);
    born.destroy({ctx, f});
}

TEST(PairBorn, ZeroAtCutoff) {
    const double A = 1000.0, rho = 3.0, sigma = 1.0, C = 50.0, D = 20.0;
    const double rc = 8.0;
    const double r = rc - 1e-6;
    const double U = U_un(A, rho, sigma, C, D, r) - U_un(A, rho, sigma, C, D, rc);

    auto f = dimer(r);
    f.meta.get(kEnergy) = 0;
    for (int i = 0; i < 2; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;

    cpu::Context ctx;
    cpu::PairBorn born(A, rho, sigma, C, D, rc);
    born.build({ctx, f});
    born.launch({ctx, f, 0, 0});
    EXPECT_NEAR(f.meta.get(kEnergy), U, 1e-12);
    EXPECT_LT(std::fabs(f.meta.get(kEnergy)), 1e-3);
    born.destroy({ctx, f});
}

TEST(PairBorn, RejectsInvalidParams) {
    EXPECT_THROW(cpu::PairBorn(1.0, 0.0, 1.0, 1.0, 0.0, 5.0),
                 std::invalid_argument);
    EXPECT_THROW(cpu::PairBorn(1.0, 0.3, 1.0, 1.0, 0.0, 0.0),
                 std::invalid_argument);
    EXPECT_THROW(cpu::PairBorn(-1.0, 0.3, 1.0, 1.0, 0.0, 5.0),
                 std::invalid_argument);
    EXPECT_THROW(cpu::PairBorn(1.0, 0.3, 1.0, -1.0, 0.0, 5.0),
                 std::invalid_argument);
}

TEST(PairBorn, NameIsBorn) {
    cpu::PairBorn born(1.0, 0.3, 1.0, 1.0, 0.0, 5.0);
    EXPECT_EQ(born.name(), "born");
    EXPECT_EQ(born.as_neighbor_consumer(), nullptr);
}
