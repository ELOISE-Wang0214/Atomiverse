/// @file tests/cpu/scf/test_scf_solver.cpp
/// @brief Concrete methods `RHF` / `UKS` converge on H2O and are stable
///        across independent instances and repeated `run()` calls.

#ifdef ATV_HAS_HF

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

namespace {

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = 0.7572;
    f.span(ATV::kAtom_x)[2] = -0.7572;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.0;
    f.span(ATV::kAtom_y)[2] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_z)[1] = -0.5860;
    f.span(ATV::kAtom_z)[2] = -0.5860;
    f.span(ATV::kAtom_type)[0] = 8;
    f.span(ATV::kAtom_type)[1] = 1;
    f.span(ATV::kAtom_type)[2] = 1;
    f.span(ATV::kAtom_mass)[0] = 15.999;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.span(ATV::kAtom_mass)[2] = 1.008;
    f.simbox = ATV::cpu::Box::cubic(50.0);
    return f;
}

} // namespace

TEST(RHF, H2OConvergesAcrossIndependentInstances) {
    auto frame = make_h2o();
    cpu::scf::RHF a;
    a.basis = cpu::scf::BasisSet::STO3G;
    const auto& result = a.run(frame);

    EXPECT_TRUE(result.converged) << "H2O / STO-3G RHF should converge";
    EXPECT_GT(result.iterations, 0);
    EXPECT_LT(result.iterations, a.max_iter);

    cpu::scf::RHF b;
    b.basis = cpu::scf::BasisSet::STO3G;
    const auto& second = b.run(frame);
    EXPECT_NEAR(result.energy, second.energy, 1e-10);
    ASSERT_EQ(result.P_alpha.shape(0), second.P_alpha.shape(0));
    for (std::size_t i = 0; i < result.P_alpha.shape(0); ++i) {
        for (std::size_t j = 0; j < result.P_alpha.shape(1); ++j) {
            EXPECT_NEAR(result.P_alpha(i, j), second.P_alpha(i, j), 1e-10);
        }
    }
}

TEST(RHF, RunIsIdempotentAcrossInvocations) {
    auto frame = make_h2o();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;

    auto& r1 = rhf.run(frame);
    const double e1 = r1.energy;
    auto P1 = r1.P_alpha;

    auto& r2 = rhf.run(frame);
    EXPECT_NEAR(e1, r2.energy, 1e-10);
    for (std::size_t i = 0; i < P1.shape(0); ++i)
        for (std::size_t j = 0; j < P1.shape(1); ++j)
            EXPECT_NEAR(P1(i, j), r2.P_alpha(i, j), 1e-10);
}

TEST(UKS, H2OConverges) {
    auto frame = make_h2o();
    // Closed-shell water: 5 alpha + 5 beta
    cpu::scf::UKS uks(/*n_alpha=*/5, /*n_beta=*/5, cpu::scf::XCType::LSDA);
    uks.basis = cpu::scf::BasisSet::STO3G;
    uks.max_iter = 80;
    const auto& result = uks.run(frame);
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.energy, 0.0);
}

#endif // ATV_HAS_HF
