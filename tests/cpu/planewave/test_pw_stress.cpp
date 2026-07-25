/// @file tests/cpu/pw/test_pw_stress.cpp
/// @brief Unit gates for plane-wave cell stress (Nielsen–Martin + FD oracle).

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/planewave/pw_evaluator.h"
#include "atomiverse/cpu/planewave/pw_solver.h"
#include "atomiverse/cpu/scf/xc_functional.h"

using namespace ATV;
using namespace ATV::cpu::pw;

namespace {

cpu::Frame make_h_cubic(double L) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.simbox = cpu::Box::cubic(L);
    return f;
}

PW cheap_cfg() {
    PW cfg;
    cfg.ecut_ry = 5.0;
    cfg.k_mesh = 1;
    cfg.tol = 1e-4;
    cfg.max_iter = 8;
    cfg.mix_beta = 0.5;
    cfg.xc = cpu::scf::XCType::LSDA; // LDA XC stress term is complete
    return cfg;
}

} // namespace

TEST(PWStress, CubicMonoatomicOffDiagonalSoft) {
    auto f = make_h_cubic(6.0);
    const auto s = pw_stress(f, cheap_cfg());
    for (int k = 0; k < 6; ++k)
        EXPECT_TRUE(std::isfinite(s[static_cast<std::size_t>(k)])) << k;
    EXPECT_LT(std::abs(s[3]), 0.5);
    EXPECT_LT(std::abs(s[4]), 0.5);
    EXPECT_LT(std::abs(s[5]), 0.5);
    EXPECT_NEAR(s[0], s[1], 0.5);
    EXPECT_NEAR(s[1], s[2], 0.5);
}

TEST(PWStress, AnalyticMatchesFdCheapH) {
    // Dual-gate A4: max |σ_an − σ_FD| ≤ 0.15 eV/Å³ on cubic H (toy PP).
    auto f = make_h_cubic(6.0);
    auto cfg = cheap_cfg();
    const auto san = pw_stress(f, cfg);
    const auto sfd = pw_stress_fd(f, cfg, 2e-3);
    double max_abs = 0.0;
    for (int v = 0; v < 6; ++v) {
        const double d = std::abs(san[static_cast<std::size_t>(v)]
                                  - sfd[static_cast<std::size_t>(v)]);
        max_abs = std::max(max_abs, d);
        EXPECT_TRUE(std::isfinite(san[static_cast<std::size_t>(v)]));
        EXPECT_TRUE(std::isfinite(sfd[static_cast<std::size_t>(v)]));
    }
    EXPECT_LE(max_abs, 0.15)
        << "an=" << san[0] << "," << san[1] << "," << san[2] << " fd=" << sfd[0]
        << "," << sfd[1] << "," << sfd[2] << " max|Δ|=" << max_abs;
}

TEST(PWStress, FdRejectsBadStep) {
    auto f = make_h_cubic(6.0);
    EXPECT_THROW(pw_stress_fd(f, cheap_cfg(), 0.0), std::invalid_argument);
    EXPECT_THROW(pw_stress_fd(f, cheap_cfg(), -1e-3), std::invalid_argument);
}

TEST(PWStress, EvaluatorWithStressWritesKStress) {
    auto frame = make_h_cubic(6.0);
    PWEvaluator pwe{cheap_cfg()};
    pwe.with_stress(true);
    EXPECT_TRUE(pwe.compute_stress());
    EXPECT_FALSE(pwe.compute_forces());

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add_ref(pwe);
    eval.enable_timing(false);
    eval.build();
    eval.run(frame);

    EXPECT_TRUE(std::isfinite(frame.meta.get(kEnergy)));
    const auto& st = frame.meta.get(kStress);
    for (int k = 0; k < 6; ++k)
        EXPECT_TRUE(std::isfinite(st[static_cast<std::size_t>(k)]));
}
