/// @file test_pw_evaluator.cpp
/// @brief PWEvaluator — Evaluate-driver integration + PWSolver parity.
///
/// Verifies the plane-wave SCF is reachable through the canonical
/// `cpu::Evaluate + Evaluator-component` shape (taxonomy rule: SCF is an
/// inner solver behind an Evaluator, never an outer driver), and that
/// the component path reproduces the `pw_scf()` shim bit-for-bit.
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/evaluate.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/pw/pw_evaluator.h"
#include "atomiverse/cpu/pw/pw_scf.h"
#include "atomiverse/facade/stage.h"
#include <cmath>
#include <gtest/gtest.h>
#include <type_traits>

using namespace ATV;
using namespace ATV::cpu::pw;

// Taxonomy guards: PWSolver is an inner solver (plain class), PWEvaluator
// is the Plan-facing component.
static_assert(!std::is_base_of_v<cpu::ComponentBase, PWSolver>,
              "PWSolver must not be a Plan component");
static_assert(std::is_base_of_v<cpu::ComponentBase, PWEvaluator>,
              "PWEvaluator must be a cpu::ComponentBase");

namespace {
cpu::Frame make_si_fcc(double a) {
    cpu::Frame f(2);
    f.atoms.type[0] = 14;
    f.atoms.type[1] = 14;
    f.atoms.mass[0] = 28.0855;
    f.atoms.mass[1] = 28.0855;
    f.atoms.x[0] = 0;
    f.atoms.y[0] = 0;
    f.atoms.z[0] = 0;
    f.atoms.x[1] = a / 4.0;
    f.atoms.y[1] = a / 4.0;
    f.atoms.z[1] = a / 4.0;
    double h = a / 2.0;
    Mat3 cell;
    cell(0, 0) = 0;
    cell(0, 1) = h;
    cell(0, 2) = h;
    cell(1, 0) = h;
    cell(1, 1) = 0;
    cell(1, 2) = h;
    cell(2, 0) = h;
    cell(2, 1) = h;
    cell(2, 2) = 0;
    f.box = cell;
    return f;
}

PWConfig small_config() {
    PWConfig cfg;
    cfg.ecut_ry = 10.0;
    cfg.k_mesh = 2;
    cfg.tol = 1e-6;
    cfg.max_iter = 20;
    return cfg;
}
} // namespace

TEST(PWEvaluator, RoleAndCapabilities) {
    PWEvaluator pwe{small_config()};
    EXPECT_EQ(pwe.role(), ComponentRole::Evaluator);
    EXPECT_TRUE(has_capability(pwe.capabilities(), Capability::Energy));
    EXPECT_FALSE(has_capability(pwe.capabilities(), Capability::Force));
}

TEST(PWEvaluator, EvaluateDriverWritesFrameEnergyInEv) {
    auto frame = make_si_fcc(5.43);

    PWEvaluator pwe{small_config()};
    cpu::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate).add(pwe);
    eval.build();
    eval.run(frame);

    const auto& r = pwe.last_result();
    EXPECT_TRUE(std::isfinite(r.energy));
    EXPECT_NEAR(r.energy_ha, r.energy * kRy2Hartree, 1e-12);
    // frame.energy follows the project convention (eV).
    EXPECT_NEAR(frame.scalars.energy, r.energy * kRy2eV, 1e-9);
}

TEST(PWEvaluator, MatchesPwScfShim) {
    auto cfg = small_config();

    auto frame_shim = make_si_fcc(5.43);
    auto r_shim = pw_scf(frame_shim, cfg.upf_path, cfg.ecut_ry, cfg.k_mesh,
                         cfg.tol, cfg.max_iter);

    auto frame_eval = make_si_fcc(5.43);
    PWEvaluator pwe{cfg};
    cpu::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate).add(pwe);
    eval.build();
    eval.run(frame_eval);

    EXPECT_TRUE(std::isfinite(r_shim.energy));
    EXPECT_NEAR(pwe.last_result().energy, r_shim.energy, 1e-10)
        << "Evaluate + PWEvaluator must reproduce the pw_scf() shim";
    EXPECT_EQ(pwe.last_result().iterations, r_shim.iterations);
}
