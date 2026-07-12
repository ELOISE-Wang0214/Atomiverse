/// @file tests/cpu/test_cpu_driver_role_validation.cpp
/// @brief Phase 1.5 — assert that the Evaluate / MD / Minimize presets throw at
///        build() time when their required role + capability set is not
///        satisfied, and build successfully when it is.
///
/// The category bases these tests used to instantiate (EvaluateDriver /
/// DynamicsDriver / OptimizeDriver) no longer exist — `validate_components()`
/// moved down into the presets, which is where the rule was always enforced.
///
/// Acceptance: ac-003 in `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#include <gtest/gtest.h>

#include <string>

#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/minimizer/minimize.h"
#include "atomiverse/facade/stage.h"

#include "cpu/driver_test_helpers.h"

using atv_test::StubBareInfra;
using atv_test::StubEnergyEvaluator;
using atv_test::StubForceEvaluator;
using atv_test::StubIntegrator;
using atv_test::StubOptimizer;

namespace {

bool message_mentions(const std::runtime_error& e, const std::string& needle) {
    return std::string(e.what()).find(needle) != std::string::npos;
}

} // namespace

// ─── Dynamics (cpu::md::MD) ──────────────────────────────────────────

TEST(CPUDriverRoleValidation, DynamicsThrowsWithoutIntegrator) {
    using namespace ATV;
    cpu::md::MD d;
    StubForceEvaluator force;
    d.on(md::Phase::force).add(force);
    try {
        d.build();
        FAIL() << "expected throw — MD missing Integrator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Integrator")) << e.what();
    }
}

TEST(CPUDriverRoleValidation, DynamicsThrowsWithoutForceEvaluator) {
    using namespace ATV;
    cpu::md::MD d;
    StubIntegrator integ;
    StubEnergyEvaluator energy_only;
    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(energy_only);
    try {
        d.build();
        FAIL() << "expected throw — MD missing Force-capable "
                  "evaluator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Force")) << e.what();
    }
}

TEST(CPUDriverRoleValidation, DynamicsBuildsWithIntegratorAndForce) {
    using namespace ATV;
    cpu::md::MD d;
    StubIntegrator integ;
    StubForceEvaluator force;
    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(force);
    EXPECT_NO_THROW(d.build());
}

// ─── Optimize (cpu::minimizer::Minimize) ─────────────────────────────

TEST(CPUDriverRoleValidation, OptimizeThrowsWithoutOptimizer) {
    using namespace ATV;
    cpu::minimizer::Minimize o;
    StubForceEvaluator force;
    o.on(minimizer::Phase::force).add(force);
    try {
        o.build();
        FAIL() << "expected throw — Minimize missing Optimizer";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Optimizer")) << e.what();
    }
}

TEST(CPUDriverRoleValidation, OptimizeThrowsWithoutEnergyOrForceEvaluator) {
    using namespace ATV;
    cpu::minimizer::Minimize o;
    StubOptimizer opt;
    StubBareInfra infra;
    o.on(minimizer::Phase::search_direction).add(opt);
    o.on(minimizer::Phase::end_of_step).add(infra);
    try {
        o.build();
        FAIL()
            << "expected throw — Minimize missing Energy/Force evaluator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Energy")
                    || message_mentions(e, "Force"))
            << e.what();
    }
}

TEST(CPUDriverRoleValidation, OptimizeBuildsWithOptimizerAndEnergy) {
    using namespace ATV;
    cpu::minimizer::Minimize o;
    StubOptimizer opt;
    StubEnergyEvaluator energy;
    o.on(minimizer::Phase::search_direction).add(opt);
    o.on(minimizer::Phase::force).add(energy);
    EXPECT_NO_THROW(o.build());
}

// ─── Evaluate (cpu::Evaluate) ────────────────────────────────────────

TEST(CPUDriverRoleValidation, EvaluateThrowsWithoutEvaluator) {
    using namespace ATV;
    cpu::Evaluate e;
    StubBareInfra infra;
    e.on(evaluate::Phase::prepare).add(infra);
    try {
        e.build();
        FAIL() << "expected throw — Evaluate missing Evaluator";
    } catch (const std::runtime_error& ex) {
        EXPECT_TRUE(message_mentions(ex, "Evaluator")) << ex.what();
    }
}

TEST(CPUDriverRoleValidation, EvaluateBuildsWithEnergyEvaluator) {
    using namespace ATV;
    cpu::Evaluate e;
    StubEnergyEvaluator energy;
    e.on(evaluate::Phase::evaluate).add(energy);
    EXPECT_NO_THROW(e.build());
}
