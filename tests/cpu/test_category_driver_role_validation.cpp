/// @file tests/cpu/test_category_driver_role_validation.cpp
/// @brief Phase 1.5 — assert that EvaluateDriver / DynamicsDriver /
///        OptimizeDriver throw at build() time when their required role +
///        capability set is not satisfied, and build successfully when it
///        is.
///
/// Acceptance: ac-003 in `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#include <gtest/gtest.h>

#include <string>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/evaluate_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/optimize_driver.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace {

// ─── Stubs ──────────────────────────────────────────────────────────

struct StubIntegrator : ATV::cpu::ComponentBase {
    std::string name() const override { return "stub_integrator"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Integrator; }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

struct StubForceEvaluator : ATV::cpu::ComponentBase {
    std::string name() const override { return "stub_force"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Evaluator; }
    ATV::Capability   capabilities() const override {
        return ATV::Capability::Energy | ATV::Capability::Force;
    }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

struct StubEnergyEvaluator : ATV::cpu::ComponentBase {
    std::string name() const override { return "stub_energy"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Evaluator; }
    ATV::Capability   capabilities() const override { return ATV::Capability::Energy; }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

struct StubOptimizer : ATV::cpu::ComponentBase {
    std::string name() const override { return "stub_optimizer"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Optimizer; }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

struct StubBareInfra : ATV::cpu::ComponentBase {
    std::string name() const override { return "stub_bare"; }
    // role() defaults to None. Should not satisfy any required slot.
    void launch(const ATV::cpu::LaunchContext&) override {}
};

bool message_mentions(const std::runtime_error& e, const std::string& needle) {
    return std::string(e.what()).find(needle) != std::string::npos;
}

} // namespace

// ─── DynamicsDriver ─────────────────────────────────────────────────

TEST(CategoryDriverRoleValidation, DynamicsThrowsWithoutIntegrator) {
    using namespace ATV;
    cpu::DynamicsDriver d;
    StubForceEvaluator force;
    d.on(md::Phase::force).add(force);
    try {
        d.build();
        FAIL() << "expected throw — DynamicsDriver missing Integrator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Integrator")) << e.what();
    }
}

TEST(CategoryDriverRoleValidation, DynamicsThrowsWithoutForceEvaluator) {
    using namespace ATV;
    cpu::DynamicsDriver d;
    StubIntegrator integ;
    StubEnergyEvaluator energy_only;
    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(energy_only);
    try {
        d.build();
        FAIL() << "expected throw — DynamicsDriver missing Force-capable evaluator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Force")) << e.what();
    }
}

TEST(CategoryDriverRoleValidation, DynamicsBuildsWithIntegratorAndForce) {
    using namespace ATV;
    cpu::DynamicsDriver d;
    StubIntegrator integ;
    StubForceEvaluator force;
    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(force);
    EXPECT_NO_THROW(d.build());
}

// ─── OptimizeDriver ─────────────────────────────────────────────────

TEST(CategoryDriverRoleValidation, OptimizeThrowsWithoutOptimizer) {
    using namespace ATV;
    cpu::OptimizeDriver o;
    StubForceEvaluator force;
    o.on(minimizer::Phase::force).add(force);
    try {
        o.build();
        FAIL() << "expected throw — OptimizeDriver missing Optimizer";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Optimizer")) << e.what();
    }
}

TEST(CategoryDriverRoleValidation, OptimizeThrowsWithoutEnergyOrForceEvaluator) {
    using namespace ATV;
    cpu::OptimizeDriver o;
    StubOptimizer opt;
    StubBareInfra infra;
    o.on(minimizer::Phase::search_direction).add(opt);
    o.on(minimizer::Phase::end_of_step).add(infra);
    try {
        o.build();
        FAIL() << "expected throw — OptimizeDriver missing Energy/Force evaluator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Energy") ||
                    message_mentions(e, "Force"))
            << e.what();
    }
}

TEST(CategoryDriverRoleValidation, OptimizeBuildsWithOptimizerAndEnergy) {
    using namespace ATV;
    cpu::OptimizeDriver o;
    StubOptimizer opt;
    StubEnergyEvaluator energy;
    o.on(minimizer::Phase::search_direction).add(opt);
    o.on(minimizer::Phase::force).add(energy);
    EXPECT_NO_THROW(o.build());
}

// ─── EvaluateDriver ─────────────────────────────────────────────────

TEST(CategoryDriverRoleValidation, EvaluateThrowsWithoutEvaluator) {
    using namespace ATV;
    cpu::EvaluateDriver e;
    StubBareInfra infra;
    e.on(evaluate::Phase::prepare).add(infra);
    try {
        e.build();
        FAIL() << "expected throw — EvaluateDriver missing Evaluator";
    } catch (const std::runtime_error& ex) {
        EXPECT_TRUE(message_mentions(ex, "Evaluator")) << ex.what();
    }
}

TEST(CategoryDriverRoleValidation, EvaluateBuildsWithEnergyEvaluator) {
    using namespace ATV;
    cpu::EvaluateDriver e;
    StubEnergyEvaluator energy;
    e.on(evaluate::Phase::evaluate).add(energy);
    EXPECT_NO_THROW(e.build());
}
