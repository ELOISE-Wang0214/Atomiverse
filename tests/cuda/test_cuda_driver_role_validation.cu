/// @file tests/cuda/test_cuda_driver_role_validation.cu
/// @brief Assert that the three CUDA category bases — cuda::EvaluateDriver
///        / DynamicsDriver / OptimizeDriver — throw at build() time when
///        their required role + capability set is not satisfied, and
///        build successfully when it is.
///
/// Mirrors tests/cpu/test_cpu_driver_role_validation.cpp.
/// Acceptance: ac-008 in
/// `.claude/specs/cuda-driver-taxonomy-01-foundations.acceptance.md`.
///
/// `build()` only validates the component set and compiles the Plan — it
/// performs no device work — so these tests need no CUDA device.

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/cuda/dynamics_driver.cuh"
#include "atomiverse/cuda/evaluate_driver.cuh"
#include "atomiverse/cuda/optimize_driver.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace {

bool message_mentions(const std::runtime_error& e, const std::string& needle) {
    return std::string(e.what()).find(needle) != std::string::npos;
}

// ─── CUDA component stubs ───────────────────────────────────────────

/// @brief Integrator-role stub. No default bindings — bound explicitly.
class StubIntegrator : public ATV::cuda::ComponentBase {
  public:
    std::string name() const override { return "stub_integrator"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Integrator;
    }
    void launch(const ATV::cuda::LaunchContext&) override {}
};

/// @brief Evaluator-role stub advertising `Energy | Force`.
class StubForceEvaluator : public ATV::cuda::ComponentBase {
  public:
    std::string name() const override { return "stub_force"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Evaluator;
    }
    ATV::Capability capabilities() const override {
        return ATV::Capability::Energy | ATV::Capability::Force;
    }
    void launch(const ATV::cuda::LaunchContext&) override {}
};

/// @brief Evaluator-role stub advertising only `Energy` (no `Force`).
class StubEnergyEvaluator : public ATV::cuda::ComponentBase {
  public:
    std::string name() const override { return "stub_energy"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Evaluator;
    }
    ATV::Capability capabilities() const override {
        return ATV::Capability::Energy;
    }
    void launch(const ATV::cuda::LaunchContext&) override {}
};

/// @brief Optimizer-role stub.
class StubOptimizer : public ATV::cuda::ComponentBase {
  public:
    std::string name() const override { return "stub_optimizer"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Optimizer;
    }
    void launch(const ATV::cuda::LaunchContext&) override {}
};

/// @brief `role()` defaults to `None` — structural-only stub.
class StubBareInfra : public ATV::cuda::ComponentBase {
  public:
    std::string name() const override { return "stub_bare"; }
    void launch(const ATV::cuda::LaunchContext&) override {}
};

} // namespace

// ─── DynamicsDriver ─────────────────────────────────────────────────

TEST(CUDADriverRoleValidation, DynamicsThrowsWithoutIntegrator) {
    using namespace ATV;
    cuda::DynamicsDriver d;
    StubForceEvaluator force;
    d.on(md::Phase::force).add(force);
    try {
        d.build();
        FAIL() << "expected throw — DynamicsDriver missing Integrator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Integrator")) << e.what();
    }
}

TEST(CUDADriverRoleValidation, DynamicsThrowsWithoutForceEvaluator) {
    using namespace ATV;
    cuda::DynamicsDriver d;
    StubIntegrator integ;
    StubEnergyEvaluator energy_only;
    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(energy_only);
    try {
        d.build();
        FAIL() << "expected throw — DynamicsDriver missing Force-capable "
                  "evaluator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Force")) << e.what();
    }
}

TEST(CUDADriverRoleValidation, DynamicsBuildsWithIntegratorAndForce) {
    using namespace ATV;
    cuda::DynamicsDriver d;
    StubIntegrator integ;
    StubForceEvaluator force;
    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(force);
    EXPECT_NO_THROW(d.build());
}

// ─── OptimizeDriver ─────────────────────────────────────────────────

TEST(CUDADriverRoleValidation, OptimizeThrowsWithoutOptimizer) {
    using namespace ATV;
    cuda::OptimizeDriver o;
    StubForceEvaluator force;
    o.on(minimizer::Phase::force).add(force);
    try {
        o.build();
        FAIL() << "expected throw — OptimizeDriver missing Optimizer";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Optimizer")) << e.what();
    }
}

TEST(CUDADriverRoleValidation, OptimizeThrowsWithoutEnergyOrForceEvaluator) {
    using namespace ATV;
    cuda::OptimizeDriver o;
    StubOptimizer opt;
    StubBareInfra infra;
    o.on(minimizer::Phase::search_direction).add(opt);
    o.on(minimizer::Phase::end_of_step).add(infra);
    try {
        o.build();
        FAIL()
            << "expected throw — OptimizeDriver missing Energy/Force evaluator";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Energy")
                    || message_mentions(e, "Force"))
            << e.what();
    }
}

TEST(CUDADriverRoleValidation, OptimizeBuildsWithOptimizerAndEnergy) {
    using namespace ATV;
    cuda::OptimizeDriver o;
    StubOptimizer opt;
    StubEnergyEvaluator energy;
    o.on(minimizer::Phase::search_direction).add(opt);
    o.on(minimizer::Phase::force).add(energy);
    EXPECT_NO_THROW(o.build());
}

// ─── EvaluateDriver ─────────────────────────────────────────────────

TEST(CUDADriverRoleValidation, EvaluateThrowsWithoutEvaluator) {
    using namespace ATV;
    cuda::EvaluateDriver e;
    StubBareInfra infra;
    e.on(evaluate::Phase::prepare).add(infra);
    try {
        e.build();
        FAIL() << "expected throw — EvaluateDriver missing Evaluator";
    } catch (const std::runtime_error& ex) {
        EXPECT_TRUE(message_mentions(ex, "Evaluator")) << ex.what();
    }
}

TEST(CUDADriverRoleValidation, EvaluateBuildsWithEnergyEvaluator) {
    using namespace ATV;
    cuda::EvaluateDriver e;
    StubEnergyEvaluator energy;
    e.on(evaluate::Phase::evaluate).add(energy);
    EXPECT_NO_THROW(e.build());
}
