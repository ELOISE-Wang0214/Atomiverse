/// @file tests/cpu/test_install_infra_hook.cpp
/// @brief Phase 1 (driver-infra-auto-injection) — assert the foundation:
///        (a) `ComponentBase::reads()` defaults to `Capability::None`;
///        (b) `ComponentBase::needs_neighbor_list()` defaults to `false`;
///        (c) a CPU category-driver subclass can override
///            `install_infra()` to inject a component, and the hook fires
///            *before* `validate_components()` so an injection can satisfy
///            an otherwise-unsatisfied validation rule.
///
/// Acceptance: T1 + T2 of `.claude/specs/driver-infra-auto-injection.md`.

#include <gtest/gtest.h>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

#include "cpu/driver_test_helpers.h"

using atv_test::make_trivial_frame;
using atv_test::StubForceEvaluator;
using atv_test::StubIntegrator;

namespace {

// ─── Foundation virtuals: defaults ──────────────────────────────────

TEST(InstallInfraHook, ReadsDefaultsToNone) {
    StubIntegrator c;
    EXPECT_EQ(c.reads(), ATV::Capability::None);
}

TEST(InstallInfraHook, NeedsNeighborListDefaultsToFalse) {
    StubIntegrator c;
    EXPECT_FALSE(c.needs_neighbor_list());
}

// ─── install_infra() hook fires before validate_components() ────────

/// @brief Counts how often `install_infra()` is called per `build()`.
class HookCounterDynamics : public ATV::cpu::DynamicsDriver {
  public:
    int install_calls = 0;
    StubIntegrator integ;

  protected:
    void install_infra() override {
        ++install_calls;
        // Inject the integrator from inside the hook — driver's
        // `validate_components()` must see this component on the
        // pending list, otherwise validation throws "missing Integrator".
        this->on(ATV::md::Phase::initial_integrate).add(integ);
    }
};

TEST(InstallInfraHook, HookFiresBeforeValidateComponents) {
    HookCounterDynamics d;
    // Only register the force evaluator from outside; the integrator
    // arrives via install_infra(). If validate_components() ran first,
    // the build would throw "missing Integrator".
    StubForceEvaluator force;
    d.on(ATV::md::Phase::force).add(force);

    EXPECT_NO_THROW(d.build());
    EXPECT_EQ(d.install_calls, 1);
}

TEST(InstallInfraHook, HookFiresExactlyOncePerBuild) {
    HookCounterDynamics d;
    StubForceEvaluator force;
    d.on(ATV::md::Phase::force).add(force);

    d.build();
    EXPECT_EQ(d.install_calls, 1);
}

// ─── Default install_infra is no-op (existing tests still pass) ─────

TEST(InstallInfraHook, DefaultInstallInfraIsNoOp) {
    // The plain `cpu::DynamicsDriver` (no subclass override) must still
    // throw on missing Integrator — proves the base `install_infra()` does
    // not silently inject anything.
    ATV::cpu::DynamicsDriver d;
    StubForceEvaluator force;
    d.on(ATV::md::Phase::force).add(force);
    EXPECT_THROW(d.build(), std::runtime_error);
}

} // namespace
