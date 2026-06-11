/// @file tests/cpu/driver_test_helpers.h
/// @brief Shared component stubs + frame factory for the CPU driver
///        validation / harness tests.
///
/// Phase 1 of the driver-taxonomy work introduced four
/// role-validation tests that each defined the same handful of stub
/// components (Stub{Integrator,ForceEvaluator,EnergyEvaluator,
/// Optimizer,BareInfra} + a 1-atom `make_trivial_frame`). Phase 5
/// added three async-infrastructure harness tests that each defined
/// the same `NoOpForce` component to satisfy `DynamicsDriver`'s
/// "≥ 1 component with Capability::Force" rule. This header
/// consolidates all of those stubs.
///
/// Test-specific scenario stubs (`AlwaysConverged`, `StopAtStep`,
/// `ConvergeAtStep`, `PositionCapture`, `ThreadRecorder`) stay
/// per-test — they encode the test's specific assertion path and
/// don't generalize.
///
/// All stubs have inline member function definitions so multiple
/// translation units including this header don't violate ODR.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace atv_test {

// ─── Generic component stubs ────────────────────────────────────────

/// @brief No-op force evaluator. Satisfies `DynamicsDriver`'s
/// "≥ 1 component with Capability::Force" requirement so that
/// infrastructure tests can use `cpu::md::MD` purely as a tick driver.
/// Default-binds to `md::Phase::force` so callers can use plain
/// `sim.add(NoOpForce{})`.
class NoOpForce : public ATV::cpu::ComponentBase {
  public:
    std::string name() const override { return "noop_force"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Evaluator;
    }
    ATV::Capability capabilities() const override {
        return ATV::Capability::Force;
    }
    std::vector<ATV::DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::force), 0}};
    }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

/// @brief Integrator role stub with optional step counter for test
/// assertions. Bound explicitly via `on(phase).add()` — no default
/// bindings.
class StubIntegrator : public ATV::cpu::ComponentBase {
  public:
    int step_count = 0;
    std::string name() const override { return "stub_integrator"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Integrator;
    }
    void launch(const ATV::cpu::LaunchContext&) override { ++step_count; }
};

/// @brief Evaluator stub advertising `Energy | Force` — satisfies both
/// `DynamicsDriver` and `OptimizeDriver` validation. No-op launch.
class StubForceEvaluator : public ATV::cpu::ComponentBase {
  public:
    std::string name() const override { return "stub_force"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Evaluator;
    }
    ATV::Capability capabilities() const override {
        return ATV::Capability::Energy | ATV::Capability::Force;
    }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

/// @brief Evaluator stub advertising only `Energy` — accepted by
/// `EvaluateDriver` and `OptimizeDriver`, rejected by `DynamicsDriver`
/// (no `Force`). Used to test the capability-axis validation contract.
class StubEnergyEvaluator : public ATV::cpu::ComponentBase {
  public:
    std::string name() const override { return "stub_energy"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Evaluator;
    }
    ATV::Capability capabilities() const override {
        return ATV::Capability::Energy;
    }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

/// @brief Optimizer role stub with optional step counter.
class StubOptimizer : public ATV::cpu::ComponentBase {
  public:
    int step_count = 0;
    std::string name() const override { return "stub_optimizer"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Optimizer;
    }
    void launch(const ATV::cpu::LaunchContext&) override { ++step_count; }
};

/// @brief `role()` defaults to `None`. Used to confirm structural-only
/// components don't satisfy any required role/capability slot.
class StubBareInfra : public ATV::cpu::ComponentBase {
  public:
    std::string name() const override { return "stub_bare"; }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

// ─── Frame factory ──────────────────────────────────────────────────

/// @brief Single-atom unit-mass frame in a 10 Å cubic box. Sized to
/// keep the SCF / MD machinery happy without doing any meaningful
/// physics — used only to drive the validation / iteration-status
/// path through `build()` + `run()`.
inline ATV::cpu::Frame make_trivial_frame() {
    ATV::cpu::Frame f(1);
    f.atoms.x[0] = 0.0;
    f.atoms.y[0] = 0.0;
    f.atoms.z[0] = 0.0;
    f.atoms.mass[0] = 1.0;
    f.atoms.type[0] = 0;
    f.box = ATV::Mat3::cubic(10.0);
    return f;
}

} // namespace atv_test
