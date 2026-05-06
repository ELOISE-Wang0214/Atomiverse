/// @file tests/cpu/test_evaluate_driver_capabilities.cpp
/// @brief Phase 1.6 — Evaluate accepts evaluators *by capability*, not by
///        ForceProvider role. An energy-only evaluator (no Force capability)
///        is valid in EvaluateDriver and rejected by DynamicsDriver.
///
/// Acceptance: ac-004 in `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#include <gtest/gtest.h>

#include <string>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/evaluate_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace {

struct EnergyOnlyEvaluator : ATV::cpu::ComponentBase {
    int launch_count = 0;
    std::string name() const override { return "energy_only"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Evaluator; }
    ATV::Capability   capabilities() const override { return ATV::Capability::Energy; }
    void launch(const ATV::cpu::LaunchContext&) override { ++launch_count; }
};

struct StubIntegrator : ATV::cpu::ComponentBase {
    std::string name() const override { return "stub_integrator"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Integrator; }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

ATV::cpu::Frame make_trivial_frame() {
    ATV::cpu::Frame f(1);
    f.atoms.x[0] = 0.0;
    f.atoms.y[0] = 0.0;
    f.atoms.z[0] = 0.0;
    f.atoms.mass[0] = 1.0;
    f.atoms.type[0] = 0;
    f.box = ATV::Mat3::cubic(10.0);
    return f;
}

} // namespace

TEST(EvaluateCapabilities, EnergyOnlyAcceptedByEvaluate) {
    using namespace ATV;
    cpu::EvaluateDriver e;
    EnergyOnlyEvaluator energy;
    e.on(evaluate::Phase::evaluate).add(energy);
    EXPECT_NO_THROW(e.build());

    auto frame = make_trivial_frame();
    e.run(frame);
    EXPECT_EQ(energy.launch_count, 1)
        << "Evaluate should run the evaluator exactly once";
}

TEST(EvaluateCapabilities, EnergyOnlyRejectedByDynamics) {
    using namespace ATV;
    cpu::DynamicsDriver d;
    StubIntegrator integ;
    EnergyOnlyEvaluator energy;
    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(energy);
    try {
        d.build();
        FAIL() << "Dynamics requires Force capability — energy-only must reject";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Force"), std::string::npos)
            << e.what();
    }
}
