/// @file tests/cpu/test_driver_binding_kind.cpp
/// @brief arch-refactor-01-bindings — `DefaultBinding` gains a `PhaseKind`
///        discriminant and `Driver::add()` rejects cross-workflow default
///        bindings at registration time.
///
/// Covers:
///   - `cpu::md::MD::add(AM1Evaluator{})` throws (the evaluator loses its
///     hardcoded `default_bindings()` and falls onto the existing
///     "no default bindings" rejection path);
///   - kind-gate errors name BOTH workflows: Evaluate/Minimize rejecting an
///     MD-kind `VelocityVerlet` (rvalue `add()` branch), MD rejecting an
///     Evaluate-kind probe component (lvalue `add()` branch);
///   - same-kind `add()` is unchanged (MD + VelocityVerlet,
///     Minimize + LBFGS);
///   - plan-shape regression anchors for all five driver categories
///     (Evaluate / MD / Minimize / Hessian / MC): the compiled Plan's
///     per-phase component-name sequences must not change.
///
/// Tests use the FINAL API (`ATV::PhaseKind`, `ATV::make_binding`) directly.
/// On the pre-refactor tree this file fails to COMPILE — that is the
/// intended RED state.
///
/// Spec: `.claude/specs/arch-refactor-01-bindings.md` (Files E).

#include <gtest/gtest.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/mc/mc.h"
#include "atomiverse/cpu/mc/translate_move.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/minimizer/minimize.h"
#include "atomiverse/cpu/semiempirical/am1_evaluator.h"
#include "atomiverse/cpu/vib/hessian.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/plan.h"
#include "atomiverse/facade/stage.h"

#include "cpu/driver_test_helpers.h"

namespace {

bool message_mentions(const std::runtime_error& e, const std::string& needle) {
    return std::string(e.what()).find(needle) != std::string::npos;
}

/// @brief Component-name sequence of one compiled Plan phase slot.
std::vector<std::string>
slot_names(const ATV::Plan<ATV::cpu::ComponentBase>& plan, std::size_t phase) {
    std::vector<std::string> names;
    for (const auto& op : plan.slots[phase].ops)
        names.push_back(op.component->name());
    return names;
}

/// @brief Local Evaluate-kind probe: declares a default binding on
/// `evaluate::Phase::evaluate` through the final `make_binding` API. Used
/// to drive the kind gate on the LVALUE `add()` branch of an MD driver.
class EvaluateKindProbe : public ATV::cpu::ComponentBase {
  public:
    std::string name() const override { return "evaluate_kind_probe"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::Evaluator;
    }
    ATV::Capability capabilities() const override {
        return ATV::Capability::Energy;
    }
    std::vector<ATV::DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::evaluate::Phase::evaluate, 0)};
    }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

// ─── Violating evaluator: AM1Evaluator loses default_bindings() ────────
//
// Regression target: pre-refactor, `md.add(AM1Evaluator{})` silently landed
// the evaluator on `md::Phase::post_integrate` (the erased int `1` of
// `evaluate::Phase::evaluate`). Post-refactor the override is deleted, so
// `add()` takes the existing "has no default bindings. Use
// .on(phase).add()." rejection path.

TEST(DriverBindingKind, MdAddAm1EvaluatorThrows) {
    ATV::cpu::md::MD md;
    EXPECT_THROW(md.add(ATV::cpu::semiempirical::AM1Evaluator{}),
                 std::runtime_error);
}

// ─── Kind gate: error names BOTH workflows ──────────────────────────────

TEST(DriverBindingKind, EvaluateAddVelocityVerletThrowsWithBothKinds) {
    ATV::cpu::Evaluate ev;
    try {
        // rvalue `add()` branch.
        ev.add(ATV::cpu::md::VelocityVerlet{0.01});
        FAIL() << "expected throw — MD-kind default bindings registered "
                  "into an Evaluate driver";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "MD")) << e.what();
        EXPECT_TRUE(message_mentions(e, "Evaluate")) << e.what();
    }
}

TEST(DriverBindingKind, MinimizeAddVelocityVerletThrowsWithBothKinds) {
    ATV::cpu::minimizer::Minimize opt;
    try {
        // rvalue `add()` branch.
        opt.add(ATV::cpu::md::VelocityVerlet{0.01});
        FAIL() << "expected throw — MD-kind default bindings registered "
                  "into a Minimize driver";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "MD")) << e.what();
        EXPECT_TRUE(message_mentions(e, "Minimizer")) << e.what();
    }
}

TEST(DriverBindingKind, MdAddEvaluateKindProbeThrows) {
    ATV::cpu::md::MD md;
    EvaluateKindProbe probe;
    try {
        // lvalue `add()` branch (caller retains ownership).
        md.add(probe);
        FAIL() << "expected throw — Evaluate-kind default bindings "
                  "registered into an MD driver";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "Evaluate")) << e.what();
        EXPECT_TRUE(message_mentions(e, "MD")) << e.what();
    }
}

// ─── Same-kind add() unchanged ──────────────────────────────────────────

TEST(DriverBindingKind, MdAddSameKindVelocityVerletNoThrow) {
    ATV::cpu::md::MD md;
    EXPECT_NO_THROW(md.add(ATV::cpu::md::VelocityVerlet{0.01}));
}

TEST(DriverBindingKind, MinimizeAddSameKindLbfgsNoThrow) {
    ATV::cpu::minimizer::Minimize opt;
    EXPECT_NO_THROW(opt.add(ATV::cpu::minimizer::LBFGS{5, 0.01, 0.05}));
}

// ─── Plan-shape regression anchors (five driver categories) ─────────────
//
// The sequences below were MEASURED on the pre-refactor tree (dev,
// 2026-07-10) by compiling each minimal driver and walking
// `driver.plan()`. The kind discriminant is consumed inside `add()` and
// erased before `PlanBuilder`; the compiled Plan must be byte-identical.
// None of the minimal setups registers a pair evaluator, so no halo
// component is auto-injected under either an MPI or a serial build —
// the sequences are build-config invariant.

TEST(DriverBindingKind, EvaluatePlanShapeUnchanged) {
    ATV::cpu::Evaluate ev;
    atv_test::StubEnergyEvaluator eval;
    ev.on(ATV::evaluate::Phase::evaluate).add(eval);
    ev.build();

    const auto& plan = ev.plan();
    ASSERT_EQ(plan.slots.size(),
              static_cast<std::size_t>(ATV::evaluate::Phase::COUNT));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(
                                   ATV::evaluate::Phase::prepare)),
              (std::vector<std::string>{"zero_forces", "zero_energy"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(
                                   ATV::evaluate::Phase::evaluate)),
              (std::vector<std::string>{"stub_energy"}));
    EXPECT_TRUE(slot_names(plan, static_cast<std::size_t>(
                                     ATV::evaluate::Phase::finalize))
                    .empty());
}

TEST(DriverBindingKind, MdPlanShapeUnchanged) {
    ATV::cpu::md::MD sim;
    sim.set_integrator(ATV::cpu::md::VelocityVerlet{0.01});
    atv_test::NoOpForce force;
    sim.add(force);
    sim.build();

    const auto& plan = sim.plan();
    ASSERT_EQ(plan.slots.size(),
              static_cast<std::size_t>(ATV::md::Phase::COUNT));
    using P = ATV::md::Phase;
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::initial_integrate)),
              (std::vector<std::string>{"velocity_verlet"}));
    EXPECT_TRUE(
        slot_names(plan, static_cast<std::size_t>(P::post_integrate)).empty());
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::pre_exchange)),
              (std::vector<std::string>{"migrate"}));
    EXPECT_TRUE(
        slot_names(plan, static_cast<std::size_t>(P::pre_neighbor)).empty());
    EXPECT_TRUE(
        slot_names(plan, static_cast<std::size_t>(P::neighbor)).empty());
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::pre_force)),
              (std::vector<std::string>{"zero_forces", "zero_energy",
                                        "zero_stress"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::force)),
              (std::vector<std::string>{"noop_force"}));
    EXPECT_TRUE(
        slot_names(plan, static_cast<std::size_t>(P::post_force)).empty());
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::final_integrate)),
              (std::vector<std::string>{"velocity_verlet"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::end_of_step)),
              (std::vector<std::string>{"zero_kinetic", "kinetic_energy",
                                        "reduce", "publish"}));
}

TEST(DriverBindingKind, MinimizePlanShapeUnchanged) {
    ATV::cpu::minimizer::Minimize opt;
    opt.add(ATV::cpu::minimizer::LBFGS{5, 0.01, 0.05});
    atv_test::StubForceEvaluator force;
    opt.on(ATV::minimizer::Phase::force).add(force);
    opt.build();

    const auto& plan = opt.plan();
    ASSERT_EQ(plan.slots.size(),
              static_cast<std::size_t>(ATV::minimizer::Phase::COUNT));
    using P = ATV::minimizer::Phase;
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::search_direction)),
              (std::vector<std::string>{"lbfgs"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::line_search)),
              (std::vector<std::string>{"lbfgs"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::force)),
              (std::vector<std::string>{"zero_forces", "zero_energy", "migrate",
                                        "stub_force", "reduce"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::convergence_check)),
              (std::vector<std::string>{"lbfgs"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::end_of_step)),
              (std::vector<std::string>{"publish"}));
}

TEST(DriverBindingKind, HessianPlanShapeUnchanged) {
    ATV::cpu::vib::Hessian hess;
    atv_test::StubForceEvaluator force;
    hess.on(ATV::vib::Phase::force).add(force);
    hess.build();

    const auto& plan = hess.plan();
    ASSERT_EQ(plan.slots.size(),
              static_cast<std::size_t>(ATV::vib::Phase::COUNT));
    using P = ATV::vib::Phase;
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::setup)),
              (std::vector<std::string>{"vib_reference_geometry"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::displace)),
              (std::vector<std::string>{"vib_displace_step"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::force)),
              (std::vector<std::string>{"stub_force"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::accumulate)),
              (std::vector<std::string>{"vib_hessian_accumulator"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::diagonalize)),
              (std::vector<std::string>{"vib_symmetrize", "vib_mass_weighting",
                                        "hessian_diagonalizer"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::thermochem)),
              (std::vector<std::string>{"vib_thermochem"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::end_of_step)),
              (std::vector<std::string>{"vib_frequency_report",
                                        "vib_normal_mode_writer"}));
}

TEST(DriverBindingKind, McPlanShapeUnchanged) {
    ATV::cpu::mc::MC mc{300.0, 12345u};
    mc.set_mover(ATV::cpu::mc::TranslateMove{0.2});
    atv_test::StubEnergyEvaluator energy;
    mc.on(ATV::sampling::Phase::energy).add(energy);
    mc.build();

    const auto& plan = mc.plan();
    ASSERT_EQ(plan.slots.size(),
              static_cast<std::size_t>(ATV::sampling::Phase::COUNT));
    using P = ATV::sampling::Phase;
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::propose)),
              (std::vector<std::string>{"translate"}));
    EXPECT_EQ(slot_names(plan, static_cast<std::size_t>(P::energy)),
              (std::vector<std::string>{"zero_energy", "stub_energy"}));
    EXPECT_TRUE(
        slot_names(plan, static_cast<std::size_t>(P::end_of_step)).empty());
}

} // namespace
