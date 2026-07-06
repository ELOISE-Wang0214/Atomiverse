/// @file tests/cpu/test_mc_scaffold.cpp
/// @brief Reference-free unit tests for the Monte Carlo sampling scaffold
///        (sub-spec mc-sampling-01-scaffold): sampling::Phase shape, the
///        ComponentRole::Mover tag, the pure-data mc::Transaction, the two
///        new cpu::Context handles, and SamplingDriver::validate_components.
///
/// No physics and no accept/reject dynamics are exercised here — there is no
/// Mover yet (that lands in mc-sampling-02-move). These tests only pin the
/// dormant machinery so the fifth (Sampling) driver category is green.

#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/driver/sampling_driver.h"
#include "atomiverse/cpu/mc/transaction.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;

// ── sampling::Phase execution skeleton ────────────────────────────────────

TEST(MCScaffold, SamplingPhaseShape) {
    EXPECT_EQ(static_cast<int>(sampling::Phase::propose), 0);
    EXPECT_EQ(static_cast<int>(sampling::Phase::energy), 1);
    EXPECT_EQ(static_cast<int>(sampling::Phase::end_of_step), 2);
    EXPECT_EQ(static_cast<int>(sampling::Phase::COUNT), 3);

    const auto names = sampling::phase_names(sampling::Phase{});
    ASSERT_EQ(names.size(), 3u);
    EXPECT_EQ(names[0], "propose");
    EXPECT_EQ(names[1], "energy");
    EXPECT_EQ(names[2], "end_of_step");
}

// ── ComponentRole::Mover structural tag ───────────────────────────────────

TEST(MCScaffold, MoverRoleIsDistinct) {
    EXPECT_NE(ComponentRole::Mover, ComponentRole::None);
    EXPECT_NE(ComponentRole::Mover, ComponentRole::Barostat);
    EXPECT_NE(ComponentRole::Mover, ComponentRole::Integrator);
    EXPECT_NE(ComponentRole::Mover, ComponentRole::Evaluator);
}

// ── mc::Transaction: pure-data round-trip + validation ─────────────────────

TEST(MCScaffold, TransactionRoundTrip) {
    std::vector<LocalInt> moved{2, 5};
    std::vector<StorageReal> ox{1, 2}, oy{3, 4}, oz{5, 6};
    std::vector<StorageReal> nx{7, 8}, ny{9, 10}, nz{11, 12};
    cpu::mc::Transaction t{moved, ox, oy, oz, nx, ny, nz};

    ASSERT_EQ(t.size(), 2u);
    EXPECT_EQ(t.moved[0], LocalInt(2));
    EXPECT_EQ(t.moved[1], LocalInt(5));
    EXPECT_EQ(t.old_x[0], StorageReal(1));
    EXPECT_EQ(t.old_z[1], StorageReal(6));
    EXPECT_EQ(t.new_x[0], StorageReal(7));
    EXPECT_EQ(t.new_z[1], StorageReal(12));
}

TEST(MCScaffold, TransactionRejectsMismatchedLengths) {
    std::vector<LocalInt> moved{0};
    std::vector<StorageReal> ok{1};
    std::vector<StorageReal> bad{1, 2}; // one element too many
    EXPECT_THROW((cpu::mc::Transaction{moved, ok, ok, ok, bad, ok, ok}),
                 std::invalid_argument);
}

TEST(MCScaffold, TransactionClearEmpties) {
    cpu::mc::Transaction t;
    t.moved = {1};
    t.old_x = {1};
    t.old_y = {2};
    t.old_z = {3};
    t.new_x = {4};
    t.new_y = {5};
    t.new_z = {6};
    ASSERT_EQ(t.size(), 1u);
    t.clear();
    EXPECT_EQ(t.size(), 0u);
    EXPECT_TRUE(t.new_z.empty());
}

// ── cpu::Context: the two new handles default to null ──────────────────────

TEST(MCScaffold, ContextHandlesDefaultNull) {
    cpu::Context ctx;
    EXPECT_EQ(ctx.proposal, nullptr);
    EXPECT_EQ(ctx.rng, nullptr);
}

// ── SamplingDriver::validate_components gate ───────────────────────────────

namespace {

/// @brief Minimal Mover-role stub (no proposal logic) for validation tests.
struct StubMover : cpu::ComponentBase {
    std::string name() const override { return "stub_mover"; }
    ComponentRole role() const override { return ComponentRole::Mover; }
    void launch(const cpu::LaunchContext&) override {}
};

/// @brief Minimal Energy-writing evaluator stub for validation tests.
struct StubEnergy : cpu::ComponentBase {
    std::string name() const override { return "stub_energy"; }
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    Capability capabilities() const override { return Capability::Energy; }
    void launch(const cpu::LaunchContext&) override {}
};

/// @brief Exposes the protected validate_components() for direct testing.
class TestSampling : public cpu::SamplingDriver {
  public:
    TestSampling() : cpu::SamplingDriver(/*temperature_K=*/300.0, /*seed=*/12345u) {}
    using cpu::SamplingDriver::validate_components;
};

} // namespace

TEST(MCScaffold, ValidateThrowsWithoutMover) {
    TestSampling d;
    d.on(sampling::Phase::energy).add(StubEnergy{});
    EXPECT_THROW(d.validate_components(), std::runtime_error);
}

TEST(MCScaffold, ValidateThrowsWithoutEnergyWriter) {
    TestSampling d;
    d.on(sampling::Phase::propose).add(StubMover{});
    EXPECT_THROW(d.validate_components(), std::runtime_error);
}

TEST(MCScaffold, ValidatePassesWithMoverAndEnergy) {
    TestSampling d;
    d.on(sampling::Phase::propose).add(StubMover{});
    d.on(sampling::Phase::energy).add(StubEnergy{});
    EXPECT_NO_THROW(d.validate_components());
}

TEST(MCScaffold, SamplingDriverRejectsNonPositiveTemperature) {
    EXPECT_THROW(
        { struct D : cpu::SamplingDriver { D() : cpu::SamplingDriver(0.0, 1u) {} } d; (void)d; },
        std::invalid_argument);
}
