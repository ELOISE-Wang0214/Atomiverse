#include <gtest/gtest.h>
#include <stdexcept>

#include "atomiverse/facade/plan.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;

// Dummy component for testing

class DummyComponent : public ComponentBase {
public:
    DummyComponent(std::string n) : name_(std::move(n)) {}

    std::string name() const override { return name_; }
    void launch(const LaunchContext&) override {}

private:
    std::string name_;
};

TEST(Plan, LayerSeparation) {
    // Plan builds on host with no CUDA API calls.
    DummyComponent comp("test");
    PlanBuilder builder;
    builder.add({0, 0, &comp});
    auto plan = builder.build(static_cast<int>(md::Phase::COUNT));

    EXPECT_EQ(plan.slots.size(), static_cast<int>(md::Phase::COUNT));
    EXPECT_EQ(plan.slots[0].ops.size(), 1u);
    EXPECT_EQ(plan.slots[0].ops[0].component->name(), "test");
}

TEST(Plan, OrderBasedScheduling) {
    DummyComponent b("B");
    DummyComponent a("A");

    PlanBuilder builder;
    int force = static_cast<int>(md::Phase::force);
    builder.add({force, 10, &b});
    builder.add({force, 1, &a});

    auto plan = builder.build(static_cast<int>(md::Phase::COUNT));
    auto& slot = plan.slots[force];
    ASSERT_EQ(slot.ops.size(), 2u);
    EXPECT_EQ(slot.ops[0].component->name(), "A");
    EXPECT_EQ(slot.ops[1].component->name(), "B");
}

TEST(Plan, InvalidPhaseThrows) {
    DummyComponent comp("bad");
    PlanBuilder builder;
    builder.add({999, 0, &comp});
    EXPECT_THROW(builder.build(static_cast<int>(md::Phase::COUNT)),
                 std::runtime_error);
}
