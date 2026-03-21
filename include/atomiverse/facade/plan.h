/// @file plan.h
/// @brief Plan + PlanBuilder -- execution plan compilation. CUDA-free.
///
/// The Plan is the compiled schedule that bridges Facade and Backend.
/// PlanBuilder resolves component bindings into sorted StageSlots.
/// Built once, executed many times.
///
/// Templated on Component type to support both CPU and CUDA backends.
///
/// @see Driver::build() for how Plans are created
/// @see Backend::step() for how Plans are executed
#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/types.h"

namespace ATV {

/// @brief One resolved operation within a StageSlot.
/// @tparam Component Component base type (cpu::ComponentBase or ATV::ComponentBase).
template <typename Component>
struct Operation {
    Component* component; ///< Non-owning pointer to the component.
    Scheduling scheduling; ///< Cached scheduling trait (blocking or deferrable).
    int frequency;         ///< Cached execution frequency (0 = every step).
};

/// @brief All operations for one phase, sorted by order.
/// @tparam Component Component base type.
template <typename Component>
struct StageSlot {
    int phase;                              ///< Phase index within the driver's Phase enum.
    std::vector<Operation<Component>> ops;  ///< Operations sorted by registration order.
};

/// @brief Compiled execution plan -- the schedule between Facade and Backend.
///
/// Contains one StageSlot per phase. Built once by PlanBuilder, then
/// executed many times by Backend::step().
///
/// @tparam Component Component base type.
template <typename Component>
struct Plan {
    std::vector<StageSlot<Component>> slots; ///< One StageSlot per phase.
};

/// @brief Compiles component bindings into a Plan.
///
/// Collects binding entries (phase + order + component), then sorts
/// operations within each phase by order value (stable sort).
///
/// @tparam Component Component base type.
/// @see Driver::build() which uses PlanBuilder internally
template <typename Component>
class PlanBuilder {
public:
    /// @brief Binding entry: associates a component with a phase and order.
    struct Entry {
        int phase;              ///< Target phase index.
        int order;              ///< Ordering within phase (lower = earlier).
        Component* component;   ///< Non-owning pointer to the component.
    };

    /// @brief Register a binding entry.
    /// @param entry Binding to register.
    void add(Entry entry) { entries_.push_back(entry); }

    /// @brief Compile all registered bindings into a Plan.
    /// @param phase_count Number of phases (typically PhaseEnum::COUNT).
    /// @return Compiled Plan with sorted operations per phase.
    /// @throws std::runtime_error If a component binds to an invalid phase.
    Plan<Component> build(int phase_count) {
        Plan<Component> plan;
        plan.slots.resize(phase_count);
        for (int p = 0; p < phase_count; ++p)
            plan.slots[p].phase = p;

        struct SE { int order; Operation<Component> op; };
        std::vector<std::vector<SE>> slot_entries(phase_count);

        for (auto& e : entries_) {
            if (e.phase < 0 || e.phase >= phase_count)
                throw std::runtime_error("Component '" + e.component->name() +
                    "' binds to invalid phase " + std::to_string(e.phase));
            Operation<Component> op{e.component, e.component->scheduling(),
                         e.component->frequency()};
            slot_entries[e.phase].push_back({e.order, op});
        }

        for (int p = 0; p < phase_count; ++p) {
            std::stable_sort(slot_entries[p].begin(), slot_entries[p].end(),
                [](const SE& a, const SE& b) { return a.order < b.order; });
            for (auto& se : slot_entries[p])
                plan.slots[p].ops.push_back(se.op);
        }

        return plan;
    }

private:
    std::vector<Entry> entries_;  ///< Accumulated binding entries.
};

} // namespace ATV
