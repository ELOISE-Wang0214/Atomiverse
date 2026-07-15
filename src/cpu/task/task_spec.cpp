/// @file cpu/task/task_spec.cpp
/// @brief TaskSpec validating constructors (the task-kind invariants).

#include "atomiverse/cpu/task/task_spec.h"

#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace ATV {
namespace cpu {
namespace task {

TaskSpec::TaskSpec(std::vector<ComponentEntry> entries, std::int64_t n_steps)
    : kind_(TaskKind::MD), entries_(std::move(entries)), n_steps_(n_steps) {
    if (n_steps_ < 1)
        throw std::invalid_argument(
            "task_spec: MD task needs n_steps >= 1, got " +
            std::to_string(n_steps_));
    // Checked narrowing (the write_back tag precedent): MD::run takes an
    // int budget, so reject here instead of truncating silently at
    // materialization.
    if (n_steps_ > std::numeric_limits<int>::max())
        throw std::invalid_argument(
            "task_spec: MD n_steps " + std::to_string(n_steps_) +
            " exceeds the int step-budget limit");
    for (const auto& e : entries_)
        if (std::holds_alternative<AM1EvaluatorParams>(e.params))
            throw std::invalid_argument(
                std::string("task_spec: MD task cannot carry the "
                            "evaluator-only component ") +
                component_name(e.params.index()));
}

TaskSpec::TaskSpec(AM1EvaluatorParams params)
    : kind_(TaskKind::SinglePoint),
      entries_{ComponentEntry{ComponentParams{std::move(params)},
                              std::nullopt}},
      n_steps_(1) {}

} // namespace task
} // namespace cpu
} // namespace ATV
