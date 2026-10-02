/// @file cuda/driver/nested_loop_adapter.h
/// @brief Host adapter: a cuda::DriverNode contributes as one Component.
///
/// CUDA-free header. Inner phases stay inside `run_nested_erased`. This
/// type is a Component, never a DriverNode. Energy/force probes follow
/// `cuda::HF` on the inner node (same rule as `cpu::NestedLoopAdapter`).
#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "atomiverse/cuda/component/component.cuh"
#include "atomiverse/cuda/component/hf.cuh"
#include "atomiverse/cuda/driver/driver_node.h"

namespace ATV {
namespace cuda {

/// @brief Wraps a nested CUDA loop so the outer Plan only sees a Component.
class NestedLoopAdapter : public Component {
  public:
    /// @brief Borrow `inner` for the adapter lifetime.
    /// @param inner Nested loop node (caller retains ownership).
    /// @param budget Iteration budget; `< 0` uses `inner.default_budget()`.
    NestedLoopAdapter(DriverNode& inner, int budget = -1)
        : inner_(&inner), budget_(resolve_budget(inner, budget)) {}

    /// @brief Take ownership of `inner`.
    /// @param inner Non-null nested loop node.
    /// @param budget Iteration budget; `< 0` uses `inner.default_budget()`.
    /// @throws std::invalid_argument if `inner` is null.
    NestedLoopAdapter(std::unique_ptr<DriverNode> inner, int budget = -1)
        : owned_(std::move(inner)), inner_(owned_.get()) {
        if (!inner_) {
            throw std::invalid_argument(
                "cuda::NestedLoopAdapter: null DriverNode. Pass a node or "
                "use the borrow constructor.");
        }
        budget_ = resolve_budget(*inner_, budget);
    }

    /// @brief Own a concrete DriverNode rvalue that is not a Component.
    template <
        typename Node,
        std::enable_if_t<
            std::is_base_of_v<DriverNode, std::remove_reference_t<Node>>
                && !std::is_base_of_v<Component, std::remove_reference_t<Node>>
                && !std::is_lvalue_reference_v<Node>,
            int> = 0>
    NestedLoopAdapter(Node&& node, int budget = -1)
        : NestedLoopAdapter(std::unique_ptr<DriverNode>(
                                std::make_unique<std::remove_reference_t<Node>>(
                                    std::forward<Node>(node))),
                            budget) {}

    NestedLoopAdapter(const NestedLoopAdapter&) = delete;
    NestedLoopAdapter& operator=(const NestedLoopAdapter&) = delete;
    NestedLoopAdapter(NestedLoopAdapter&&) noexcept = default;
    NestedLoopAdapter& operator=(NestedLoopAdapter&&) noexcept = default;

    /// @brief Component name.
    /// @return `"nested_loop_adapter"`.
    std::string name() const override { return "nested_loop_adapter"; }

    /// @brief Energy probe. True when the inner node is an HF method that
    ///        provides energy.
    bool provides_energy() const override {
        if (const auto* hf = dynamic_cast<const HF*>(inner_))
            return hf->provides_energy();
        return false;
    }
    /// @brief Force probe. Follows the inner HF method (`with_gradient`).
    bool provides_force() const override {
        if (const auto* hf = dynamic_cast<const HF*>(inner_))
            return hf->provides_force();
        return false;
    }

    /// @brief The wrapped node.
    DriverNode* inner() noexcept { return inner_; }
    /// @brief The wrapped node.
    const DriverNode* inner() const noexcept { return inner_; }

    void build(const BuildContext& ctx) override {
        inner_->build_on(ctx.runtime, ctx.frame);
    }
    void launch(const LaunchContext& ctx) override {
        inner_->run_nested_erased(ctx.runtime, ctx.frame, budget_);
    }
    void destroy(const DestroyContext& ctx) override {
        inner_->destroy_on(ctx.runtime, ctx.frame);
    }

  private:
    static int resolve_budget(DriverNode& inner, int budget) {
        return budget < 0 ? inner.default_budget() : budget;
    }

    std::unique_ptr<DriverNode> owned_;
    DriverNode* inner_ = nullptr;
    int budget_ = 0;
};

} // namespace cuda
} // namespace ATV
