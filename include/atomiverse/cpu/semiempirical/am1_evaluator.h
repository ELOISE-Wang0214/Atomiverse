/// @file cpu/semiempirical/am1_evaluator.h
/// @brief AM1 Evaluator component for Atomiverse drivers.
#pragma once

#include <string>
#include <utility>

#include "atomiverse/cpu/component/semiempirical.h"
#include "atomiverse/cpu/semiempirical/am1_solver.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

/// @brief AM1 semiempirical Evaluator (owns concrete `AM1` method).
class AM1Evaluator : public cpu::Semiempirical {
  public:
    AM1Evaluator() = default;
    /// @brief Own a concrete AM1 method (parameters already on the method).
    explicit AM1Evaluator(AM1 method) : solver_(std::move(method)) {}
    /// @brief Construct method with SCF parameters (same args as `AM1`).
    AM1Evaluator(int charge, int multiplicity, int max_iter = 200,
                 double density_tol = 1.0e-8, double energy_tol = 1.0e-10)
        : solver_(charge, multiplicity, max_iter, density_tol, energy_tol) {}

    std::string name() const override { return "am1_evaluator"; }

    AM1Evaluator& with_gradient(bool enabled) & noexcept {
        compute_force_ = enabled;
        return *this;
    }
    AM1Evaluator&& with_gradient(bool enabled) && noexcept {
        compute_force_ = enabled;
        return std::move(*this);
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    const AM1Result& last_result() const { return solver_.last_result(); }
    const AM1& method() const { return solver_; }
    AM1& method() { return solver_; }
    bool gradient_enabled() const noexcept { return compute_force_; }

  private:
    AM1 solver_;
    bool compute_force_ = true;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
