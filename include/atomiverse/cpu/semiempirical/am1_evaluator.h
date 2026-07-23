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
    explicit AM1Evaluator(AM1 method) : solver_(std::move(method)) {}
    explicit AM1Evaluator(AM1Config config) {
        solver_.charge = config.charge;
        solver_.multiplicity = config.multiplicity;
        solver_.max_iter = config.max_iter;
        solver_.density_tol = config.density_tol;
        solver_.energy_tol = config.energy_tol;
    }
    AM1Config config() const {
        AM1Config c;
        c.charge = solver_.charge;
        c.multiplicity = solver_.multiplicity;
        c.max_iter = solver_.max_iter;
        c.density_tol = solver_.density_tol;
        c.energy_tol = solver_.energy_tol;
        return c;
    }

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
