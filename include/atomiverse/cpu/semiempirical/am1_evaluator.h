/// @file cpu/semiempirical/am1_evaluator.h
/// @brief AM1 Evaluator component for Atomiverse drivers.
///
/// Mirrors the existing SCF/PW/GPW shape: `AM1Evaluator` is the Component that
/// plugs into outer drivers, while `AM1Solver` is the inner solver it owns.
#pragma once

#include <string>
#include <utility>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/semiempirical/am1_solver.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

/// @brief AM1 semiempirical Evaluator component (wraps an `AM1Solver`).
class AM1Evaluator : public cpu::ComponentBase {
  public:
    /// @brief Construct with an AM1 configuration.
    /// @param config Solver configuration (charge, SCF tolerance, ...).
    explicit AM1Evaluator(AM1Config config = {}) : solver_(std::move(config)) {}

    /// @brief Component name.
    /// @return "am1_evaluator".
    std::string name() const override { return "am1_evaluator"; }
    /// @brief Structural role of this component.
    /// @return ComponentRole::Evaluator.
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    /// @brief Frame fields this component writes each launch.
    /// @return `Energy | Charge | DensityMatrix`.
    Capability capabilities() const override {
        return Capability::Energy | Capability::Charge
               | Capability::DensityMatrix;
    }
    // Deliberately no default_bindings(): like SCFEvaluator/PWEvaluator/
    // GPWEvaluator, this evaluator is reusable across workflows (Evaluate,
    // MD, Minimize), so the target phase must be chosen explicitly with
    // .on(phase).add() — a hardcoded evaluate::Phase would misbind when
    // added to another driver.

    /// @brief Build the owned solver against the reference frame (one-time).
    /// @param ctx Build context (runtime services + frame).
    void build(const cpu::BuildContext& ctx) override;
    /// @brief Run the AM1 SCF cycle; writes energy [eV] and per-atom
    ///        charges [e] onto the frame.
    /// @param ctx Launch context (runtime services + frame + stage/step).
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Converged result of the most recent launch.
    /// @return The solver's last `AM1Result` (energies, charges, density).
    const AM1Result& last_result() const { return solver_.last_result(); }
    /// @brief Solver configuration in effect.
    /// @return The `AM1Config` this evaluator was constructed with.
    const AM1Config& config() const { return solver_.config(); }

  private:
    AM1Solver solver_;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
