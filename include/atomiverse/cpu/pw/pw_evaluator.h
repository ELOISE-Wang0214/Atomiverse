/// @file cpu/pw/pw_evaluator.h
/// @brief PWEvaluator — plane-wave DFT energy component for CPU drivers.
///
/// Mirrors `cpu::scf::SCF`: holds a `PWSolver` inner solver
/// built once; each `launch()` runs the plane-wave SCF at the frame's
/// current geometry and writes the total energy to `frame.meta.get(kEnergy)`
/// in eV. Energy-only (no plane-wave forces yet) — plugs into
/// `cpu::Evaluate` for single-point work:
/// @code
///   cpu::pw::PWConfig cfg;
///   cfg.ecut_ry = 30.0;
///   cpu::pw::PWEvaluator pwe{cfg};
///   cpu::Evaluate eval;
///   eval.on(ATV::evaluate::Phase::evaluate).add_ref(pwe);
///   eval.build();
///   eval.run(frame);
///   // frame.energy [eV]; pwe.last_result().energy [Ry].
/// @endcode
#pragma once

#include "atomiverse/cpu/component/plane_wave.h"
#include "atomiverse/cpu/pw/pw_solver.h"

namespace ATV {
namespace cpu {
namespace pw {

/// @brief Plane-wave DFT total-energy evaluator component.
///
/// `role() == Evaluator`, `capabilities() == Energy`. Writes
/// `frame.meta.get(kEnergy)` [eV]; the Ry / Hartree breakdown stays on
/// `last_result()`.
class PWEvaluator : public cpu::PlaneWave {
  public:
    /// @brief Construct from a plane-wave configuration.
    /// @param config Cutoff, k-mesh, convergence, mixing parameters.
    explicit PWEvaluator(PWConfig config) : solver_(std::move(config)) {}

    std::string name() const override { return "pw_evaluator"; }
    /* role via Potential/kind base */

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Most recent plane-wave SCF result.
    const PWResult& last_result() const { return solver_.last_result(); }

    /// @brief Plane-wave configuration (for inspection).
    const PWConfig& config() const { return solver_.config(); }

  private:
    PWSolver solver_;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
