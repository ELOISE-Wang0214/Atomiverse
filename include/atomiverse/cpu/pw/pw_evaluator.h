/// @file cpu/pw/pw_evaluator.h
/// @brief PWEvaluator — plane-wave DFT energy (+ optional FD forces) component.
///
/// Mirrors `cpu::scf::SCF`: holds a `PWSolver` inner solver built once;
/// each `launch()` runs the plane-wave SCF and writes energy [eV]. With
/// `with_gradient(true)` (default **false** — FD forces cost 6N SCFs),
/// ionic forces are filled via `pw_forces_fd` and written [eV/Å].
#pragma once

#include <stdexcept>
#include <utility>

#include "atomiverse/cpu/component/plane_wave.h"
#include "atomiverse/cpu/pw/pw_solver.h"

namespace ATV {
namespace cpu {
namespace pw {

/// @brief Plane-wave DFT total-energy evaluator component.
///
/// Kind `PlaneWave`. Writes `frame.meta.get(kEnergy)` [eV]; Ry/Hartree stay
/// on `last_result()`. Forces (when enabled) are central finite differences
/// of the SCF energy until analytic PW forces land.
class PWEvaluator : public cpu::PlaneWave {
  public:
    /// @brief Construct from a plane-wave configuration.
    /// @param config Cutoff, k-mesh, convergence, mixing parameters.
    explicit PWEvaluator(PWConfig config) : solver_(std::move(config)) {}

    std::string name() const override { return "pw_evaluator"; }

    /// @brief Enable/disable ionic force evaluation (FD of total energy).
    /// @param enabled When true, `launch()` also writes per-atom forces [eV/Å].
    PWEvaluator& with_gradient(bool enabled) & noexcept {
        compute_forces_ = enabled;
        return *this;
    }
    PWEvaluator&& with_gradient(bool enabled) && noexcept {
        compute_forces_ = enabled;
        return std::move(*this);
    }

    /// @brief True when forces will be written on launch.
    bool compute_forces() const noexcept { return compute_forces_; }

    /// @brief FD step for forces [Å]. Default 1e-3.
    PWEvaluator& set_fd_step(double h_ang) & {
        if (h_ang <= 0.0)
            throw std::invalid_argument("PWEvaluator: fd step must be > 0");
        fd_step_ = h_ang;
        return *this;
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Most recent plane-wave SCF result.
    const PWResult& last_result() const { return solver_.last_result(); }

    /// @brief Plane-wave configuration (for inspection).
    const PWConfig& config() const { return solver_.config(); }

  private:
    PWSolver solver_;
    bool compute_forces_ = false;
    double fd_step_ = 1e-3;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
