/// @file cuda/minimizer/cg.cuh
/// @brief cuda::minimizer::CG — Polak-Ribière+ conjugate-gradient
///        Optimizer component.
///
/// CUDA mirror of `cpu::minimizer::CG`: a hot-swappable alternative to
/// `cuda::minimizer::LBFGS` / `SteepestDescent` inside
/// `cuda::minimizer::Minimize`. Same outer contract as the CPU class —
/// binds `search_direction` (computes and applies the step) and
/// `convergence_check` (reports `outer_status.converged` to the driver).
///
/// Algorithm (Polak-Ribière+):
///   d_0 = -g_0
///   For k ≥ 1:
///     β = max(0, (g_k · (g_k − g_{k−1})) / (g_{k−1} · g_{k−1}))
///     d_k = -g_k + β · d_{k−1}
///   x_{k+1} = x_k + α · d_k
///
/// Execution-model classification: the `search_direction` phase is
/// **GPU-pure compute** — β is produced and consumed entirely on the
/// device (a fused dot-product reduction feeds a single-thread β kernel
/// whose result the step kernel reads by pointer, the same
/// device-scalar chaining as `cuda::md::NoseHoover`'s scale factor).
/// The `convergence_check` phase is the shared `ConvCheck`
/// host-dependent bridge (one pinned D2H + sync per iteration).
///
/// Convergence rule is identical to `LBFGS` (and to both CPU
/// optimizers): max/RMS force, max/RMS displacement, energy change.
///
/// @note Like its CPU mirror, CG applies a FIXED step `x += α₀·d` with
/// no per-atom displacement cap (`LBFGS` caps at 0.1 Å, CG does not).
/// α₀ must therefore be sized to the stiffness of the landscape:
/// `α₀ · |F|_max` should stay well below the shortest interaction
/// length scale, or the trajectory diverges. For stiff starts prefer
/// `LBFGS` or shrink α₀.
///
/// References:
///   Polak & Ribière, ESAIM Math. Modelling Numer. Anal. 3, 35 (1969)
///   Nocedal & Wright, Numerical Optimization, Ch. 5
///
/// @see cuda::minimizer::Minimize for the driver this component plugs into
/// @see cpu::minimizer::CG for the CPU mirror
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component/component.cuh"
#include "atomiverse/cuda/minimizer/conv_check.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace minimizer {

/// @brief Polak-Ribière+ conjugate-gradient Optimizer.
///
/// Optimizer state (previous positions / gradient / direction) lives in
/// planar 3N `work_t` device buffers allocated in `build()` and released
/// in `destroy()` — no allocation in `launch()`. Buffers are zeroed at
/// build, which reproduces the CPU step-0 behaviour exactly: forces are
/// zero at step 0, `g_prev · g_prev = 0` forces β = 0, and the first
/// moving step is plain steepest descent.
class CG : public ComponentBase {
  public:
    /// @brief Construct CG optimizer (CPU-default tolerances).
    /// @param alpha0    Step size [Angstrom^2/eV] (default 0.01).
    /// @param ftol      Max force convergence tolerance [eV/Angstrom] (default 0.05).
    /// @param rms_ftol  RMS force tolerance [eV/Angstrom] (default 0.03).
    /// @param max_dtol  Max displacement tolerance [Angstrom] (default 0.002).
    /// @param rms_dtol  RMS displacement tolerance [Angstrom] (default 0.001).
    /// @param etol      Energy change tolerance [eV] (default 1e-5).
    explicit CG(double alpha0 = 0.01, double ftol = 0.05,
                double rms_ftol = 0.03, double max_dtol = 0.002,
                double rms_dtol = 0.001, double etol = 1e-5)
        : alpha0_(alpha0), tols_{ftol, rms_ftol, max_dtol, rms_dtol, etol} {}

    /// @brief Component name.
    /// @return "cg".
    std::string name() const override { return "cg"; }

    /// @brief Structural role of this component.
    /// @return ComponentRole::Optimizer.
    ComponentRole role() const override { return ComponentRole::Optimizer; }

    /// @brief Default bindings: search_direction + convergence_check.
    std::vector<DefaultBinding> default_bindings() const override {
        return {
            {static_cast<int>(ATV::minimizer::Phase::search_direction), 0},
            {static_cast<int>(ATV::minimizer::Phase::convergence_check), 0}};
    }

    /// @brief Allocate + zero the optimizer state buffers.
    void build(const BuildContext& ctx) override;

    /// @brief Dispatch the bound phase (search_direction / convergence_check).
    void launch(const LaunchContext& ctx) override;

    /// @brief Release all device buffers.
    void destroy(const DestroyContext& ctx) override;

    /// @brief Check if optimization has converged.
    bool converged() const { return converged_; }

    /// @brief Max force magnitude at the latest convergence check.
    double max_force() const { return conv_.max_force(); }

    /// @brief Per-step optimization trace.
    const std::vector<OptRecord>& trace() const { return trace_; }

  private:
    double alpha0_; ///< Step size [Angstrom^2/eV].
    ConvTols tols_; ///< Five-criteria convergence tolerances.

    int n_ = 0;     ///< Atom count.
    int step_ = 0;  ///< Optimizer step counter.
    bool converged_ = false;

    // Planar [x:0..n, y:n..2n, z:2n..3n] device state, work_t precision
    // (position/gradient differences of nearly-equal values would lose
    // too much precision in real_t storage).
    work_t* d_xprev_ = nullptr; ///< Pre-step positions (3N).
    work_t* d_gprev_ = nullptr; ///< Previous gradient (3N).
    work_t* d_dprev_ = nullptr; ///< Previous search direction (3N).
    work_t* d_pr_ = nullptr;    ///< Device {num, den, beta} scalars.

    ConvCheck conv_;
    std::vector<OptRecord> trace_;
};

} // namespace minimizer
} // namespace cuda
} // namespace ATV
