/// @file cuda/minimizer/conv_check.cuh
/// @brief ConvCheck — shared multi-criteria convergence machinery for the
///        CUDA Optimizer components (`cuda::minimizer::CG` / `LBFGS`).
///
/// Both optimizers report convergence under the SAME five criteria as
/// their CPU counterparts (max/RMS force, max/RMS displacement, energy
/// change — see `cpu::minimizer::LBFGS`). Per the CLAUDE.md
/// "no duplication across methods" rule the shared skeleton lives here
/// once, owned by each Optimizer component (it is NOT a Plan component
/// itself — exactly as the CPU optimizers carry their convergence logic
/// inside their own `convergence_check`-phase branch).
///
/// Execution-model classification: **host-dependent compute**. The
/// metric reduction is a GPU-pure fused kernel, but the convergence
/// verdict must reach the host every iteration — `cuda::OptimizeDriver`
/// polls the host-side `Context::outer_status` after each `step_once()`
/// to break the outer loop. `check()` therefore issues one pinned D2H +
/// `cudaStreamSynchronize` per minimizer iteration. This is the
/// documented exception for the minimizer category (an iteration is
/// dominated by the force evaluation; the deferred `AsyncBridge` ring
/// is not applicable because its K-step lag would overrun the
/// convergence point).
///
/// @see cpu/minimizer/lbfgs.h for the CPU convergence rule this mirrors
/// @see cuda::OptimizeDriver for the outer loop that polls outer_status
#pragma once

#include <vector>

#include "atomiverse/cpu/minimizer/lbfgs.h" // cpu::minimizer::OptRecord
#include "atomiverse/cuda/component/component.cuh"
#include "atomiverse/cuda/frame/types.cuh"

namespace ATV {
namespace cuda {
namespace minimizer {

/// @brief Per-step optimization record — shared with the CPU backend so
///        diagnostics tooling reads one type.
using OptRecord = ATV::cpu::minimizer::OptRecord;

/// @brief The five convergence tolerances (CPU-default units).
struct ConvTols {
    double ftol = 0.05;      ///< Max force tolerance [eV/Angstrom].
    double rms_ftol = 0.03;  ///< RMS force tolerance [eV/Angstrom].
    double max_dtol = 0.002; ///< Max displacement tolerance [Angstrom].
    double rms_dtol = 0.001; ///< RMS displacement tolerance [Angstrom].
    double etol = 1e-5;      ///< Energy change tolerance [eV].
};

/// @brief Fused force/displacement metric reduction + host-side
///        five-criteria convergence verdict.
///
/// Lifecycle mirrors the owning component: `build()` allocates the
/// device accumulators + pinned host mirror, `check()` runs once per
/// `convergence_check` phase, `destroy()` releases.
///
/// Copyable only before `build()` (raw device pointers; components are
/// copied into the Plan at registration time, before any allocation).
class ConvCheck {
  public:
    /// @brief Allocate device metric accumulators + pinned host mirror.
    /// @param n_atoms Atom count (fixes the RMS normalizations).
    void build(int n_atoms);

    /// @brief Release all resources allocated in build().
    void destroy();

    /// @brief Reduce metrics, sync, evaluate the five criteria, record.
    ///
    /// Launches the fused metric kernel, copies the four metrics + the
    /// frame potential energy D2H into pinned memory, synchronizes the
    /// compute stream, evaluates the CPU-identical convergence rule,
    /// appends an OptRecord to @p trace and — on convergence — sets
    /// `ctx.runtime.outer_status.{converged, reason}`.
    ///
    /// Displacement and energy-change criteria are skipped at
    /// `step == 0`, exactly as on the CPU. Because `check()` ends with
    /// a stream synchronize, any D2H copies the caller enqueued on the
    /// same stream beforehand are complete when it returns (LBFGS
    /// batches its line-search scalars this way).
    ///
    /// @param ctx     Launch context (frame, stream, outer_status).
    /// @param x_prev  Pre-step positions, planar [x:0..n, y:n..2n,
    ///                z:2n..3n] device buffer (saved by the optimizer's
    ///                search_direction step). May be stale at step 0 —
    ///                it is not read then.
    /// @param step    Optimizer step counter (post search_direction
    ///                increment, mirroring the CPU components).
    /// @param tols    The five tolerances.
    /// @param reason  outer_status.reason on convergence.
    /// @param trace   Per-step diagnostics sink (host).
    /// @return true iff all five criteria are satisfied.
    bool check(const LaunchContext& ctx, const work_t* x_prev, int step,
               const ConvTols& tols, const char* reason,
               std::vector<OptRecord>& trace);

    /// @brief Potential energy read back by the most recent check() [eV].
    double last_energy() const { return last_energy_; }

    /// @brief Max force magnitude from the most recent check() [eV/Angstrom].
    double max_force() const { return max_force_; }

  private:
    int n_ = 0;                  ///< Atom count.
    work_t* d_metrics_ = nullptr; ///< Device {max_f, sum_f2, max_disp, sum_d2}.
    work_t* h_metrics_ = nullptr; ///< Pinned mirror of d_metrics_.
    double* h_energy_ = nullptr;  ///< Pinned mirror of frame.scalars->energy.
    double energy_prev_ = 0;      ///< Energy at the previous check [eV].
    double last_energy_ = 0;      ///< Energy at the latest check [eV].
    double max_force_ = 0;        ///< Max |F| at the latest check [eV/Angstrom].
};

} // namespace minimizer
} // namespace cuda
} // namespace ATV
