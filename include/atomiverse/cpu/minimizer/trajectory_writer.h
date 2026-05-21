/// @file cpu/minimizer/trajectory_writer.h
/// @brief TrajectoryWriter -- append geometry frames to an XYZ trajectory.
///
/// Per-step recorder for geometry optimization or AIMD runs. Binds by default
/// to `minimizer::Phase::end_of_step`; compose on any phase via `on().add()`.
///
/// Usage:
/// @code
///   cpu::minimizer::Minimize opt;
///   opt.add(cpu::minimizer::LBFGS{/*...*/});
///   opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
///   opt.on(minimizer::Phase::end_of_step)
///      .add(cpu::minimizer::TrajectoryWriter{"aspirin_traj.xyz", /*stride=*/1});
///   opt.build();
///   opt.run(frame, 200);
/// @endcode
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Append the current geometry to an XYZ file each step.
///
/// File is opened on first launch (or per-call if the backend is stateless)
/// and written one frame at a time through molrs. Respects `stride` to thin
/// long trajectories.
///
/// MPI-aware: only rank 0 writes.
class TrajectoryWriter : public cpu::ComponentBase {
  public:
    /// @brief Construct.
    /// @param filename Output XYZ path (appended; existing file is overwritten
    ///                 on the first launched step).
    /// @param stride   Write every `stride` steps. Must be >= 1.
    explicit TrajectoryWriter(std::string filename, int stride = 1);

    std::string name() const override { return "trajectory_writer"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::minimizer::Phase::end_of_step), 0}};
    }

    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Number of frames written so far.
    int n_frames_written() const { return n_written_; }

  private:
    std::string filename_;
    int stride_;
    int n_written_ = 0;
    bool truncated_ = false;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
