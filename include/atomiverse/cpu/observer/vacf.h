/// @file cpu/observer/vacf.h
/// @brief CPU async velocity-autocorrelation-function observer.
///
/// The engine measures (accumulates per-frame velocities from the async snapshot
/// ring); molrs analyzes (`compute::VACF`, FFT-ACF). No analysis math in C++.
/// The VACF is the VDOS / Green-Kubo-diffusion input, complementing the Einstein
/// route in cpu::MSD.
///
/// @see cpu::MSD (observer/msd.h) for the async-observer pattern.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/async/async_observer.h"
#include "atomiverse/cpu/async/frame_snapshot.h"

namespace ATV {

/// @brief MD async velocity-autocorrelation observer.
///
/// Accumulates velocities from each captured `FrameSnapshot` into a row-major
/// `[n_frames, 3*n_atoms]` buffer; `finalize()` delegates to molrs `compute::VACF`
/// and writes the DOF-averaged VACF(t). Requires a constant atom count across the
/// captured window (single-rank / no migration).
class VACF : public cpu::md::AsyncObserver {
  public:
    /// @param dt         time between captured frames [fs].
    /// @param resolution max lag (0 = full: n_frames - 1).
    /// @param freq       capture every `freq` steps (0 = every snapshot).
    /// @param output_file VACF(t) output path.
    explicit VACF(double dt, int resolution = 0, int freq = 0,
                  std::string output_file = "vacf.dat");

    std::string name() const override { return "vacf_async"; }
    ComponentRole role() const override { return ComponentRole::Observer; }
    int frequency() const override { return freq_; }

    cpu::CaptureMask capture_mask() const override {
        cpu::CaptureMask m;
        m.velocities = true;
        return m;
    }

    void on_snapshot(const cpu::FrameSnapshot& snap) override;
    void finalize() override;

    /// @brief The last computed VACF(t) curve. For inspection/tests.
    const std::vector<double>& curve() const { return vacf_curve_; }
    int n_frames() const { return n_frames_; }

  private:
    double dt_;
    int resolution_;
    int freq_;
    std::string output_file_;
    int n_frames_ = 0;
    int n_atoms_ = 0;
    std::vector<double> vel_; ///< flat [n_frames, 3*n_atoms], blocked vx|vy|vz.
    std::vector<double> vacf_curve_;
};

} // namespace ATV
