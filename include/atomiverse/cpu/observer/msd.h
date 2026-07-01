/// @file cpu/observer/msd.h
/// @brief CPU async mean-squared-displacement observer.
///
/// LAMMPS `compute msd`-shaped: the C++ engine *measures* (stages each captured
/// `FrameSnapshot` into a molrs-backed trajectory on the analysis worker
/// thread); molrs *analyzes* (`compute::MSD`, Direct mode, first frame =
/// reference). No analysis math lives in C++ — see molrs/src/compute/msd/.
///
/// Backend-uniform: regardless of the producing driver/backend, snapshots are
/// delivered to the CPU analysis executor's ring buffer (CUDA frames after
/// async D2H), so this observer is CPU-only by construction.
///
/// @see cpu::RDF (observer/rdf.h) for the async-observer pattern.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "atomiverse/cpu/async/async_observer.h"
#include "atomiverse/cpu/async/frame_snapshot.h"
#include "atomiverse/cpu/io/molrec.h"

namespace ATV {

/// @brief MD async mean-squared-displacement observer.
///
/// Accumulates positions from each captured `FrameSnapshot` into a molrs
/// trajectory; `finalize()` delegates to molrs `compute::MSD` and writes
/// MSD(t) to file.
///
/// Usage:
/// @code{.cpp}
/// cpu::md::MD sim;
/// sim.set_integrator(cpu::md::VelocityVerlet{0.001});
/// sim.set_pair(cpu::PairLJCut{0.0104, 3.4, 6.8});
/// MSD msd{/*dt_between_frames=*/0.1, /*freq=*/100, "msd.dat"};
/// sim.add(msd);
/// sim.build();  sim.run(frame, 100'000);
/// @endcode
class MSD : public cpu::md::AsyncObserver {
  public:
    /// @param dt          time between captured frames [fs] (= timestep * freq).
    /// @param freq        capture every `freq` steps (0 = every snapshot).
    /// @param output_file MSD(t) output path.
    explicit MSD(double dt, int freq = 0, std::string output_file = "msd.dat");

    std::string name() const override { return "msd_async"; }
    ComponentRole role() const override { return ComponentRole::Observer; }
    int frequency() const override { return freq_; }

    cpu::CaptureMask capture_mask() const override {
        cpu::CaptureMask m;
        m.positions = true;
        return m;
    }

    void on_snapshot(const cpu::FrameSnapshot& snap) override;
    void finalize() override;

    /// @brief The last computed MSD(t) curve (index 0 = 0). For inspection/tests.
    const std::vector<double>& curve() const { return msd_curve_; }
    /// @brief Einstein diffusion coefficient D from the last finalize() (NaN if
    /// not yet computed / too few frames).
    double diffusion() const { return d_; }
    int n_frames() const { return n_frames_; }
    const std::string& output_file() const { return output_file_; }

  private:
    double dt_;
    int freq_;
    std::string output_file_;
    cpu::MolRec traj_; ///< molrs-backed trajectory accumulator (move-only).
    int n_frames_ = 0;
    std::vector<double> msd_curve_;
    double d_ = 0.0; ///< Einstein diffusion coefficient from the last finalize().
};

} // namespace ATV
