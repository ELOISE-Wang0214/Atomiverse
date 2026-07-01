/// @file cpu/observer/vacf.cpp
/// @brief CPU async VACF observer — accumulates velocities, delegates the
/// autocorrelation to molrs.

#include "atomiverse/cpu/observer/vacf.h"

#include <cstddef>
#include <fstream>
#include <utility>

#include "atomiverse/cpu/io/molrec.h" // cpu::compute_vacf (molrs bridge)

namespace ATV {

VACF::VACF(double dt, int resolution, int freq, std::string output_file)
    : dt_(dt), resolution_(resolution), freq_(freq),
      output_file_(std::move(output_file)) {}

void VACF::on_snapshot(const cpu::FrameSnapshot& snap) {
    const int n = snap.n_owned;
    if (n <= 0) return;
    if (n_frames_ == 0) n_atoms_ = n;
    if (n != n_atoms_) return; // constant atom count required (single-rank MVP)

    // Row-major [n_frames, 3*n_atoms], blocked per frame: vx[0..n) vy[..) vz[..).
    // Snapshot velocity fields are double (widened from StorageReal at capture).
    vel_.insert(vel_.end(), snap.vx.begin(), snap.vx.begin() + n);
    vel_.insert(vel_.end(), snap.vy.begin(), snap.vy.begin() + n);
    vel_.insert(vel_.end(), snap.vz.begin(), snap.vz.begin() + n);
    ++n_frames_;
}

void VACF::finalize() {
    if (n_frames_ < 2 || n_atoms_ <= 0) return;
    const int n_dof = 3 * n_atoms_;
    const int res = resolution_ > 0 ? resolution_ : (n_frames_ - 1);
    vacf_curve_ = cpu::compute_vacf(vel_, n_frames_, n_dof, dt_, res); // molrs

    std::ofstream os(output_file_);
    os << "# VACF via molrs compute::VACF (DOF-averaged, FFT-ACF)\n";
    os << "# t[fs]  VACF   (frames=" << n_frames_ << ", atoms=" << n_atoms_
       << ")\n";
    for (std::size_t i = 0; i < vacf_curve_.size(); ++i)
        os << (dt_ * static_cast<double>(i)) << "  " << vacf_curve_[i] << '\n';
}

} // namespace ATV
