/// @file cpu/observer/msd.cpp
/// @brief CPU async MSD observer — stages snapshots into a molrs trajectory and
/// delegates the mean-squared-displacement compute to molrs.

#include "atomiverse/cpu/observer/msd.h"

#include <cstddef>
#include <fstream>
#include <utility>

namespace ATV {

MSD::MSD(double dt, int freq, std::string output_file)
    : dt_(dt), freq_(freq), output_file_(std::move(output_file)) {}

void MSD::on_snapshot(const cpu::FrameSnapshot& snap) {
    const int n = snap.n_owned;
    if (n <= 0) return;

    // Element-agnostic: MSD reads only positions. Dummy Z so molrs can build the
    // frame's element column; the value is never used by the MSD kernel.
    std::vector<std::int32_t> type_id(static_cast<std::size_t>(n), 1);

    std::array<double, 9> box{};
    const auto& cell = snap.box.cell();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box[static_cast<std::size_t>(i * 3 + j)] =
                static_cast<double>(cell(i, j));

    traj_.stage_raw(type_id, snap.x, snap.y, snap.z, box);
    traj_.commit_frame();
    ++n_frames_;
}

void MSD::finalize() {
    msd_curve_ = traj_.compute_msd(); // all analysis math in molrs
    // Einstein D from the windowed-MSD slope over the 20-80% lag window (avoids
    // the ballistic short-time and noisy long-time tails). molrs owns the math.
    d_ = traj_.compute_diffusion(dt_, /*dims=*/3, /*fit_lo=*/0.2, /*fit_hi=*/0.8);

    std::ofstream os(output_file_);
    os << "# MSD via molrs compute::MSD (Direct mode, frame-0 reference)\n";
    os << "# Einstein D (molrs, 20-80% lag fit) = " << d_
       << " Angstrom^2/fs   (frames=" << n_frames_ << ")\n";
    os << "# t[fs]  MSD[Angstrom^2]\n";
    for (std::size_t i = 0; i < msd_curve_.size(); ++i)
        os << (dt_ * static_cast<double>(i)) << "  " << msd_curve_[i] << '\n';
}

} // namespace ATV
