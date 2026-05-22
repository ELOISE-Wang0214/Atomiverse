/// @file cpu/minimizer/trajectory_writer.cpp
/// @brief TrajectoryWriter -- molrs-backed XYZ append.

#include "atomiverse/cpu/minimizer/trajectory_writer.h"

#include <cstdio>

#include "atomiverse/cpu/molrec.h"
#include "atomiverse/cpu/recorder/xyz.h"

namespace ATV {
namespace cpu {
namespace minimizer {

TrajectoryWriter::TrajectoryWriter(std::string filename, int stride)
    : filename_(std::move(filename)), stride_(stride < 1 ? 1 : stride) {}

void TrajectoryWriter::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_rank != 0)
        return;
    if ((ctx.step % stride_) != 0)
        return;

    // First write: truncate any pre-existing file so we do not append to it.
    if (!truncated_) {
        std::remove(filename_.c_str());
        truncated_ = true;
    }

    MolRec rec;
    rec.set_geometry(ctx.frame);
    rec.add_scalar("step", static_cast<double>(ctx.step));
    rec.add_scalar("energy_eV", ctx.frame.scalars.energy);
    rec.commit_frame();

    XYZ writer(filename_);
    writer.append(rec);
    ++n_written_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
