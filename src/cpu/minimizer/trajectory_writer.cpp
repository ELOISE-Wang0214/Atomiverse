/// @file cpu/minimizer/trajectory_writer.cpp
/// @brief TrajectoryWriter -- molrs-backed XYZ append.
///        MPI-aware: collective owner gather, rank-0 I/O.

#include "atomiverse/cpu/minimizer/trajectory_writer.h"

#include <spdlog/fmt/fmt.h>

#include "atomiverse/cpu/comm/gather_for_io.h"

namespace ATV {
namespace cpu {
namespace minimizer {

TrajectoryWriter::TrajectoryWriter(std::string filename, int stride)
    : filename_(std::move(filename)), stride_(stride < 1 ? 1 : stride),
      writer_(filename_) {}

void TrajectoryWriter::launch(const cpu::LaunchContext& ctx) {
    // Stride gate first: ctx.step is identical on every rank, so the
    // collective below is still entered in lockstep.
    if ((ctx.step % stride_) != 0)
        return;

    // Multi-rank: collective owner gather (every rank participates before
    // the rank gate) so the written frame holds the complete tag-ordered
    // geometry; returns nullptr on non-root ranks.
    const cpu::Frame* src =
        md::mpi::gather_for_root_io(ctx.frame, ctx.runtime, gathered_);
    if (!src)
        return;

    // Per-frame step + energy on the extxyz comment line (the pre-refactor
    // MolRec scalars, restored through the molrs meta hook).
    const cpu::XYZMeta meta{
        {"step", std::to_string(ctx.step)},
        {"energy_eV", fmt::format("{:.10g}", ctx.frame.scalars.energy)}};

    // First write truncates any pre-existing file; later writes append.
    if (!truncated_) {
        writer_.write(*src, meta);
        truncated_ = true;
    } else {
        writer_.append(*src, meta);
    }
    ++n_written_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
