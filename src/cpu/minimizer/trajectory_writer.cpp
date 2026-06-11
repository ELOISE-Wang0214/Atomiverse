/// @file cpu/minimizer/trajectory_writer.cpp
/// @brief TrajectoryWriter -- molrs-backed XYZ append.
///        MPI-aware: collective owner gather, rank-0 I/O.

#include "atomiverse/cpu/minimizer/trajectory_writer.h"

#include <cstdio>

#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/comm/gather.h"
#include "atomiverse/cpu/io/molrec.h"
#include "atomiverse/cpu/recorder/xyz.h"

namespace ATV {
namespace cpu {
namespace minimizer {

TrajectoryWriter::TrajectoryWriter(std::string filename, int stride)
    : filename_(std::move(filename)), stride_(stride < 1 ? 1 : stride) {}

void TrajectoryWriter::launch(const cpu::LaunchContext& ctx) {
    // Stride gate first: ctx.step is identical on every rank, so the
    // collective below is still entered in lockstep.
    if ((ctx.step % stride_) != 0)
        return;

    // Multi-rank: gather the owner rows to rank 0 (collective — every rank
    // participates before the rank gate) so the written frame holds the
    // complete tag-ordered geometry, not one rank's subdomain + stale tail.
    const auto* dd = ctx.runtime.dd;
    const bool multi = dd && ctx.runtime.mpi_size > 1;
    cpu::Frame gathered{0};
    const cpu::Frame* src = &ctx.frame;
    if (multi) {
        const auto g = md::mpi::gather_to_root(ctx.frame, *dd);
        if (ctx.runtime.mpi_rank == 0) {
            gathered = md::mpi::to_frame(g);
            src = &gathered;
        }
    }
    if (ctx.runtime.mpi_rank != 0)
        return;

    // First write: truncate any pre-existing file so we do not append to it.
    if (!truncated_) {
        std::remove(filename_.c_str());
        truncated_ = true;
    }

    MolRec rec;
    rec.set_geometry(*src);
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
