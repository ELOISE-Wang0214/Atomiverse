/// @file cpu/minimizer/final_geometry.cpp
/// @brief FinalGeometry -- rewrite a single-frame XYZ each step.
///        MPI-aware: collective owner gather, rank-0 I/O.

#include "atomiverse/cpu/minimizer/final_geometry.h"

#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/comm/gather.h"
#include "atomiverse/cpu/molrec.h"
#include "atomiverse/cpu/recorder/xyz.h"

namespace ATV {
namespace cpu {
namespace minimizer {

FinalGeometry::FinalGeometry(std::string filename)
    : filename_(std::move(filename)) {}

void FinalGeometry::launch(const cpu::LaunchContext& ctx) {
    // Multi-rank: gather the owner rows to rank 0 (collective — every rank
    // participates before the rank gate) so the snapshot holds the complete
    // tag-ordered geometry, not one rank's subdomain + stale tail.
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

    MolRec rec;
    rec.set_geometry(*src);
    rec.add_scalar("step", static_cast<double>(ctx.step));
    rec.add_scalar("energy_eV", ctx.frame.scalars.energy);
    rec.commit_frame();

    // XYZ::write truncates any existing file; every call replaces the prior
    // contents with a single-frame snapshot.
    XYZ writer(filename_);
    writer.write(rec);
    ++n_writes_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
