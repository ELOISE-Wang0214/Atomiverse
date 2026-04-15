/// @file cpu/minimizer/final_geometry.cpp
/// @brief FinalGeometry -- rewrite a single-frame XYZ each step.

#include "atomiverse/cpu/minimizer/final_geometry.h"

#include "atomiverse/cpu/molrec.h"
#include "atomiverse/cpu/recorder/xyz.h"

namespace ATV {
namespace cpu {
namespace minimizer {

FinalGeometry::FinalGeometry(std::string filename)
    : filename_(std::move(filename)) {}

void FinalGeometry::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_rank != 0) return;

    MolRec rec;
    rec.set_geometry(ctx.frame);
    rec.add_scalar("step", static_cast<double>(ctx.step));
    rec.add_scalar("energy_eV", ctx.frame.energy);
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
