/// @file cpu/backend.cpp
/// @brief cpu::Backend implementation.

#include "atomiverse/cpu/backend.h"
#include "atomiverse/core/log.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {

void Backend::init() {
#ifdef _OPENMP
    if (context_.n_threads > 0)
        omp_set_num_threads(context_.n_threads);
    else
        context_.n_threads = omp_get_max_threads();
#else
    context_.n_threads = 1;
#endif

#ifdef ATV_HAS_MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) {
        MPI_Comm_rank(MPI_COMM_WORLD, &context_.mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &context_.mpi_size);
    }
    ATV_INFO("CPU backend initialized: {} threads, MPI rank {}/{}",
             context_.n_threads, context_.mpi_rank, context_.mpi_size);
#else
    ATV_INFO("CPU backend initialized: {} threads", context_.n_threads);
#endif
}

void Backend::build_components() {
    BuildContext ctx{context_, *frame_};
    for (auto& slot : plan_->slots)
        for (auto& op : slot.ops) {
            ATV_DEBUG("building component '{}'", op.component->name());
            op.component->build(ctx);
        }
    ATV_INFO("all components built");
}

void Backend::step(int step_idx) {
    for (auto& slot : plan_->slots) {
        if (slot.ops.empty()) continue;
        for (auto& op : slot.ops) {
            if (op.frequency > 0 && step_idx % op.frequency != 0)
                continue;
            LaunchContext ctx{context_, *frame_, slot.phase, step_idx};
            op.component->launch(ctx);
        }
    }
}

void Backend::run(int n_steps) {
    ATV_INFO("running {} steps ({} atoms)", n_steps, frame_->atoms.size());
    for (int s = 0; s < n_steps; ++s) step(s);
    ATV_INFO("run complete");
}

void Backend::destroy_components() {
    DestroyContext ctx{context_, *frame_};
    for (int s = static_cast<int>(plan_->slots.size()) - 1; s >= 0; --s) {
        auto& slot = plan_->slots[s];
        for (int o = static_cast<int>(slot.ops.size()) - 1; o >= 0; --o) {
            ATV_DEBUG("destroying component '{}'",
                      slot.ops[o].component->name());
            slot.ops[o].component->destroy(ctx);
        }
    }
}

void Backend::shutdown() {
    ATV_DEBUG("CPU backend shutting down");
    frame_ = nullptr;
    plan_ = nullptr;
}

} // namespace cpu
} // namespace ATV
