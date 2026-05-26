/// @file cpu/backend.cpp
/// @brief cpu::Backend implementation.

#include "atomiverse/cpu/backend.h"
#include "atomiverse/core/log.h"
#include "atomiverse/cpu/analysis_executor.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#include <cstdio>
#include <cstdlib>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {

Backend::Backend() = default;
Backend::~Backend() { shutdown(); } // idempotent; releases cart_comm + executor
Backend::Backend(Backend&&) noexcept = default;
Backend& Backend::operator=(Backend&&) noexcept = default;

void Backend::init() {
#ifdef _OPENMP
    if (context_.n_threads > 0)
        omp_set_num_threads(context_.n_threads);
    else
        context_.n_threads = omp_get_max_threads();
#else
    context_.n_threads = 1;
#endif

    // Phase 5 Track B.6 — warn the user if they will pay NUMA / cache
    // penalties for letting the OS migrate threads. We check once per
    // process via a function-local static (the cpu::Backend may be
    // re-instantiated when Driver::shutdown() destroys it between
    // run() calls; the warning should fire once per process, not
    // once per backend init). Harmless on single-socket machines,
    // load-bearing on 2-socket / NUMA boxes.
    static bool s_omp_bind_warned = false;
    if (!s_omp_bind_warned && context_.n_threads > 1) {
        const char* bind = std::getenv("OMP_PROC_BIND");
        const bool unset = (bind == nullptr) || std::string(bind).empty();
        // Treat "false" as "unset" too — that's the OMP value users
        // explicitly opt out with, but it's still bad for NUMA.
        const bool unbinding = unset || std::string(bind) == "false"
                               || std::string(bind) == "FALSE";
        if (unbinding) {
            int my_rank = 0;
#ifdef ATV_HAS_MPI
            int mpi_init_check = 0;
            MPI_Initialized(&mpi_init_check);
            if (mpi_init_check)
                MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
            if (my_rank == 0) {
                std::fprintf(
                    stderr,
                    "[atv] warning: OMP_PROC_BIND is unset and "
                    "n_threads=%d > 1; consider 'export OMP_PROC_BIND=close "
                    "OMP_PLACES=cores' for stable hot-path throughput "
                    "(see benchmarks/md/reference/polyethylene/README.md "
                    "'Thread binding').\n",
                    context_.n_threads);
                std::fflush(stderr);
            }
            s_omp_bind_warned = true;
        }
    }

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

    if (!dd_)
        dd_ = std::make_unique<md::mpi::DomainDecomposition>();
    context_.dd = dd_.get();

#ifdef ATV_HAS_MPI
    if (context_.mpi_size > 1) {
        dd_->init_topology(MPI_COMM_WORLD);
    }
#endif

    if (!executor_) {
        executor_ = std::make_unique<AnalysisExecutor>();
        AnalysisExecutor::Config cfg;
        cfg.n_threads = 1;
        cfg.ring_size = 3;
        cfg.policy = BackpressurePolicy::drop_oldest;
        executor_->init(cfg, context_);
    }
}

void Backend::adopt_dd(
    std::unique_ptr<md::mpi::DomainDecomposition> dd) {
    dd_ = std::move(dd);
    context_.dd = dd_.get();
}

std::unique_ptr<md::mpi::DomainDecomposition> Backend::release_dd() {
    context_.dd = nullptr;
    return std::move(dd_);
}

void Backend::set_frame(cpu::Frame& frame) {
    frame_ = &frame;
    if (!dd_)
        return;

#ifdef ATV_HAS_MPI
    if (context_.mpi_size > 1 && dd_->cart_comm != MPI_COMM_NULL) {
        dd_->partition_atoms(frame);
        return;
    }
#endif
    dd_->bind(frame);
    dd_->compute_sublo_subhi(frame.box);
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
        if (slot.ops.empty())
            continue;
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
    for (int s = 0; s < n_steps; ++s)
        step(s);
    ATV_INFO("run complete");
}

void Backend::drain_analysis() {
    if (executor_)
        executor_->quiesce();
}

void Backend::destroy_components() {
    // Drain and join analysis workers first so any in-flight `on_snapshot`
    // completes against live subscribers, and `finalize()` runs on the main
    // thread before any component enters its `destroy()` path.
    if (executor_)
        executor_->shutdown();

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
    if (executor_) {
        executor_->shutdown(); // idempotent
        executor_.reset();
    }
    context_.executor = nullptr;

#ifdef ATV_HAS_MPI
    if (dd_) {
        dd_->release_comm();
    }
#endif

    frame_ = nullptr;
    plan_ = nullptr;
}

} // namespace cpu
} // namespace ATV
