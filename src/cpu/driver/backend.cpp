/// @file cpu/driver/backend.cpp
/// @brief cpu::Backend implementation.

#include "atomiverse/cpu/driver/backend.h"
#include "atomiverse/core/log.h"
#include "atomiverse/cpu/async/analysis_executor.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/driver/step_timing.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>

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
                    "(see benchmarks/md/pe/validation/reference/SOURCES.md "
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
#endif

    // Rank-aware logging auto-init: gates console to rank 0 (so 64-rank runs
    // don't spam stdout) and honors ATV_LOG_* env knobs. Skipped if the
    // application already called `ATV::log_init(...)` explicitly (that wins).
    if (!ATV::log_initialized()) {
        LogConfig lc;
        lc.mpi_rank = context_.mpi_rank;
        lc.mpi_size = context_.mpi_size;
        ATV::log_init(lc);
    }

#ifdef ATV_HAS_MPI
    ATV_INFO("CPU backend initialized: {} threads, MPI rank {}/{}",
             context_.n_threads, context_.mpi_rank, context_.mpi_size);
#else
    ATV_INFO("CPU backend initialized: {} threads", context_.n_threads);
#endif

    // Per-component timing is ON BY DEFAULT (near-zero overhead: a couple of
    // clock reads per component per step) so every driver can emit the built-in
    // end-of-run time+memory summary, the way LAMMPS always reports its loop
    // breakdown. Opt OUT with ATV_TIMER=0/off for the rare latency-critical
    // micro-measurement. Publish the registry so components can record
    // sub-scopes through Context::timer.
    // Re-apply the remembered preference (default true) so a pre-run
    // enable_timing(false) sticks across this init; ATV_TIMER overrides.
    timing_.set_enabled(timing_pref_);
    if (const char* tv = std::getenv("ATV_TIMER")) {
        const std::string v = tv;
        if (v == "0" || v == "false" || v == "FALSE" || v == "off")
            timing_.set_enabled(false);
        else if (v == "1" || v == "true" || v == "TRUE" || v == "on")
            timing_.set_enabled(true);
    }
    context_.timer = &timing_;

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
        cfg.ring_size = 3; // global ring depth K
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
        // First call: dd has no prior partition (n_local==0) → do the full
        // initial scan that classifies the global frame into per-rank slabs.
        //
        // Subsequent calls (e.g. `adopt_dd` across per-stage driver
        // rebuilds in polyethylene_tension): dd already holds the previous
        // stage's partition state. Re-running `partition_atoms` here would
        // iterate only `[0, n_local)` and DROP any atom that has drifted
        // across a slab boundary since the last Migrate — those
        // atoms vanish (no rank picks them up) and the global count
        // silently decreases (observed: 10000 → 9978 after stage 1, then
        // catastrophic NPT response in stage 2 step ~100-300 once the
        // density mismatch destabilizes the barostat). Skip the
        // re-partition; refresh slab bounds against the (possibly
        // box-changed) current box and let the next Migrate absorb
        // any small drift exactly as it does mid-stage.
        if (dd_->n_local == 0) {
            // Re-derive the processor grid box-aware now that the box is
            // known (init() ran before any frame was bound, so dims were
            // chosen for an unknown box). Must precede the initial
            // partition: builder-side pre-filtered frames use the same
            // `choose_dims`, and the two grids must agree on ownership.
            dd_->init_topology(MPI_COMM_WORLD, frame.box.length(0),
                               frame.box.length(1), frame.box.length(2));
            dd_->partition_atoms(frame);
        } else {
            dd_->compute_sublo_subhi(frame.box);
            // Carry-over preserves the partition (no re-scan — see above), but
            // the allocation bound must still grow if this frame is larger than
            // the one the partition was sized for (e.g. a reused backend
            // rebound to a bigger system). Without this, HaloSend /
            // Migrate reserve their per-instance buffers against a stale,
            // too-small `capacity` and can overflow. Grow-only: never shrink
            // below the existing ghost slack.
            if (frame.atoms.size() > dd_->capacity)
                dd_->capacity = frame.atoms.size();
        }
        release_tag_maps_if_unused(frame);
        return;
    }
#endif
    dd_->bind(frame);
    dd_->compute_sublo_subhi(frame.box);
    release_tag_maps_if_unused(frame);
}

void Backend::release_tag_maps_if_unused(cpu::Frame& frame) {
    if (frame.atoms.has_bonded())
        return;
    // No bonded topology → nothing resolves global tags at run time
    // (ResolveView and the add_* topology loaders are the only readers).
    // The loaders may have populated `atoms.map` before topology turned
    // out to be empty (materialize / lammps_data / benches call map_set
    // unconditionally); drop it and the DD mirror. ~50 B/atom each —
    // ~400 MB dead weight on an 8M-atom single-rank LJ run.
    frame.atoms.map_release();
    if (dd_) {
        std::unordered_map<TagInt, int>{}.swap(dd_->tag_to_local);
    }
}

void Backend::build_components() {
    BuildContext ctx{context_, *frame_};
    // Register one timing handle per op (component name + phase + category
    // rollup), precomputed so the hot step loop never re-classifies / allocates.
    op_handle_.assign(plan_->slots.size(), {});
    for (std::size_t s = 0; s < plan_->slots.size(); ++s) {
        auto& slot = plan_->slots[s];
        op_handle_[s].reserve(slot.ops.size());
        for (auto& op : slot.ops) {
            const std::string nm = op.component->name();
            ATV_DEBUG("building component '{}'", nm);
            op.component->build(ctx);
            op_handle_[s].push_back(timing_.register_op(
                TimerKey{nm, slot.phase,
                         static_cast<int>(md::StepTiming::classify(nm)), ""}));
        }
    }
    ATV_INFO("all components built");
}

void Backend::step(int step_idx) {
    const bool timing = timing_.enabled();
    const bool have_handles = op_handle_.size() == plan_->slots.size();
    const auto step_t0 = timing ? std::chrono::steady_clock::now()
                                : std::chrono::steady_clock::time_point{};
    for (std::size_t s = 0; s < plan_->slots.size(); ++s) {
        auto& slot = plan_->slots[s];
        if (slot.ops.empty())
            continue;
        for (std::size_t o = 0; o < slot.ops.size(); ++o) {
            auto& op = slot.ops[o];
            if (op.frequency > 0 && step_idx % op.frequency != 0)
                continue;
            LaunchContext ctx{context_, *frame_, slot.phase, step_idx};
            if (timing && have_handles) {
                const auto t0 = std::chrono::steady_clock::now();
                op.component->launch(ctx);
                const auto t1 = std::chrono::steady_clock::now();
                timing_.record(op_handle_[s][o],
                               std::chrono::duration<double>(t1 - t0).count());
            } else {
                op.component->launch(ctx); // zero-overhead path when off
            }
        }
    }
    if (timing) {
        const auto step_t1 = std::chrono::steady_clock::now();
        timing_.set_loop(
            std::chrono::duration<double>(step_t1 - step_t0).count(),
            dd_ ? dd_->n_local : static_cast<int>(frame_->atoms.size()));
        timing_.tick_step();
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

std::vector<std::pair<std::string, std::size_t>>
Backend::collect_component_memory() const {
    std::vector<std::pair<std::string, std::size_t>> out;
    // Unique components across all plan slots (a component may bind to several
    // phases) — dedup by pointer so its memory is counted once.
    if (plan_) {
        std::vector<const cpu::ComponentBase*> seen;
        for (const auto& slot : plan_->slots) {
            for (const auto& op : slot.ops) {
                const cpu::ComponentBase* c = op.component;
                if (!c
                    || std::find(seen.begin(), seen.end(), c) != seen.end())
                    continue;
                seen.push_back(c);
                const std::size_t b = c->memory_bytes();
                if (b)
                    out.emplace_back(c->name(), b);
            }
        }
    }
    // Per-atom SoA + topology + tag map (the dominant "Atoms" bucket).
    if (frame_)
        out.emplace_back("[atoms]", frame_->atoms.memory_bytes());
    // Domain-decomposition tag→local map + ghost-tag list (Comm bucket).
    if (dd_) {
        std::size_t b = dd_->tag_to_local.size() * 32
                        + dd_->tag_to_local.bucket_count() * sizeof(void*)
                        + dd_->ghost_tags.capacity() * sizeof(std::int64_t);
        out.emplace_back("[domain]", b);
    }
    return out;
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
