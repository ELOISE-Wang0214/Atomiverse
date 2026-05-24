/// @file driver.cu
/// @brief CUDA Backend implementation + Driver explicit instantiations.

#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/backend.cuh"
#include "atomiverse/cuda/backend_traits.cuh"
#include "atomiverse/cuda/error.cuh"
#include "atomiverse/cuda/log_callback.h"
#include "atomiverse/cuda/timer.cuh"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/facade/driver_impl.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

// ── Backend ────────────────────────────────────────────────────────

void Backend::init() {
    cuda::register_log_callback();
    cudaStream_t main, defer;
    ATV_CUDA_CHECK(cudaStreamCreate(&main));
    ATV_CUDA_CHECK(cudaStreamCreate(&defer));
    context_.stream = main;
    context_.defer_stream = defer;
}

void Backend::materialize(const cpu::Frame& frame) {
    frame_ = cpu::to_cuda(frame);
}

void Backend::build_components() {
    BuildContext ctx{context_, frame_};
    for (auto& slot : plan_->slots)
        for (auto& op : slot.ops) {
            if (timer_)
                timer_->impl()->begin(static_cast<cudaStream_t>(context_.stream));
            op.component->build(ctx);
            if (timer_) {
                timer_->impl()->end(static_cast<cudaStream_t>(context_.stream),
                            slot.phase, op.component->name() + "/build");
            }
        }
    cudaStreamSynchronize(static_cast<cudaStream_t>(context_.stream));
    if (timer_)
        timer_->impl()->flush();
}

void Backend::step(int step_idx) {
    for (auto& slot : plan_->slots) {
        if (slot.ops.empty())
            continue;
        for (auto& op : slot.ops) {
            if (op.frequency > 0 && step_idx % op.frequency != 0)
                continue;
            auto stream = static_cast<cudaStream_t>(
                (op.scheduling == Scheduling::deferrable)
                    ? context_.defer_stream
                    : context_.stream);
            LaunchContext ctx{context_, frame_, slot.phase, step_idx};
            if (timer_)
                timer_->impl()->begin(stream);
            op.component->launch(ctx);
            if (timer_)
                timer_->impl()->end(stream, slot.phase, op.component->name());
        }
        // Per-phase syncs intentionally removed
        // (cuda-md-review-cleanup-03-launch-purity ac-005):
        //   - Components within a phase already serialize via CUDA's
        //     in-stream FIFO ordering. The previous unconditional
        //     `cudaStreamSynchronize(stream) + cudaStreamSynchronize(defer)`
        //     drained both streams 6× per step (per-phase × main + defer),
        //     ~80 µs/step pure host wait on Hopper for no semantic gain.
        //   - Cross-stream dependencies (defer_stream reads main_stream's
        //     buffers) are explicitly ordered by AsyncBridge::launch via
        //     cudaEventRecord(main) + cudaStreamWaitEvent(defer) per
        //     deferred op, so cross-stream reads see consistent values
        //     without per-phase sync.
        //   - End-of-run sync happens in Backend::sync() / retrieve();
        //     individual components issue event syncs when they need
        //     specific host-visible values (NPT, NeighborList).
        if (timer_)
            timer_->impl()->flush();
    }
    if (timer_)
        timer_->step_done();
}

void Backend::run(int n_steps) {
    for (int s = 0; s < n_steps; ++s)
        step(s);
}

void Backend::destroy_components() {
    DestroyContext ctx{context_, frame_};
    for (int s = static_cast<int>(plan_->slots.size()) - 1; s >= 0; --s) {
        auto& slot = plan_->slots[s];
        for (int o = static_cast<int>(slot.ops.size()) - 1; o >= 0; --o)
            slot.ops[o].component->destroy(ctx);
    }
}

void Backend::retrieve(cpu::Frame& frame) { cuda::to_cpu(frame_, frame); }

void Backend::sync() {
    cudaStreamSynchronize(static_cast<cudaStream_t>(context_.stream));
    cudaStreamSynchronize(static_cast<cudaStream_t>(context_.defer_stream));
}

void Backend::shutdown() {
    cuda::destroy(frame_);
    frame_ = {};
    auto main = static_cast<cudaStream_t>(context_.stream);
    auto defer = static_cast<cudaStream_t>(context_.defer_stream);
    if (main)
        cudaStreamDestroy(main);
    if (defer)
        cudaStreamDestroy(defer);
    context_ = {};
    plan_ = nullptr;
    cuda::unregister_log_callback();
}

} // namespace cuda

// ── Explicit instantiations (CUDA) ───────────────────────────────

template class Driver<md::Phase, cuda::device_tag>;
template class Driver<minimizer::Phase, cuda::device_tag>;
template class Driver<evaluate::Phase, cuda::device_tag>;

} // namespace ATV
