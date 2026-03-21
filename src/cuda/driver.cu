/// @file driver.cu
/// @brief CUDA Backend implementation + Driver explicit instantiations.

#include "atomiverse/cuda/backend.cuh"
#include "atomiverse/cuda/backend_traits.cuh"
#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/facade/driver_impl.h"
#include "atomiverse/facade/stage.h"
#include "atomiverse/cuda/timer.cuh"

namespace ATV {

// ── Backend ────────────────────────────────────────────────────────

void Backend::init() {
    cudaStream_t main, defer;
    cudaStreamCreate(&main);
    cudaStreamCreate(&defer);
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
            if (timer_) timer_->begin(static_cast<cudaStream_t>(context_.stream));
            op.component->build(ctx);
            if (timer_) {
                timer_->end(static_cast<cudaStream_t>(context_.stream),
                            slot.phase, op.component->name() + "/build");
            }
        }
    cudaStreamSynchronize(static_cast<cudaStream_t>(context_.stream));
    if (timer_) timer_->flush();
}

void Backend::step(int step_idx) {
    for (auto& slot : plan_->slots) {
        if (slot.ops.empty()) continue;
        for (auto& op : slot.ops) {
            if (op.frequency > 0 && step_idx % op.frequency != 0)
                continue;
            auto stream = static_cast<cudaStream_t>(
                (op.scheduling == Scheduling::deferrable)
                    ? context_.defer_stream : context_.stream);
            LaunchContext ctx{context_, frame_, slot.phase, step_idx};
            if (timer_) timer_->begin(stream);
            op.component->launch(ctx);
            if (timer_) timer_->end(stream, slot.phase, op.component->name());
        }
        cudaStreamSynchronize(static_cast<cudaStream_t>(context_.stream));
        cudaStreamSynchronize(static_cast<cudaStream_t>(context_.defer_stream));
        if (timer_) timer_->flush();
    }
    if (timer_) timer_->step_done();
}

void Backend::run(int n_steps) {
    for (int s = 0; s < n_steps; ++s) step(s);
}

void Backend::destroy_components() {
    DestroyContext ctx{context_, frame_};
    for (int s = static_cast<int>(plan_->slots.size()) - 1; s >= 0; --s) {
        auto& slot = plan_->slots[s];
        for (int o = static_cast<int>(slot.ops.size()) - 1; o >= 0; --o)
            slot.ops[o].component->destroy(ctx);
    }
}

void Backend::retrieve(cpu::Frame& frame) {
    cuda::to_cpu(frame_, frame);
}

void Backend::sync() {
    cudaStreamSynchronize(static_cast<cudaStream_t>(context_.stream));
    cudaStreamSynchronize(static_cast<cudaStream_t>(context_.defer_stream));
}

void Backend::shutdown() {
    cuda::destroy(frame_);
    frame_ = {};
    auto main = static_cast<cudaStream_t>(context_.stream);
    auto defer = static_cast<cudaStream_t>(context_.defer_stream);
    if (main) cudaStreamDestroy(main);
    if (defer) cudaStreamDestroy(defer);
    context_ = {};
    plan_ = nullptr;
}

// ── Explicit instantiations (CUDA) ───────────────────────────────

template class Driver<md::Phase, cuda::device_tag>;
template class Driver<minimizer::Phase, cuda::device_tag>;
template class Driver<scf::Phase, cuda::device_tag>;

} // namespace ATV
