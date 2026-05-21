/// @file cuda/async_bridge.cu
/// @brief AsyncBridge sidecar boilerplate: event creation, defer-stream ordering, sync.

#include "atomiverse/cuda/async_bridge.cuh"

#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

void AsyncBridge::build(const BuildContext& ctx) {
    cudaEvent_t ev;
    cudaEventCreateWithFlags(&ev, cudaEventDisableTiming);
    event_ = ev;
    on_build(ctx);
}

void AsyncBridge::launch(const LaunchContext& ctx) {
    if (!can_run(ctx))
        return;

    auto main_stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto defer_stream = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
    auto event = static_cast<cudaEvent_t>(event_);

    // Order defer stream after all main-stream physics work (non-blocking)
    cudaEventRecord(event, main_stream);
    cudaStreamWaitEvent(defer_stream, event, 0);

    // Subclass issues all cudaMemcpyAsync calls on the defer stream
    issue_d2h(ctx);

    // Snapshot boundary: wait only for the defer stream
    cudaStreamSynchronize(defer_stream);

    // Subclass calls the CPU-side observer / recorder
    on_data_ready(ctx);
}

void AsyncBridge::destroy(const DestroyContext& ctx) {
    on_destroy(ctx);
    if (event_) {
        cudaEventDestroy(static_cast<cudaEvent_t>(event_));
        event_ = nullptr;
    }
}

} // namespace ATV
