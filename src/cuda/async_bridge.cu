/// @file cuda/async_bridge.cu
/// @brief AsyncBridge sidecar boilerplate: event ordering + on-device-completion
///        host callback via `cudaLaunchHostFunc`.
///
/// Hot-path contract: `launch()` issues exactly one D2H per subclass
/// (via `issue_d2h(ctx)` on the defer stream), records a host callback
/// to run `on_data_ready(ctx)` when the D2H completes, and returns
/// immediately. No `cudaStreamSynchronize` — the main thread never
/// waits on the defer stream during the per-step loop. The previous
/// `cudaStreamSynchronize(defer_stream)` was the sole AsyncBridge sync
/// blocking the `Backend::step()` zero-sync gate.

#include "atomiverse/cuda/async_bridge.cuh"

#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

namespace ATV {
namespace cuda {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

namespace {

/// @brief Per-call callback payload. Heap-allocated by `launch()` and
///        freed by `cuda_callback_` after the user hook returns.
///
/// `LaunchContext` is non-trivially-copyable (it holds reference
/// members), so we capture its scalar fields and the addresses of its
/// reference targets here. Lifetime: from the `cudaLaunchHostFunc`
/// enqueue in `launch()` until the worker thread invokes
/// `cuda_callback_`. The defer stream's FIFO ordering serializes
/// callbacks, so each instance has at most ~K outstanding callbacks
/// queued (where K = how many steps the GPU is ahead of the host) —
/// the heap traffic is bounded and not on the GPU critical path
/// (~50 ns / step on Hopper).
struct CallbackPayload {
    AsyncBridge* self;
    cuda::Context* runtime;
    cuda::Frame* frame;
    int stage;
    int step;
};

} // anonymous namespace

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

    // Snapshot boundary: enqueue a host callback that runs after the
    // D2H above completes. The callback fires on a CUDA worker thread;
    // the main thread is free to issue the next step's kernels in
    // parallel. This is the link 03 zero-sync replacement for the
    // previous `cudaStreamSynchronize(defer_stream)` + inline
    // `on_data_ready(ctx)` pattern.
    auto* payload = new CallbackPayload{this, &ctx.runtime, &ctx.frame,
                                        ctx.stage, ctx.step};
    cudaLaunchHostFunc(defer_stream, &AsyncBridge::cuda_callback_, payload);
}

void AsyncBridge::destroy(const DestroyContext& ctx) {
    // Drain any in-flight defer-stream callbacks before tearing down
    // pinned buffers — the worker-thread callback would dereference
    // freed memory otherwise. This sync is at teardown only (not on
    // the per-step path), so it does not affect the zero-sync hot-path
    // contract.
    auto defer_stream = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
    cudaStreamSynchronize(defer_stream);

    on_destroy(ctx);
    if (event_) {
        cudaEventDestroy(static_cast<cudaEvent_t>(event_));
        event_ = nullptr;
    }
}

void AsyncBridge::cuda_callback_(void* user_data) {
    // CUDA worker thread. Reconstruct the LaunchContext from the
    // captured pointer/scalar fields and invoke the virtual user hook.
    // Then free the payload — single owner pattern.
    auto* payload = static_cast<CallbackPayload*>(user_data);
    LaunchContext ctx{*payload->runtime, *payload->frame, payload->stage,
                      payload->step};
    payload->self->on_data_ready(ctx);
    delete payload;
}

} // namespace cuda
} // namespace ATV
