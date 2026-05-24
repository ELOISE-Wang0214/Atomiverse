/// @file cuda/async_bridge.cuh
/// @brief AsyncBridge -- base class for all CUDA async-D2H observer/recorder components.
///
/// Encapsulates the sidecar pattern shared by every CUDA observer and recorder:
///
/// @code
///   build():   cudaEventCreate + on_build(ctx)
///   launch():  EventRecord(main) → StreamWaitEvent(defer) → issue_d2h(ctx)
///              → cudaLaunchHostFunc(defer, &cuda_callback_)
///              → [worker thread: on_data_ready(ctx) → free payload]
///   destroy(): cudaStreamSynchronize(defer)  [drain pending callbacks]
///              → on_destroy(ctx) → cudaEventDestroy
/// @endcode
///
/// Subclasses implement four hooks:
///   - on_build()       allocate pinned buffers, D2H build-time data, open files
///   - issue_d2h()      cudaMemcpyAsync into pinned buffers on the defer stream
///   - on_data_ready()  call the CPU-side observer or recorder function
///   - on_destroy()     free pinned buffers, close files
///
/// The CUDA event and the sidecar stream-ordering logic are handled entirely
/// by this base class.  Subclasses never touch cudaEventCreate/Destroy.
///
/// CPU-side logic MUST live in cpu/observer/ or cpu/recorder/ — not in subclasses.
///
/// @see cpu/observer/rdf.h, cpu/recorder/thermo.h, cpu/recorder/trajectory.h
#pragma once

#include "atomiverse/cuda/component.cuh"

namespace ATV {
namespace cuda {

/// @brief Base class for CUDA async-D2H observer and recorder components.
///
/// All observers and recorders that follow the defer-stream sidecar pattern
/// must extend this class instead of ComponentBase.
class AsyncBridge : public ComponentBase {
  public:
    /// @brief All AsyncBridge subclasses are deferrable by definition.
    Scheduling scheduling() const override { return Scheduling::deferrable; }

    /// @brief Create CUDA event; call on_build().
    /// @param ctx Build context.
    void build(const BuildContext& ctx) final;

    /// @brief Sidecar boilerplate: event-order defer stream, issue
    ///        D2H, dispatch on_data_ready via cudaLaunchHostFunc.
    ///
    /// Calls can_run(ctx) first; returns immediately if false. The
    /// host hook runs on a CUDA worker thread after the D2H completes;
    /// the main thread is free to queue subsequent steps in parallel.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) final;

    /// @brief Call on_destroy(); destroy CUDA event.
    /// @param ctx Destroy context.
    void destroy(const DestroyContext& ctx) final;

  protected:
    /// @brief Allocate pinned host buffers; D2H build-time data; open files.
    ///
    /// The CUDA event is already created before this is called.
    /// @param ctx Build context.
    virtual void on_build(const BuildContext& ctx) = 0;

    /// @brief Issue cudaMemcpyAsync calls on the defer stream.
    ///
    /// The defer stream has already been ordered after the main stream via
    /// cudaStreamWaitEvent before this is called.  Subclasses must use
    /// ctx.runtime.defer_stream for all D2H transfers.
    ///
    /// @param ctx Launch context.
    virtual void issue_d2h(const LaunchContext& ctx) = 0;

    /// @brief CPU-side implementation invoked after the D2H synchronization.
    ///
    /// At this point all pinned buffers contain valid data.  Subclasses call
    /// the corresponding cpu/observer/ or cpu/recorder/ function.
    ///
    /// @param ctx Launch context.
    virtual void on_data_ready(const LaunchContext& ctx) = 0;

    /// @brief Free pinned buffers; close files; release component resources.
    ///
    /// Called before the CUDA event is destroyed.
    /// @param ctx Destroy context.
    virtual void on_destroy(const DestroyContext& ctx) = 0;

    /// @brief Guard hook — return false to skip one launch iteration entirely.
    ///
    /// Default: always run.  Override when the component requires certain frame
    /// data to be present (e.g. pairs.n_pairs, scf.P).
    ///
    /// @param ctx Launch context.
    /// @return true if launch should proceed.
    virtual bool can_run(const LaunchContext& ctx) const { return true; }

    void* event_ =
        nullptr; ///< CUDA event (cudaEvent_t, as void*) owned by this base.

  public:
    /// @brief Stream-callback dispatcher — internal, do not call directly.
    ///
    /// Routed via `cudaLaunchHostFunc` on the defer stream after the
    /// subclass's `issue_d2h(ctx)` completes. Runs on a CUDA worker
    /// thread (not the main thread). Reads the heap-allocated callback
    /// args struct (owned per-call), constructs a `LaunchContext`, and
    /// invokes the subclass's `on_data_ready()` hook. Frees the args.
    ///
    /// Exposed `public` so the C-callable callback signature
    /// (`cudaHostFn_t`) can reach the static dispatch without friend
    /// declarations. Not part of the user-extension API.
    static void cuda_callback_(void* user_data);
};

} // namespace cuda
} // namespace ATV
