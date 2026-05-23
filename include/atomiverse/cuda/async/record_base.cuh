/// @file cuda/async/record_base.cuh
/// @brief `cuda::AsyncRecordBase<Snapshot>` -- zero-block base class for
///        CUDA-side recorder/thermo/observer components.
///
/// Provides the **non-blocking** sidecar pattern that every CUDA-side
/// recorder must follow: the GPU stages a snapshot into a device-resident
/// buffer on the compute stream, the base class enqueues a
/// `cudaMemcpyAsync` D2H copy on the **defer stream** (NOT the compute
/// stream), records a `cudaEvent`, and submits a closure to an in-process
/// host worker thread that calls `cudaEventSynchronize` off the hot path
/// and then invokes the user-supplied `work_fn(host_snapshot, step)`.
///
/// The hot-path contract is:
///   * `submit_async(...)` does NOT call `cudaStreamSynchronize` on the
///     caller's thread (the MD step loop).
///   * The compute stream is never blocked by this base class.
///   * `cudaEventSynchronize` runs on the worker thread, not the hot path.
///
/// This contract is verified by `tests/cuda/async/test_hot_path_no_sync.cu`
/// via the opt-in counter in `cuda/async/sync_guard.cuh`.
///
/// The host worker is owned by `AsyncRecordBase` itself — a minimum-viable
/// single-thread worker that processes a queue of `std::function<void()>`.
/// `destroy()` drains the queue and joins the worker before freeing the
/// device snapshot buffer and the event, so the recorder's CSV is complete
/// when the simulation returns.
///
/// Concrete subclasses (e.g. the CUDA-side `CudaCsvRecorder` in
/// `benchmarks/md/polyethylene_tension.cpp`) launch a tiny pack-snapshot
/// kernel on the compute stream that fills `device_snapshot_`, then call
/// `submit_async(step, [](const Snapshot& s, int64_t step) { ... })`.
///
/// @see cuda/async/sync_guard.cuh for the sync-counter wrapper.
/// @see .claude/specs/msu-pe-cuda-components-07-async-correctness.md
#pragma once

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>

#include <cuda_runtime.h>

#include "atomiverse/cuda/component.cuh"

namespace ATV {
namespace cuda {

namespace detail {

/// @brief Minimum-viable single-thread worker that processes a FIFO queue
///        of `std::function<void()>` closures.
///
/// Owned by `AsyncRecordBase` (one worker per recorder instance). Lifetime:
///   * Constructed by `AsyncRecordBase::build()` -- spawns the worker thread.
///   * `submit(fn)` pushes a closure onto the queue; worker dequeues and runs.
///   * `join()` is called by `AsyncRecordBase::destroy()` to drain the queue
///     and stop the thread cleanly.
///
/// Defined in `src/cuda/async/record_base.cu` so the std::thread / mutex /
/// condition_variable machinery lives in a single translation unit.
class HostWorker {
  public:
    HostWorker();
    ~HostWorker();

    HostWorker(const HostWorker&) = delete;
    HostWorker& operator=(const HostWorker&) = delete;
    HostWorker(HostWorker&&) = delete;
    HostWorker& operator=(HostWorker&&) = delete;

    /// @brief Lazily spawn the worker thread on the first call.
    void start();

    /// @brief Push a closure onto the worker queue.
    /// @param fn Closure invoked on the worker thread (must be self-contained).
    void submit(std::function<void()> fn);

    /// @brief Drain the queue and join the worker thread. Idempotent.
    void join();

  private:
    void loop();

    std::thread thread_;
    std::deque<std::function<void()>> queue_;
    std::mutex mtx_;
    std::condition_variable cv_;
    std::atomic<bool> stop_{false};
    bool started_{false};
    bool joined_{false};
};

} // namespace detail

/// @brief Driver-agnostic base for CUDA recorders/observers that need
///        host-thread access to per-step device state.
///
/// @tparam Snapshot Plain-POD struct holding the fields the recorder
///         needs on the host. Must be trivially copyable so a single
///         `cudaMemcpyAsync` can stage it from device to host.
///
/// **Zero-block contract:** `submit_async(...)` returns to the caller
/// without calling `cudaStreamSynchronize`. The D2H copy uses the defer
/// stream (`ctx.runtime.defer_stream` resolved at `build()` time), so the
/// compute stream is never blocked.
///
/// **Lifecycle:**
///   * `build(ctx)` -- allocates `device_snapshot_` via `cudaMalloc` (owning),
///     creates `copy_done_event_` via `cudaEventCreateWithFlags(disable_timing)`
///     (owning), captures `defer_stream_` from `ctx.runtime.defer_stream`
///     (non-owning), starts the host worker.
///   * `submit_async(step, work_fn)` -- enqueue D2H + event record + worker job.
///   * `destroy(ctx)` -- joins the worker, frees `device_snapshot_`, destroys
///     the event, resets pointer/handle members to nullptr.
template <typename Snapshot>
class AsyncRecordBase : public ComponentBase {
  public:
    AsyncRecordBase() = default;

    /// @brief Allocate device snapshot buffer (owning) + cudaEvents (owning);
    ///        start host worker.
    /// @param ctx Build context — supplies the compute + defer stream
    ///        pointers (non-owning handles borrowed from `cuda::Backend`).
    ///
    /// Two events are owned by this base:
    ///   * `pack_done_event_` — recorded on the compute stream right
    ///     before the D2H copy is enqueued on the defer stream; lets
    ///     the defer stream `cudaStreamWaitEvent` on it to guarantee
    ///     the snapshot was fully staged before the D2H reads it.
    ///   * `copy_done_event_` — recorded on the defer stream after the
    ///     D2H; the worker thread `cudaEventSynchronize`s on this one.
    void build(const BuildContext& ctx) override {
        if (device_snapshot_ != nullptr || copy_done_event_ != nullptr) {
            throw std::runtime_error(
                "cuda::AsyncRecordBase::build called twice without destroy");
        }
        cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&device_snapshot_),
                                     sizeof(Snapshot));
        if (err != cudaSuccess) {
            throw std::runtime_error(
                std::string("cuda::AsyncRecordBase: cudaMalloc failed: ")
                + cudaGetErrorString(err));
        }
        err = cudaEventCreateWithFlags(&copy_done_event_,
                                       cudaEventDisableTiming);
        if (err != cudaSuccess) {
            cudaFree(device_snapshot_);
            device_snapshot_ = nullptr;
            throw std::runtime_error(
                std::string("cuda::AsyncRecordBase: cudaEventCreate "
                            "(copy_done) failed: ")
                + cudaGetErrorString(err));
        }
        err = cudaEventCreateWithFlags(&pack_done_event_,
                                       cudaEventDisableTiming);
        if (err != cudaSuccess) {
            cudaFree(device_snapshot_);
            device_snapshot_ = nullptr;
            cudaEventDestroy(copy_done_event_);
            copy_done_event_ = nullptr;
            throw std::runtime_error(
                std::string("cuda::AsyncRecordBase: cudaEventCreate "
                            "(pack_done) failed: ")
                + cudaGetErrorString(err));
        }
        compute_stream_ = static_cast<cudaStream_t>(ctx.runtime.stream);
        defer_stream_ = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        worker_.start();
    }

    /// @brief Tear down: drain worker, free device buffer, destroy events.
    /// @param ctx Destroy context (unused; signature mandated by base).
    void destroy(const DestroyContext& /*ctx*/) override {
        worker_.join();
        if (device_snapshot_) {
            cudaFree(device_snapshot_);
            device_snapshot_ = nullptr;
        }
        if (copy_done_event_) {
            cudaEventDestroy(copy_done_event_);
            copy_done_event_ = nullptr;
        }
        if (pack_done_event_) {
            cudaEventDestroy(pack_done_event_);
            pack_done_event_ = nullptr;
        }
        compute_stream_ = nullptr;
        defer_stream_ = nullptr;
    }

    /// @brief Convenience: copy @p n_bytes from a device-resident source
    ///        into the snapshot buffer on the compute stream. Subclasses
    ///        call this once per snapshot field they want staged, then
    ///        call `submit_async`. Wraps `cudaMemcpyAsync` D2D so the
    ///        subclass never names the CUDA runtime API directly.
    /// @param dst Pointer inside `device_snapshot_` to write into.
    /// @param src Device-resident source.
    /// @param n_bytes Byte count.
    /// @param compute_stream Stream on which the simulation's prior step
    ///        kernels ran. Pass `ctx.runtime.stream` from `launch`.
    void stage_device_field(void* dst, const void* src, std::size_t n_bytes,
                            cudaStream_t compute_stream) {
        cudaMemcpyAsync(dst, src, n_bytes, cudaMemcpyDeviceToDevice,
                        compute_stream);
    }

    /// @brief Hot-path entry. Enqueues a non-blocking D2H of
    ///        `device_snapshot_` onto the defer stream, records the
    ///        completion event, and submits a closure to the host worker
    ///        that synchronizes on the event and invokes @p work_fn.
    /// @param step Logical simulation step index (passed verbatim to @p work_fn).
    /// @param work_fn Host-side analysis/recording function. Runs on the
    ///        worker thread, never on the hot path.
    void submit_async(
        std::int64_t step,
        std::function<void(const Snapshot&, std::int64_t)> work_fn) {
        // 1. Order the defer stream after the snapshot-pack on the
        //    compute stream so the D2H sees fully-staged bytes. Record
        //    a pack-done event on the compute stream; have the defer
        //    stream wait for it without blocking the host.
        cudaEventRecord(pack_done_event_, compute_stream_);
        cudaStreamWaitEvent(defer_stream_, pack_done_event_, 0);
        // 2. Enqueue D2H of the device snapshot onto the defer stream.
        //    `&host_snapshot_` is a member -- safe as long as `*this`
        //    outlives the worker job, which is guaranteed by
        //    `destroy()` joining the worker before tearing down.
        cudaMemcpyAsync(&host_snapshot_, device_snapshot_, sizeof(Snapshot),
                        cudaMemcpyDeviceToHost, defer_stream_);
        // 3. Record completion event on the defer stream.
        cudaEventRecord(copy_done_event_, defer_stream_);
        // 4. Submit closure to worker. Capture by value the small POD fields;
        //    the snapshot reference points to a member kept alive by step 2's
        //    join-before-free ordering in destroy().
        cudaEvent_t ev = copy_done_event_;
        Snapshot* snap_ptr = &host_snapshot_;
        worker_.submit([ev, snap_ptr, step, work_fn = std::move(work_fn)]() {
            cudaEventSynchronize(ev);
            work_fn(*snap_ptr, step);
        });
        // submit_async returns immediately: no cudaStreamSynchronize on the
        // caller's thread. This is the load-bearing invariant verified by
        // tests/cuda/async/test_hot_path_no_sync.cu.
    }

  protected:
    /// @brief Device-resident staging buffer (owning, freed in `destroy()`).
    ///        Subclasses' `launch()` kernels write into it on the compute
    ///        stream before calling `submit_async`.
    Snapshot* device_snapshot_ = nullptr;

    /// @brief Host-resident snapshot mirror filled by the async D2H.
    ///        Read on the worker thread (after `cudaEventSynchronize`),
    ///        never on the hot path.
    Snapshot host_snapshot_ = {};

    /// @brief cudaEvent signaling D2H completion (owning, destroyed in
    ///        `destroy()`). Recorded on the defer stream after the
    ///        async D2H copy; worker thread polls it.
    cudaEvent_t copy_done_event_ = nullptr;

    /// @brief cudaEvent signaling that the snapshot-pack on the compute
    ///        stream has finished (owning, destroyed in `destroy()`).
    ///        The defer stream `cudaStreamWaitEvent`s on it before
    ///        enqueuing the D2H so the host snapshot sees fully-staged
    ///        bytes.
    cudaEvent_t pack_done_event_ = nullptr;

    /// @brief Compute-stream handle resolved from `ctx.runtime.stream`
    ///        at `build()`. **Non-owning** — `cuda::Backend` owns the stream.
    cudaStream_t compute_stream_ = nullptr;

    /// @brief Defer-stream handle resolved from `ctx.runtime.defer_stream`
    ///        at `build()`. **Non-owning** — `cuda::Backend` owns the stream.
    cudaStream_t defer_stream_ = nullptr;

  private:
    detail::HostWorker worker_;
};

} // namespace cuda
} // namespace ATV
