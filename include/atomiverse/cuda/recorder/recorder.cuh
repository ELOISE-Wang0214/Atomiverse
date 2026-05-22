/// @file cuda/recorder/recorder.cuh
/// @brief Recorder -- header-only XYZ trajectory writer (async D2H bridge).
///
/// @see cpu::trajectory_append (cpu/recorder/trajectory.h)
#pragma once

#include <cstddef>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <cuda_runtime.h>

#include "atomiverse/cpu/recorder/trajectory.h"
#include "atomiverse/cuda/async_bridge.cuh"

namespace ATV {
namespace cuda {

/// @brief XYZ trajectory recorder.
///
/// @code{.cpp}
/// sim.on(md::Phase::end_of_step).add(Recorder{"traj.xyz", 10});
/// @endcode
class Recorder : public AsyncBridge {
  public:
    Recorder(std::string filename, int every = 1)
        : filename_(std::move(filename)), every_(every) {}

    ~Recorder() override = default;
    Recorder(Recorder&&) noexcept = default;
    Recorder& operator=(Recorder&&) noexcept = default;

    std::string name() const override { return "recorder"; }
    int frequency() const override { return every_; }

  protected:
    void on_build(const BuildContext& ctx) override {
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        n_atoms_ = frame.n_atoms();
        auto n = static_cast<size_t>(n_atoms_);

        cudaMallocHost(&hx_, n * sizeof(double));
        cudaMallocHost(&hy_, n * sizeof(double));
        cudaMallocHost(&hz_, n * sizeof(double));

        htype_.resize(n);
        cudaMemcpyAsync(htype_.data(), frame.atoms.type, n * sizeof(int),
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);

        file_.open(filename_);
    }

    void issue_d2h(const LaunchContext& ctx) override {
        auto ds = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        auto n = static_cast<size_t>(n_atoms_);
        cudaMemcpyAsync(hx_, ctx.frame.atoms.x, n * sizeof(double),
                        cudaMemcpyDeviceToHost, ds);
        cudaMemcpyAsync(hy_, ctx.frame.atoms.y, n * sizeof(double),
                        cudaMemcpyDeviceToHost, ds);
        cudaMemcpyAsync(hz_, ctx.frame.atoms.z, n * sizeof(double),
                        cudaMemcpyDeviceToHost, ds);
    }

    void on_data_ready(const LaunchContext& ctx) override {
        cpu::trajectory_append(file_, htype_.data(), hx_, hy_, hz_, n_atoms_,
                               ctx.step);
    }

    void on_destroy(const DestroyContext&) override {
        if (file_.is_open())
            file_.close();
        if (hx_) {
            cudaFreeHost(hx_);
            hx_ = nullptr;
        }
        if (hy_) {
            cudaFreeHost(hy_);
            hy_ = nullptr;
        }
        if (hz_) {
            cudaFreeHost(hz_);
            hz_ = nullptr;
        }
    }

  private:
    std::string filename_;
    int every_;
    int n_atoms_ = 0;
    double* hx_ = nullptr;
    double* hy_ = nullptr;
    double* hz_ = nullptr;
    std::vector<int> htype_;
    std::ofstream file_;
};

} // namespace cuda
} // namespace ATV
