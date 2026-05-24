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
#include "atomiverse/cuda/types.cuh"

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
    /// @brief Recorder — writes XYZ trajectory frames.
    ComponentRole role() const override { return ComponentRole::Recorder; }
    int frequency() const override { return every_; }

  protected:
    void on_build(const BuildContext& ctx) override {
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        n_atoms_ = frame.n_atoms();
        auto n = static_cast<size_t>(n_atoms_);

        // frame.atoms.x/y/z are real_t device buffers -- size by real_t.
        cudaMallocHost(&hx_, n * sizeof(real_t));
        cudaMallocHost(&hy_, n * sizeof(real_t));
        cudaMallocHost(&hz_, n * sizeof(real_t));

        dx_.resize(n);
        dy_.resize(n);
        dz_.resize(n);

        htype_.resize(n);
        cudaMemcpyAsync(htype_.data(), frame.atoms.type, n * sizeof(int),
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);

        file_.open(filename_);
    }

    void issue_d2h(const LaunchContext& ctx) override {
        auto ds = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        auto n = static_cast<size_t>(n_atoms_);
        cudaMemcpyAsync(hx_, ctx.frame.atoms.x, n * sizeof(real_t),
                        cudaMemcpyDeviceToHost, ds);
        cudaMemcpyAsync(hy_, ctx.frame.atoms.y, n * sizeof(real_t),
                        cudaMemcpyDeviceToHost, ds);
        cudaMemcpyAsync(hz_, ctx.frame.atoms.z, n * sizeof(real_t),
                        cudaMemcpyDeviceToHost, ds);
    }

    void on_data_ready(const LaunchContext& ctx) override {
        // Widen real_t staging buffers to double for the CPU writer.
        for (int i = 0; i < n_atoms_; ++i) {
            dx_[i] = static_cast<double>(hx_[i]);
            dy_[i] = static_cast<double>(hy_[i]);
            dz_[i] = static_cast<double>(hz_[i]);
        }
        cpu::trajectory_append(file_, htype_.data(), dx_.data(), dy_.data(),
                               dz_.data(), n_atoms_, ctx.step);
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
    real_t* hx_ = nullptr;
    real_t* hy_ = nullptr;
    real_t* hz_ = nullptr;
    std::vector<double> dx_, dy_, dz_;
    std::vector<int> htype_;
    std::ofstream file_;
};

} // namespace cuda
} // namespace ATV
