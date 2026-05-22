/// @file cuda/recorder/density_recorder.cuh
/// @brief DensityRecorder -- header-only electron density grid recorder (async D2H bridge).
///
/// @see cpu::eval_density_grid (cpu/observer/density.h)
/// @see cpu::zarr_append_density (cpu/recorder/zarr.h)
#pragma once

#include <string>
#include <utility>
#include <vector>

#include <cuda_runtime.h>

#include "atomiverse/cpu/observer/density.h"
#include "atomiverse/cuda/async_bridge.cuh"

#ifdef ATV_HAS_ZARR
#include "atomiverse/cpu/recorder/zarr.h"
#endif

#ifdef ATV_HAS_RHF

namespace ATV {
namespace cuda {

/// @brief Evaluate and store electron density rho(r) on a uniform grid.
///
/// @code{.cpp}
/// solver.on(scf::Phase::end_of_step).add(
///     DensityRecorder{"density.zarr", 20, 20, 20});
/// @endcode
class DensityRecorder : public AsyncBridge {
  public:
    DensityRecorder(std::string store_path, int nx, int ny, int nz,
                    double padding = double(2.0), int every = 1)
        : store_path_(std::move(store_path)), nx_(nx), ny_(ny), nz_(nz),
          padding_(padding), every_(every) {}

    ~DensityRecorder() override = default;
    DensityRecorder(DensityRecorder&&) noexcept = default;
    DensityRecorder& operator=(DensityRecorder&&) noexcept = default;

    std::string name() const override { return "density_recorder"; }
    int frequency() const override { return every_; }

  protected:
    void on_build(const BuildContext& ctx) override {
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        n_atoms_ = frame.n_atoms();
        n_basis_ = frame.scf.n;

        const auto na = static_cast<size_t>(n_atoms_);

        htype_.resize(na);
        hx_.resize(na);
        hy_.resize(na);
        hz_.resize(na);
        cudaMemcpyAsync(htype_.data(), frame.atoms.type, na * sizeof(int),
                        cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(hx_.data(), frame.atoms.x, na * sizeof(double),
                        cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(hy_.data(), frame.atoms.y, na * sizeof(double),
                        cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(hz_.data(), frame.atoms.z, na * sizeof(double),
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);

        const auto nb = static_cast<size_t>(n_basis_);
        cudaMallocHost(&h_P_, nb * nb * sizeof(double));
    }

    bool can_run(const LaunchContext& ctx) const override {
        return ctx.frame.scf.P != nullptr && n_basis_ > 0;
    }

    void issue_d2h(const LaunchContext& ctx) override {
        auto ds = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        const size_t nb = static_cast<size_t>(n_basis_);
        cudaMemcpyAsync(h_P_, ctx.frame.scf.P, nb * nb * sizeof(double),
                        cudaMemcpyDeviceToHost, ds);
    }

    void on_data_ready(const LaunchContext& ctx) override {
        auto grid = cpu::eval_density_grid(htype_, hx_, hy_, hz_, h_P_,
                                           n_basis_, nx_, ny_, nz_, padding_);
        if (grid.empty())
            return;
#ifdef ATV_HAS_ZARR
        cpu::zarr_append_density(store_path_, grid, ctx.step);
#endif
    }

    void on_destroy(const DestroyContext&) override {
        if (h_P_) {
            cudaFreeHost(h_P_);
            h_P_ = nullptr;
        }
    }

  private:
    std::string store_path_;
    int nx_, ny_, nz_, every_;
    double padding_;

    int n_atoms_ = 0;
    int n_basis_ = 0;

    double* h_P_ = nullptr;
    std::vector<int> htype_;
    std::vector<double> hx_, hy_, hz_;
};

} // namespace cuda
} // namespace ATV

#endif // ATV_HAS_RHF
