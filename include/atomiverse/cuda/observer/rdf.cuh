/// @file cuda/observer/rdf.cuh
/// @brief RDF -- header-only radial distribution function observer (async D2H bridge).
///
/// Accumulates pair distances into histogram bins over multiple steps.
/// Normalization and file output happen in destroy() via cpu::rdf_write.
///
/// @see cpu::rdf_accumulate, cpu::rdf_write (cpu/observer/rdf.h)
#pragma once

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <cuda_runtime.h>

#include "atomiverse/cpu/observer/rdf.h"
#include "atomiverse/cuda/async_bridge.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/types.cuh"

namespace ATV {
namespace cuda {

/// @brief Radial distribution function g(r) observer.
///
/// @code{.cpp}
/// sim.on(md::Phase::end_of_step).add(RDF{3.0, 100, "rdf.dat"});
/// @endcode
class RDF : public AsyncBridge {
  public:
    RDF(double r_max, int n_bins, std::string filename)
        : r_max_(r_max), n_bins_(n_bins), filename_(std::move(filename)) {}

    ~RDF() override = default;
    RDF(RDF&&) noexcept = default;
    RDF& operator=(RDF&&) noexcept = default;

    std::string name() const override { return "rdf"; }

  protected:
    void on_build(const BuildContext& ctx) override {
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        n_atoms_ = frame.n_atoms();
        n_samples_ = 0;

        // D2H box to compute volume (sync OK at build time).
        double box_data[9];
        cudaMemcpyAsync(box_data, frame.box, 9 * sizeof(double),
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);

        double ax = box_data[0], ay = box_data[1], az = box_data[2];
        double bx = box_data[3], by = box_data[4], bz = box_data[5];
        double cx = box_data[6], cy = box_data[7], cz = box_data[8];
        box_vol_ = std::abs(ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx)
                            + az * (bx * cy - by * cx));

        // D2H pair capacity to size pinned buffer. The pair list is owned by
        // NeighborList and published on the context.
        int cap = 0;
        const PairList* pairs = ctx.runtime.pair_list;
        if (pairs && pairs->capacity) {
            cudaMemcpyAsync(&cap, pairs->capacity, sizeof(int),
                            cudaMemcpyDeviceToHost, stream);
            cudaStreamSynchronize(stream);
        }
        pair_capacity_ = static_cast<int>(cap);

        // pairs->dist is a real_t device buffer -- size by real_t.
        if (pair_capacity_ > 0) {
            cudaMallocHost(&h_dist_, pair_capacity_ * sizeof(real_t));
            dist_.resize(static_cast<size_t>(pair_capacity_));
        }
        cudaMallocHost(&h_n_pairs_, sizeof(int));
        *h_n_pairs_ = int(0);

        histogram_.assign(n_bins_, 0);
    }

    bool can_run(const LaunchContext& ctx) const override {
        const PairList* pairs = ctx.runtime.pair_list;
        return pairs != nullptr && pairs->n_pairs != nullptr
               && pair_capacity_ > 0;
    }

    void issue_d2h(const LaunchContext& ctx) override {
        auto ds = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        const PairList* pairs = ctx.runtime.pair_list;
        cudaMemcpyAsync(h_n_pairs_, pairs->n_pairs, sizeof(int),
                        cudaMemcpyDeviceToHost, ds);
        cudaMemcpyAsync(h_dist_, pairs->dist, pair_capacity_ * sizeof(real_t),
                        cudaMemcpyDeviceToHost, ds);
    }

    void on_data_ready(const LaunchContext&) override {
        int n_pairs = *h_n_pairs_;
        if (n_pairs <= 0)
            return;
        int np = static_cast<int>(n_pairs);
        if (np > pair_capacity_)
            np = pair_capacity_;
        // Widen real_t staging buffer to double for the CPU accumulator.
        for (int i = 0; i < np; ++i)
            dist_[i] = static_cast<double>(h_dist_[i]);
        double bin_width = r_max_ / double(n_bins_);
        cpu::rdf_accumulate(histogram_, dist_.data(), np, bin_width, n_bins_);
        n_samples_++;
    }

    void on_destroy(const DestroyContext&) override {
        cpu::rdf_write(filename_, histogram_, n_samples_, n_atoms_, box_vol_,
                       r_max_, n_bins_);
        if (h_dist_) {
            cudaFreeHost(h_dist_);
            h_dist_ = nullptr;
        }
        if (h_n_pairs_) {
            cudaFreeHost(h_n_pairs_);
            h_n_pairs_ = nullptr;
        }
    }

  private:
    double r_max_;
    int n_bins_;
    std::string filename_;

    int n_atoms_ = 0;
    double box_vol_ = 0;
    int n_samples_ = 0;
    int pair_capacity_ = 0;
    std::vector<long long> histogram_;
    real_t* h_dist_ = nullptr;
    std::vector<double> dist_;
    int* h_n_pairs_ = nullptr;
};

} // namespace cuda
} // namespace ATV
