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
#include "atomiverse/cuda/async/async_bridge.cuh"
#include "atomiverse/cuda/neighbor/neighbor_list.cuh"
#include "atomiverse/cuda/frame/types.cuh"

namespace ATV {
namespace cuda {

/// @brief Radial distribution function g(r) observer.
///
/// Freezes the per-step pair-distance buffer into a device-ring slot on the
/// compute stream (after the neighbor build that produced it, before the next
/// step rebuilds it), so the histogram accumulates a self-consistent step-N
/// snapshot rather than a buffer racing the next neighbor rebuild.
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
    /// @brief Observer — accumulates pair-distance histogram (no I/O per step).
    ComponentRole role() const override { return ComponentRole::Observer; }

  protected:
    void on_build(const BuildContext& ctx) override {
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        n_atoms_ = frame.n_atoms();
        n_samples_ = 0;

        // RDF histograms pair distances from the flat PairList, which the LJ
        // force path no longer needs — request it so NeighborList re-enables
        // the per-step export. (All build()s run before any launch(), so this
        // is visible to NeighborList::launch from the first step.)
        ctx.runtime.flat_pairs_requested = true;

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

        // D2H pair capacity to size device-ring + pinned buffers. The pair list
        // is owned by NeighborList and published on the context.
        int cap = 0;
        const PairList* pairs = ctx.runtime.pair_list;
        if (pairs && pairs->capacity) {
            cudaMemcpyAsync(&cap, pairs->capacity, sizeof(int),
                            cudaMemcpyDeviceToHost, stream);
            cudaStreamSynchronize(stream);
        }
        pair_capacity_ = static_cast<int>(cap);

        // pairs->dist_sq is a real_t device buffer -- size by real_t.
        // (Per cuda-md-review-cleanup-06-launch-config-streams the pair list
        // stores r² rather than r; this observer takes the sqrt at read time.)
        if (pair_capacity_ > 0) {
            for (int i = 0; i < kRingDepth; ++i) {
                cudaMalloc(&d_dsq_[i], pair_capacity_ * sizeof(real_t));
                cudaMalloc(&d_np_[i], sizeof(int));
                cudaMallocHost(&h_dist_sq_[i], pair_capacity_ * sizeof(real_t));
                cudaMallocHost(&h_n_pairs_[i], sizeof(int));
                *h_n_pairs_[i] = int(0);
            }
            dist_.resize(static_cast<size_t>(pair_capacity_));
        }

        histogram_.assign(n_bins_, 0);
    }

    bool can_run(const LaunchContext& ctx) const override {
        const PairList* pairs = ctx.runtime.pair_list;
        return pairs != nullptr && pairs->n_pairs != nullptr
               && pair_capacity_ > 0;
    }

    void freeze_slot(const LaunchContext& ctx, int slot) override {
        // D2D freeze of the step-N pair buffer on the compute stream, before
        // the next step's neighbor rebuild overwrites it.
        auto ms = static_cast<cudaStream_t>(ctx.runtime.stream);
        const PairList* pairs = ctx.runtime.pair_list;
        cudaMemcpyAsync(d_np_[slot], pairs->n_pairs, sizeof(int),
                        cudaMemcpyDeviceToDevice, ms);
        cudaMemcpyAsync(d_dsq_[slot], pairs->dist_sq,
                        pair_capacity_ * sizeof(real_t),
                        cudaMemcpyDeviceToDevice, ms);
    }

    void d2h_slot(const LaunchContext& ctx, int slot) override {
        auto ds = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        cudaMemcpyAsync(h_n_pairs_[slot], d_np_[slot], sizeof(int),
                        cudaMemcpyDeviceToHost, ds);
        cudaMemcpyAsync(h_dist_sq_[slot], d_dsq_[slot],
                        pair_capacity_ * sizeof(real_t),
                        cudaMemcpyDeviceToHost, ds);
    }

    void on_data_ready(const LaunchContext&, int slot) override {
        int n_pairs = *h_n_pairs_[slot];
        if (n_pairs <= 0)
            return;
        int np = static_cast<int>(n_pairs);
        if (np > pair_capacity_)
            np = pair_capacity_;
        // Widen real_t staging buffer to double for the CPU accumulator and
        // take the sqrt that the export kernel skipped (per
        // cuda-md-review-cleanup-06-launch-config-streams).
        for (int i = 0; i < np; ++i)
            dist_[i] = std::sqrt(static_cast<double>(h_dist_sq_[slot][i]));
        double bin_width = r_max_ / double(n_bins_);
        cpu::rdf_accumulate(histogram_, dist_.data(), np, bin_width, n_bins_);
        n_samples_++;
    }

    void on_destroy(const DestroyContext&) override {
        cpu::rdf_write(filename_, histogram_, n_samples_, n_atoms_, box_vol_,
                       r_max_, n_bins_);
        for (int i = 0; i < kRingDepth; ++i) {
            if (d_dsq_[i]) { cudaFree(d_dsq_[i]); d_dsq_[i] = nullptr; }
            if (d_np_[i]) { cudaFree(d_np_[i]); d_np_[i] = nullptr; }
            if (h_dist_sq_[i]) { cudaFreeHost(h_dist_sq_[i]); h_dist_sq_[i] = nullptr; }
            if (h_n_pairs_[i]) { cudaFreeHost(h_n_pairs_[i]); h_n_pairs_[i] = nullptr; }
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
    real_t* d_dsq_[kRingDepth] = {}; ///< device-ring frozen r² buffer
    int* d_np_[kRingDepth] = {};     ///< device-ring frozen pair count
    real_t* h_dist_sq_[kRingDepth] = {}; ///< pinned host r² (sqrt at read time)
    int* h_n_pairs_[kRingDepth] = {};    ///< pinned host pair count
    std::vector<double> dist_;       ///< widened-to-double sqrt(r²) for histogram
};

} // namespace cuda
} // namespace ATV
