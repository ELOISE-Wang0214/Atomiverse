/// @file cpu/observer/rdf.h
/// @brief CPU-side RDF: free-function helpers + async-observer component.
///
/// The free-function pair (`rdf_accumulate`, `rdf_write`) is called by
/// `cuda::RDF` after async D2H of pair distances — the CUDA component owns
/// no histogram logic, it lives here.
///
/// `RDF` is a standalone MD async observer that computes pair distances
/// directly from a `FrameSnapshot` on the analysis worker thread and
/// accumulates them into a shared histogram, finalizing the normalized
/// g(r) to file at end-of-run.
///
/// @see .claude/specs/parallel-md-engine.md  § B.4 (async observer pattern)
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/frame_snapshot.h"
#include "atomiverse/cpu/md/async_observer.h"

namespace ATV {

/// @brief MD async radial-distribution-function observer.
///
/// Accumulates pair distances from positions in the captured `FrameSnapshot`
/// on the analysis worker thread; `finalize()` normalizes and writes g(r).
///
/// Usage:
/// @code{.cpp}
/// cpu::md::MD sim;
/// sim.add(cpu::md::VelocityVerlet{0.001});
/// sim.on(md::Phase::force).add(cpu::md::NeighborList{6.0}, -10);
/// sim.on(md::Phase::force).add(cpu::md::LJ{0.0104, 3.4, 6.8});
///
/// RDF rdf{5.5, 200, 100, "rdf.dat"};  // r_max, bins, freq
/// sim.add(rdf);
/// sim.build();
/// sim.run(frame, 1'000'000);
/// @endcode
class RDF : public cpu::md::AsyncObserver {
public:
    RDF(double r_max, int n_bins, int freq = 0,
        std::string output_file = "rdf.dat");

    std::string name() const override { return "rdf_async"; }
    ComponentRole role() const override { return ComponentRole::Observer; }
    int frequency() const override { return freq_; }

    cpu::CaptureMask capture_mask() const override {
        cpu::CaptureMask m;
        m.positions = true;
        return m;
    }

    void on_snapshot(const cpu::FrameSnapshot& snap) override;
    void finalize() override;

    const std::vector<long long>& histogram() const { return histogram_; }
    int n_samples() const { return n_samples_; }
    const std::string& output_file() const { return output_file_; }

private:
    double r_max_;
    int n_bins_;
    int freq_;
    std::string output_file_;
    double bin_width_;

    std::vector<long long> histogram_;
    std::vector<double>    dist_scratch_;
    int    n_samples_ = 0;
    int    n_atoms_last_ = 0;
    double box_vol_last_ = 0.0;
    double box_vol_sum_ = 0.0;
};

namespace cpu {

/// @brief Bin one snapshot of pair distances into a histogram.
void rdf_accumulate(std::vector<long long>& histogram,
                    const double* h_dist, int n_pairs,
                    double bin_width, int n_bins);

/// @brief Normalize accumulated histogram and write g(r) to file.
void rdf_write(const std::string& filename,
               const std::vector<long long>& histogram,
               int n_samples, int n_atoms,
               double box_vol, double r_max, int n_bins);

} // namespace cpu
} // namespace ATV
