/// @file cpu/observer/rdf.h
/// @brief CPU-side RDF: free-function helpers + MD async-observer wrapper.
///
/// The free-function pair (`rdf_accumulate`, `rdf_write`) is called by
/// `cuda::RDF` after async D2H of pair distances — the CUDA component owns
/// no histogram logic, it lives here.
///
/// `cpu::observer::RDF` is a standalone MD async observer that computes
/// pair distances directly from a `FrameSnapshot` on the analysis worker
/// thread and accumulates them into a shared histogram, finalizing the
/// normalized g(r) to file at end-of-run.
///
/// @see .claude/specs/parallel-md-engine.md  § B.4 (async observer pattern)
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/frame_snapshot.h"
#include "atomiverse/cpu/md/async_observer.h"

namespace ATV {
namespace cpu {

/// @brief Bin one snapshot of pair distances into a histogram.
///
/// @param histogram   Accumulated count per bin [n_bins] (updated in-place).
/// @param h_dist      Host-resident pair distances [n_pairs].
/// @param n_pairs     Number of valid pairs.
/// @param bin_width   Width of each histogram bin [Angstrom].
/// @param n_bins      Number of bins.
void rdf_accumulate(std::vector<long long>& histogram,
                    const double* h_dist, int n_pairs,
                    double bin_width, int n_bins);

/// @brief Normalize accumulated histogram and write g(r) to file.
///
/// Normalization: g(r) = 2 * hist[bin] / (n_samples * N * rho * 4pi r^2 dr).
/// Factor of 2 because the pair list is a half-list.
///
/// @param filename   Output file path.
/// @param histogram  Accumulated count per bin [n_bins].
/// @param n_samples  Number of snapshots accumulated.
/// @param n_atoms    Total number of atoms.
/// @param box_vol    Simulation box volume [Angstrom^3].
/// @param r_max      Maximum distance [Angstrom].
/// @param n_bins     Number of bins.
void rdf_write(const std::string& filename,
               const std::vector<long long>& histogram,
               int n_samples, int n_atoms,
               double box_vol, double r_max, int n_bins);

namespace observer {

/// @brief CPU MD async radial-distribution-function observer.
///
/// Accumulates pair distances from positions in the captured `FrameSnapshot`
/// on the analysis worker thread; `finalize()` normalizes and writes g(r).
/// Phase 1: orthorhombic minimum-image convention, O(N²) half-list per
/// snapshot. Sufficient for systems up to ~10⁴ atoms with sparse sampling;
/// larger systems should use a neighbor-list variant (Phase 2+).
///
/// Usage:
/// @code{.cpp}
/// cpu::md::MD sim;
/// sim.add(cpu::md::VelocityVerlet{0.001});
/// sim.on(md::Phase::force).add(cpu::md::NeighborList{6.0}, -10);
/// sim.on(md::Phase::force).add(cpu::md::LJ{0.0104, 3.4, 6.8});
///
/// cpu::observer::RDF rdf{5.5, 200, 100, "rdf.dat"};  // r_max, bins, freq
/// sim.add(rdf);
/// sim.build();
/// sim.run(frame, 1'000'000);
/// @endcode
class RDF : public cpu::md::AsyncObserver {
public:
    /// @brief Construct an RDF observer.
    /// @param r_max       Histogram upper bound [Angstrom].
    /// @param n_bins      Number of histogram bins.
    /// @param freq        Firing stride (0 = every step after 0).
    /// @param output_file Path for the normalized g(r) output (written on finalize).
    RDF(double r_max, int n_bins, int freq = 0,
        std::string output_file = "rdf.dat");

    std::string name() const override { return "rdf_async"; }
    int frequency() const override { return freq_; }

    cpu::CaptureMask capture_mask() const override {
        cpu::CaptureMask m;
        m.positions = true;
        return m;
    }

    void on_snapshot(const cpu::FrameSnapshot& snap) override;
    void finalize() override;

    /// @brief Accumulated counts per bin (for tests / diagnostics).
    /// Valid after run completes (finalize has no side-effect on the histogram).
    const std::vector<long long>& histogram() const { return histogram_; }

    /// @brief Number of snapshots accumulated.
    int n_samples() const { return n_samples_; }

    /// @brief Output-file path that will receive the g(r) on finalize.
    const std::string& output_file() const { return output_file_; }

private:
    double r_max_;
    int n_bins_;
    int freq_;
    std::string output_file_;
    double bin_width_;

    std::vector<long long> histogram_;
    std::vector<double>    dist_scratch_;  // reused pair-distance buffer
    int    n_samples_ = 0;
    int    n_atoms_last_ = 0;          // for normalization in finalize()
    double box_vol_last_ = 0.0;        // averaged across NVT samples
    double box_vol_sum_ = 0.0;         // running sum for NPT volume averaging
};

} // namespace observer
} // namespace cpu
} // namespace ATV
