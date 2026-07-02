/// @file cpu/observer/rdf_histogram.h
/// @brief Legacy hand-rolled RDF histogram helpers — CUDA observer only —
/// plus the single g(r) file emitter shared with the CPU observer.
///
/// The CPU `cpu::RDF` observer no longer uses the histogram pair: it streams
/// snapshots through molrs (`cpu::compute::RDFAccumulator` in
/// observer/trajectory_metrics.h). `rdf_accumulate`/`rdf_write` are retained
/// **solely** for `cuda::RDF`, which bins pair distances streamed off the
/// device (D2H) and has no molrs path yet. Kept in its own CUDA-safe header
/// (no molrs / `rust/cxx.h`) so `cuda/observer/rdf.cuh` can include it without
/// pulling the CPU observer's header.
/// TODO(roykid, 2026-07-02): delete `rdf_accumulate`/`rdf_write` (and align
/// the `long long` histogram with `BigInt`) once the CUDA observer migrates
/// to the molrs streaming accumulator.
///
/// `rdf_write_curve` is NOT legacy: it is the one emitter of the g(r) file
/// format ("# r [Angstrom]   g(r)"), used by both the CPU observer (molrs
/// curve) and legacy `rdf_write` (normalized histogram).
#pragma once

#include <string>
#include <vector>

namespace ATV {
namespace cpu {

/// @brief Write a normalized g(r) curve to file — the single emitter of the
/// g(r) file format.
/// @param filename output path. @param gr normalized g(r), one value per bin
/// (missing tail bins are written as 0). @param r_min lower edge of bin 0 [Å].
/// @param bin_width radial bin width [Å]. @param n_bins number of lines to
/// write (bin-center radii `r_min + (i + 0.5)·bin_width`).
void rdf_write_curve(const std::string& filename, const std::vector<double>& gr,
                     double r_min, double bin_width, int n_bins);

/// @brief Bin one snapshot of pair distances into a histogram.
void rdf_accumulate(std::vector<long long>& histogram, const double* h_dist,
                    int n_pairs, double bin_width, int n_bins);

/// @brief Normalize an accumulated histogram and write g(r) to file (via
/// `rdf_write_curve`).
void rdf_write(const std::string& filename,
               const std::vector<long long>& histogram, int n_samples,
               int n_atoms, double box_vol, double r_max, int n_bins);

} // namespace cpu
} // namespace ATV
