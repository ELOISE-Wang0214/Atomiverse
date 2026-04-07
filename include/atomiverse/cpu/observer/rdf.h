/// @file cpu/observer/rdf.h
/// @brief CPU-side RDF accumulation and normalization.
///
/// Called by cuda::RDF after async D2H of pair distances.
/// The CUDA component owns no histogram logic -- it is here.
#pragma once

#include <string>
#include <vector>

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

} // namespace cpu
} // namespace ATV
