/// @file cpu/observer/rdf_histogram.cpp
/// @brief Legacy hand-rolled RDF histogram helpers (CUDA observer only) + the
/// shared g(r) file emitter.
///
/// See rdf_histogram.h: `rdf_accumulate`/`rdf_write` are retained solely for
/// `cuda::RDF` (TODO(roykid, 2026-07-02): delete once the CUDA observer
/// migrates to the molrs streaming accumulator); `rdf_write_curve` is the one
/// emitter of the g(r) file format for both backends.

#include "atomiverse/cpu/observer/rdf_histogram.h"
#include "atomiverse/core/units.h"

#include <cstddef>
#include <fstream>

#include <spdlog/fmt/fmt.h>

namespace ATV {
namespace cpu {

void rdf_write_curve(const std::string& filename, const std::vector<double>& gr,
                     double r_min, double bin_width, int n_bins) {
    std::ofstream out(filename);
    out << "# r [Angstrom]   g(r)\n";
    for (int i = 0; i < n_bins; ++i) {
        const double r = r_min + (static_cast<double>(i) + 0.5) * bin_width;
        const double g = (static_cast<std::size_t>(i) < gr.size())
                             ? gr[static_cast<std::size_t>(i)]
                             : 0.0;
        out << fmt::format("{:.6f}  {:.6f}\n", r, g);
    }
}

void rdf_accumulate(std::vector<long long>& histogram, const double* h_dist,
                    int n_pairs, double bin_width, int n_bins) {
    for (int k = 0; k < n_pairs; ++k) {
        int bin = static_cast<int>(h_dist[k] / bin_width);
        if (bin >= 0 && bin < n_bins)
            histogram[bin]++;
    }
}

void rdf_write(const std::string& filename,
               const std::vector<long long>& histogram, int n_samples,
               int n_atoms, double box_vol, double r_max, int n_bins) {
    if (n_samples <= 0 || n_atoms <= 0 || box_vol <= 0.0)
        return;

    const double bin_width = r_max / static_cast<double>(n_bins);
    const double rho = static_cast<double>(n_atoms) / box_vol;

    std::vector<double> gr(static_cast<std::size_t>(n_bins), 0.0);
    for (int i = 0; i < n_bins; ++i) {
        const double r_lo = static_cast<double>(i) * bin_width;
        const double r_hi = r_lo + bin_width;
        const double shell_vol =
            (4.0 / 3.0) * kPI * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
        const double ideal = static_cast<double>(n_atoms) * rho * shell_vol;
        if (ideal > 0.0)
            gr[static_cast<std::size_t>(i)] =
                2.0
                * static_cast<double>(histogram[static_cast<std::size_t>(i)])
                / (static_cast<double>(n_samples) * ideal);
    }
    rdf_write_curve(filename, gr, /*r_min=*/0.0, bin_width, n_bins);
}

} // namespace cpu
} // namespace ATV
