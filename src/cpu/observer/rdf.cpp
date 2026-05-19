/// @file cpu/observer/rdf.cpp
/// @brief CPU RDF accumulation, file output, and `RDF` async-observer
/// implementation.

#include "atomiverse/cpu/observer/rdf.h"
#include "atomiverse/core/constants.h"

#include <algorithm>
#include <cmath>
#include <fstream>

#include <spdlog/fmt/fmt.h>

namespace ATV {

// ── RDF component ──────────────────────────────────────────────────

RDF::RDF(double r_max, int n_bins, int freq, std::string output_file)
    : r_max_(r_max),
      n_bins_(n_bins),
      freq_(freq),
      output_file_(std::move(output_file)),
      bin_width_(r_max / static_cast<double>(std::max(1, n_bins))),
      histogram_(static_cast<std::size_t>(std::max(1, n_bins)), 0LL)
{}

void RDF::on_snapshot(const cpu::FrameSnapshot& snap) {
    const int n = snap.n_owned;
    if (n < 2) return;

    const double r_max_sq = r_max_ * r_max_;

    const std::size_t worst_pairs =
        static_cast<std::size_t>(n) * static_cast<std::size_t>(n - 1) / 2;
    if (dist_scratch_.size() < worst_pairs) dist_scratch_.resize(worst_pairs);

    int n_pairs = 0;
    for (int i = 0; i < n - 1; ++i) {
        const double xi = snap.x[static_cast<std::size_t>(i)];
        const double yi = snap.y[static_cast<std::size_t>(i)];
        const double zi = snap.z[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < n; ++j) {
            double dx = xi - snap.x[static_cast<std::size_t>(j)];
            double dy = yi - snap.y[static_cast<std::size_t>(j)];
            double dz = zi - snap.z[static_cast<std::size_t>(j)];
            snap.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < r_max_sq) {
                dist_scratch_[static_cast<std::size_t>(n_pairs++)] = std::sqrt(r2);
            }
        }
    }

    cpu::rdf_accumulate(histogram_, dist_scratch_.data(), n_pairs,
                        bin_width_, n_bins_);
    ++n_samples_;
    n_atoms_last_ = n;

    const double vol = snap.box.volume();
    box_vol_sum_ += vol;
    box_vol_last_ = box_vol_sum_ / static_cast<double>(n_samples_);
}

void RDF::finalize() {
    if (n_samples_ == 0 || n_atoms_last_ == 0 || box_vol_last_ <= 0.0) return;
    cpu::rdf_write(output_file_, histogram_, n_samples_, n_atoms_last_,
                   box_vol_last_, r_max_, n_bins_);
}

} // namespace ATV

namespace ATV {
namespace cpu {

void rdf_accumulate(std::vector<long long>& histogram,
                    const double* h_dist, int n_pairs,
                    double bin_width, int n_bins)
{
    for (int k = 0; k < n_pairs; ++k) {
        int bin = static_cast<int>(h_dist[k] / bin_width);
        if (bin >= 0 && bin < n_bins)
            histogram[bin]++;
    }
}

void rdf_write(const std::string& filename,
               const std::vector<long long>& histogram,
               int n_samples, int n_atoms,
               double box_vol, double r_max, int n_bins)
{
    if (n_samples <= 0 || n_atoms <= 0 || box_vol <= 0.0) return;

    const double bin_width = r_max / static_cast<double>(n_bins);
    const double rho = static_cast<double>(n_atoms) / box_vol;
    constexpr double pi = kPI;

    std::ofstream out(filename);
    out << "# r [Angstrom]   g(r)\n";

    for (int i = 0; i < n_bins; ++i) {
        const double r = (static_cast<double>(i) + 0.5) * bin_width;
        const double r_lo = static_cast<double>(i) * bin_width;
        const double r_hi = r_lo + bin_width;
        const double shell_vol = (4.0 / 3.0) * pi *
                         (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
        const double ideal = static_cast<double>(n_atoms) * rho * shell_vol;

        const double gr = (ideal > 0.0)
            ? 2.0 * static_cast<double>(histogram[i]) /
                  (static_cast<double>(n_samples) * ideal)
            : 0.0;

        out << fmt::format("{:.6f}  {:.6f}\n", r, gr);
    }
}

} // namespace cpu
} // namespace ATV
