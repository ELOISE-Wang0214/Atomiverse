/// @file cpu/observer/rdf.cpp
/// @brief CPU RDF accumulation and file output.

#include "atomiverse/cpu/observer/rdf.h"
#include "atomiverse/core/constants.h"

#include <cmath>
#include <fstream>

#include <spdlog/fmt/fmt.h>

namespace ATV {
namespace cpu {

void rdf_accumulate(std::vector<long long>& histogram,
                    const Real* h_dist, int n_pairs,
                    Real bin_width, int n_bins)
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
               Real box_vol, Real r_max, int n_bins)
{
    if (n_samples <= 0 || n_atoms <= 0 || box_vol <= Real(0)) return;

    Real bin_width = r_max / Real(n_bins);
    Real rho = Real(n_atoms) / box_vol;
    constexpr Real pi = Real(kPI);

    std::ofstream out(filename);
    out << "# r [Angstrom]   g(r)\n";

    for (int i = 0; i < n_bins; ++i) {
        Real r    = (Real(i) + Real(0.5)) * bin_width;
        Real r_lo = Real(i) * bin_width;
        Real r_hi = r_lo + bin_width;
        Real shell_vol = (Real(4) / Real(3)) * pi *
                         (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
        Real ideal = Real(n_atoms) * rho * shell_vol;

        Real gr = (ideal > Real(0))
            ? Real(2) * Real(histogram[i]) / (Real(n_samples) * ideal)
            : Real(0);

        out << fmt::format("{:.6f}  {:.6f}\n",
                           static_cast<double>(r), static_cast<double>(gr));
    }
}

} // namespace cpu
} // namespace ATV
