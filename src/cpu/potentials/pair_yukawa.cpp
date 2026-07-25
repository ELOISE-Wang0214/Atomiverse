/// @file cpu/potentials/pair_yukawa.cpp
#include "atomiverse/cpu/potentials/pair_yukawa.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"

namespace ATV {
namespace cpu {

PairYukawa::PairYukawa(double A, double kappa, double cutoff)
    : Pair(cutoff), A_(A), kappa_(kappa) {
    if (!(cutoff > 0.0))
        throw std::invalid_argument("PairYukawa: cutoff must be > 0");
    if (!(kappa >= 0.0))
        throw std::invalid_argument("PairYukawa: kappa must be ≥ 0");
    cutoff_sq_ = cutoff * cutoff;
    const double urc = A * std::exp(-kappa * cutoff) / cutoff;
    v_shift_ = urc;
}

void PairYukawa::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = static_cast<int>(f.nlocal());

    double energy = 0.0;
    const StorageReal* __restrict x = f.span(kAtom_x);
    const StorageReal* __restrict y = f.span(kAtom_y);
    const StorageReal* __restrict z = f.span(kAtom_z);
    AccReal* __restrict fx = f.span(kAtom_fx);
    AccReal* __restrict fy = f.span(kAtom_fy);
    AccReal* __restrict fz = f.span(kAtom_fz);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double dz = z[j] - z[i];
            f.simbox.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= cutoff_sq_ || r2 == 0.0)
                continue;

            const double r = std::sqrt(r2);
            const double e = std::exp(-kappa_ * r);
            const double inv_r = 1.0 / r;
            const double v = A_ * e * inv_r - v_shift_;
            // dU/dr (unshifted) = −A e (−κ/r + 1/r²) wait:
            // U = A e^{−κr}/r → dU/dr = A e (−κ/r − 1/r²) = −A e (κ/r + 1/r²)
            const double dUdr = -A_ * e * (kappa_ * inv_r + inv_r * inv_r);
            const double fpair =
                -dUdr / r; // F on i from j: fpair * r_ij with r_ij = r_j-r_i

            fx[i] -= fpair * dx;
            fy[i] -= fpair * dy;
            fz[i] -= fpair * dz;
            fx[j] += fpair * dx;
            fy[j] += fpair * dy;
            fz[j] += fpair * dz;
            energy += v;
        }
    }
    f.meta.get(kEnergy) += energy;
}

} // namespace cpu
} // namespace ATV
