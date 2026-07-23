/// @file cpu/potentials/pair_born.cpp
#include "atomiverse/cpu/potentials/pair_born.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"

namespace ATV {
namespace cpu {

PairBorn::PairBorn(double A, double rho, double sigma, double C, double D,
                   double cutoff)
    : Pair(cutoff), A_(A), rho_(rho), sigma_(sigma), C_(C), D_(D) {
    if (!(cutoff > 0.0))
        throw std::invalid_argument("PairBorn: cutoff must be > 0");
    if (!(rho > 0.0))
        throw std::invalid_argument("PairBorn: rho must be > 0");
    if (!(A >= 0.0))
        throw std::invalid_argument("PairBorn: A must be ≥ 0");
    if (!(C >= 0.0))
        throw std::invalid_argument("PairBorn: C must be ≥ 0");
    cutoff_sq_ = cutoff * cutoff;
    // U(rc) = A exp((sigma-rc)/rho) - C/rc^6 + D/rc^8
    const double rc2 = cutoff_sq_;
    const double rc6 = rc2 * rc2 * rc2;
    const double rc8 = rc6 * rc2;
    v_shift_ = A * std::exp((sigma - cutoff) / rho) - C / rc6 + D / rc8;
}

void PairBorn::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = static_cast<int>(f.nlocal());

    double energy = 0.0;
    const StorageReal* __restrict x = f.span(kAtom_x);
    const StorageReal* __restrict y = f.span(kAtom_y);
    const StorageReal* __restrict z = f.span(kAtom_z);
    AccReal* __restrict fx = f.span(kAtom_fx);
    AccReal* __restrict fy = f.span(kAtom_fy);
    AccReal* __restrict fz = f.span(kAtom_fz);

    const double inv_rho = 1.0 / rho_;

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
            const double r6 = r2 * r2 * r2;
            const double r8 = r6 * r2;
            const double expterm = std::exp((sigma_ - r) * inv_rho);
            const double v = A_ * expterm - C_ / r6 + D_ / r8 - v_shift_;
            // dU/dr (unshifted) = -(A/rho) exp((sigma-r)/rho) + 6 C / r^7
            //                    - 8 D / r^9
            const double dUdr =
                -A_ * inv_rho * expterm + 6.0 * C_ / (r6 * r) - 8.0 * D_ / (r8 * r);
            // Convention (Yukawa/Buck/Morse): fpair = -dUdr/r and
            // F_i -= fpair * (r_j-r_i)  ⇔  F_i = (dUdr/r) * (r_j-r_i).
            const double fpair = -dUdr / r;

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
