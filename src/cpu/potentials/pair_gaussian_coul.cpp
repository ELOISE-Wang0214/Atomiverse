/// @file src/cpu/potentials/pair_gaussian_coul.cpp
#include "atomiverse/cpu/potentials/pair_gaussian_coul.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/potentials/gaussian_charge.h"

namespace ATV {
namespace cpu {

PairGaussianCoul::PairGaussianCoul(double cutoff,
                                   std::vector<double> eta_by_type,
                                   bool include_self)
    : Pair(cutoff), eta_by_type_(std::move(eta_by_type)),
      include_self_(include_self) {
    if (!(cutoff > 0.0))
        throw std::invalid_argument("PairGaussianCoul: cutoff must be > 0 [Å]");
    cutoff_sq_ = cutoff * cutoff;
}

void PairGaussianCoul::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = static_cast<int>(f.nlocal());

    const StorageReal* __restrict x = f.span(kAtom_x);
    const StorageReal* __restrict y = f.span(kAtom_y);
    const StorageReal* __restrict z = f.span(kAtom_z);
    const StorageReal* __restrict q = f.span(kAtom_charge);
    const LocalInt* __restrict type = f.span(kAtom_type);
    AccReal* __restrict fx = f.span(kAtom_fx);
    AccReal* __restrict fy = f.span(kAtom_fy);
    AccReal* __restrict fz = f.span(kAtom_fz);

    const double K = ATV::units::coulomb;
    double energy = 0.0;

    for (int i = 0; i < n; ++i) {
        const double qi = static_cast<double>(q[i]);
        const double ei = eta_of(static_cast<int>(type[i]));

        if (include_self_ && GaussianChargeKernel::is_gaussian(ei))
            energy += GaussianChargeKernel::self_energy(K, qi, ei);

        for (int j = i + 1; j < n; ++j) {
            double dx = static_cast<double>(x[j] - x[i]);
            double dy = static_cast<double>(y[j] - y[i]);
            double dz = static_cast<double>(z[j] - z[i]);
            f.simbox.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= cutoff_sq_ || r2 == 0.0)
                continue;

            const double r = std::sqrt(r2);
            const double qj = static_cast<double>(q[j]);
            const double ej = eta_of(static_cast<int>(type[j]));

            const double G = GaussianChargeKernel::greens(r, ei, ej);
            const double dG = GaussianChargeKernel::dgreens_dr(r, ei, ej);
            const double pref = K * qi * qj;
            energy += pref * G;

            // U = pref * G(r); F_i = −∇_i U.
            // r = |r_j − r_i|, ∇_i r = −(r_j−r_i)/r ⇒ ∇_i U = pref * dG/dr * (−r̂)
            // F_i = −∇_i U = pref * dG/dr * r̂ = (pref * dG/dr / r) * (r_j−r_i)
            const double fpair = pref * dG / r;
            fx[i] += fpair * dx;
            fy[i] += fpair * dy;
            fz[i] += fpair * dz;
            fx[j] -= fpair * dx;
            fy[j] -= fpair * dy;
            fz[j] -= fpair * dz;
        }
    }

    f.meta.get(kEnergy) += energy;
}

} // namespace cpu
} // namespace ATV
