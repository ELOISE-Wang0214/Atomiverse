/// @file cpu/potentials/pair_dpd.cpp
#include "atomiverse/cpu/potentials/pair_dpd.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"

namespace ATV {
namespace cpu {

PairDpd::PairDpd(double a, double gamma, double temperature, double cutoff,
                 double dt, std::uint64_t seed)
    : Pair(cutoff), a_(a), gamma_(gamma), temperature_(temperature), dt_(dt),
      seed_(seed) {
    if (!(cutoff > 0.0) || !(dt > 0.0) || !(temperature > 0.0) || !(gamma >= 0.0))
        throw std::invalid_argument("PairDpd: invalid parameters");
    cutoff_sq_ = cutoff * cutoff;
    // σ² = 2 γ k_B T  (Groot–Warren)
    sigma_ = std::sqrt(2.0 * gamma_ * units::kB_ev_per_k * temperature_);
}

void PairDpd::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = static_cast<int>(f.nlocal());
    const double rc = pair_cutoff();
    const double inv_sqrt_dt = 1.0 / std::sqrt(dt_);

    double energy = 0.0;
    const StorageReal* __restrict x = f.span(kAtom_x);
    const StorageReal* __restrict y = f.span(kAtom_y);
    const StorageReal* __restrict z = f.span(kAtom_z);
    const StorageReal* __restrict vx = f.span(kAtom_vx);
    const StorageReal* __restrict vy = f.span(kAtom_vy);
    const StorageReal* __restrict vz = f.span(kAtom_vz);
    AccReal* __restrict fx = f.span(kAtom_fx);
    AccReal* __restrict fy = f.span(kAtom_fy);
    AccReal* __restrict fz = f.span(kAtom_fz);

    // Deterministic pair RNG: seed from (global seed, step, i, j).
    auto theta = [&](int i, int j) {
        std::uint64_t s = seed_;
        s ^= static_cast<std::uint64_t>(step_) * 0x9e3779b97f4a7c15ULL;
        s ^= static_cast<std::uint64_t>(i) * 0xbf58476d1ce4e5b9ULL;
        s ^= static_cast<std::uint64_t>(j) * 0x94d049bb133111ebULL;
        s = (s ^ (s >> 30)) * 0xbf58476d1ce4e5b9ULL;
        s = (s ^ (s >> 27)) * 0x94d049bb133111ebULL;
        s = s ^ (s >> 31);
        // Map to N(0,1) via Box–Muller on two uniform bits.
        const double u1 =
            (static_cast<double>((s & 0xffffffffULL) + 1) / 4294967297.0);
        const double u2 =
            (static_cast<double>(((s >> 32) & 0xffffffffULL) + 1) / 4294967297.0);
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * units::pi * u2);
    };

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
            const double inv_r = 1.0 / r;
            const double w = 1.0 - r / rc;
            if (w <= 0.0)
                continue;

            // Conservative: F_i ← −a w r̂  (r̂ from i to j = r_ij/r)
            // Force on i from pair: −F^C along r̂_ij if F^C is the mutual repulsion.
            // Mutual: F_i = −a w r̂, F_j = +a w r̂ with r̂ = (r_j−r_i)/r.
            const double fc = a_ * w;
            double fpair = -fc * inv_r; // F_i += fpair * r_ij

            // Dissipative: −γ w² (r̂·v_ij) r̂ ; v_ij = v_j − v_i
            const double dvx = vx[j] - vx[i];
            const double dvy = vy[j] - vy[i];
            const double dvz = vz[j] - vz[i];
            const double rdotv = (dx * dvx + dy * dvy + dz * dvz) * inv_r;
            const double w2 = w * w;
            const double fd = -gamma_ * w2 * rdotv;
            fpair += fd * inv_r;

            // Random: σ w θ / √dt  along r̂
            const double fr = sigma_ * w * theta(i, j) * inv_sqrt_dt;
            fpair += fr * inv_r;

            fx[i] += fpair * dx;
            fy[i] += fpair * dy;
            fz[i] += fpair * dz;
            fx[j] -= fpair * dx;
            fy[j] -= fpair * dy;
            fz[j] -= fpair * dz;

            // U = (a/2) w² per pair (Groot–Warren soft potential)
            energy += 0.5 * a_ * w2;
        }
    }
    f.meta.get(kEnergy) += energy;
}

} // namespace cpu
} // namespace ATV
