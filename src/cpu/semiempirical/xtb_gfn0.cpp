/// @file src/cpu/semiempirical/xtb_gfn0.cpp
#include "atomiverse/cpu/semiempirical/xtb_gfn0.h"

#include <cmath>
#include <stdexcept>
#include <unordered_map>

#include "atomiverse/core/frame_schema.h"

namespace ATV {
namespace cpu {
namespace semiempirical {
namespace {

// Toy GFN0-like onsite and pair repulsion (eV, Å) — framework only.
double onsite(int Z) {
    switch (Z) {
    case 1:
        return -13.6;
    case 6:
        return -11.0;
    case 8:
        return -15.0;
    default:
        throw std::runtime_error("XtbGfn0: only H/C/O (Z=1,6,8)");
    }
}
double pair_rep(int Zi, int Zj, double r) {
    const double a = 0.3 * std::sqrt(std::abs(onsite(Zi) * onsite(Zj)));
    return a * std::exp(-2.0 * r);
}
double hop(int Zi, int Zj, double r) {
    const double h0 = 0.5 * (onsite(Zi) + onsite(Zj));
    return h0 * std::exp(-r / 1.2);
}

} // namespace

void XtbGfn0::build(const BuildContext&) {}

void XtbGfn0::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = f.nlocal();
    if (n < 1)
        throw std::invalid_argument("XtbGfn0: empty frame");
    const auto type = f.span(kAtom_type);
    const StorageReal* x = f.span(kAtom_x);
    const StorageReal* y = f.span(kAtom_y);
    const StorageReal* z = f.span(kAtom_z);
    double& E = f.meta.get(kEnergy);

    double e = 0;
    for (int i = 0; i < n; ++i)
        e += onsite(type[i]); // occupancy sketch: 1 e per atom valence toy
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = x[j] - x[i], dy = y[j] - y[i], dz = z[j] - z[i];
            f.simbox.min_image(dx, dy, dz);
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r < 1e-8)
                continue;
            // Bonding estimate: 2 * hop for filled bonding orbital (sketch)
            e += 2.0 * hop(type[i], type[j], r) + pair_rep(type[i], type[j], r);
        }
    }
    energy_ = e;
    E += e;
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
