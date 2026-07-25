/// @file cpu/dynamics/fix_dfield.cpp
#include "atomiverse/cpu/dynamics/fix_dfield.h"

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"

#include <stdexcept>

namespace ATV {
namespace cpu {
namespace md {

void FixDfield::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    if (!f.simbox.is_orthorhombic()) {
        throw std::runtime_error(
            "FixDfield: triclinic boxes are not supported in this foothold");
    }
    auto x = f.span(kAtom_x);
    auto y = f.span(kAtom_y);
    auto z = f.span(kAtom_z);
    auto q = f.span(kAtom_charge);
    auto fx = f.span(kAtom_fx);
    auto fy = f.span(kAtom_fy);
    auto fz = f.span(kAtom_fz);
    auto img = f.span(kAtom_image);
    const int n = f.nlocal();
    const double Lx = f.simbox.length(0);
    const double Ly = f.simbox.length(1);
    const double Lz = f.simbox.length(2);
    const double V = f.simbox.volume();

    double U = 0.0;
    double Pz = 0.0;
    for (int i = 0; i < n; ++i) {
        const double qi = static_cast<double>(q[i]);
        if (qi == 0.0)
            continue;
        // Unwrapped r (orthorhombic)
        const double ix = static_cast<double>(img[3 * i + 0]);
        const double iy = static_cast<double>(img[3 * i + 1]);
        const double iz = static_cast<double>(img[3 * i + 2]);
        const double rx = static_cast<double>(x[i]) + ix * Lx;
        const double ry = static_cast<double>(y[i]) + iy * Ly;
        const double rz = static_cast<double>(z[i]) + iz * Lz;
        // Foothold: E ≡ D
        fx[i] += static_cast<AccReal>(qi * Dx);
        fy[i] += static_cast<AccReal>(qi * Dy);
        fz[i] += static_cast<AccReal>(qi * Dz);
        U -= qi * (Dx * rx + Dy * ry + Dz * rz);
        Pz += qi * rz;
    }
    last_u_ = U;
    last_Pz_ = (V > 0.0) ? Pz / V : 0.0;
    f.meta.get(kEnergy) += U;
}

} // namespace md
} // namespace cpu
} // namespace ATV
