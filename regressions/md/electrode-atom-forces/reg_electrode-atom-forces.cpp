
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include <cmath>
#include <cstdio>
using namespace ATV;

#include "atomiverse/cpu/dynamics/const_pot_bo.h"
int main() {
    cpu::Frame f(3);
    double zz[3] = {5, 15, 14};
    int ty[3] = {1, 2, 3};
    double qq[3] = {0, 0, 1};
    for (int i = 0; i < 3; ++i) {
        f.span(kAtom_x)[i] = f.span(kAtom_y)[i] = 0;
        f.span(kAtom_z)[i] = zz[i];
        f.span(kAtom_type)[i] = ty[i];
        f.span(kAtom_charge)[i] = qq[i];
        f.span(kAtom_mass)[i] = 1;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::orthorhombic(10, 10, 20);
    cpu::md::ConstPotBO bo(2.0);
    bo.add_electrode(1, 0.0);
    bo.add_electrode(2, 1.0);
    bo.write_electrode_forces = true;
    bo.xy_image_range = 20;
    cpu::Context ctx;
    bo.build({ctx, f});
    bo.launch({ctx, f, 0, 0});
    double fz1 = f.span(kAtom_fz)[0], fz2 = f.span(kAtom_fz)[1];
    std::printf("electrode-atom-forces fz_bot=%.6e fz_top=%.6e Qtop=%.6e\n",
                fz1, fz2, bo.total_charge(1));
    int fail = 0;
    if (std::abs(bo.total_charge(0)) + std::abs(bo.total_charge(1)) < 1e-12) {
        std::fprintf(stderr, "zero Q\n");
        ++fail;
    }
    if (std::abs(fz1) + std::abs(fz2) < 1e-14) {
        std::fprintf(stderr, "zero electrode forces\n");
        ++fail;
    }
    return fail ? 1 : 0;
}
