
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include <cmath>
#include <cstdio>
using namespace ATV;

#include "atomiverse/cpu/dynamics/const_pot_bo.h"
int main() {
    cpu::Frame f(3);
    double zz[3] = {5, 15, 10};
    int ty[3] = {1, 2, 3};
    double qq[3] = {0, 0, 1};
    for (int i = 0; i < 3; ++i) {
        f.span(kAtom_x)[i] = f.span(kAtom_y)[i] = 0;
        f.span(kAtom_z)[i] = zz[i];
        f.span(kAtom_type)[i] = ty[i];
        f.span(kAtom_charge)[i] = qq[i];
        f.span(kAtom_mass)[i] = 1;
        f.span(kAtom_tag)[i] = i + 1;
    }
    f.simbox = cpu::Box::orthorhombic(10, 10, 20);
    cpu::md::ConstPotBO bo(2.0);
    bo.add_electrode(1, 0.0);
    bo.add_electrode(2, 1.0);
    bo.electroneutrality = true;
    bo.global_electroneutrality = true;
    bo.xy_image_range = 50;
    cpu::Context ctx;
    bo.build({ctx, f});
    bo.launch({ctx, f, 0, 0});
    double sumQ = bo.total_charge(0) + bo.total_charge(1);
    double tgt = bo.neutrality_target();
    std::printf("global-electroneutrality-bo sumQ=%.6e target=%.6e\n", sumQ,
                tgt);
    int fail = 0;
    if (std::abs(tgt + 1.0) > 1e-12) {
        std::fprintf(stderr, "bad target\n");
        ++fail;
    }
    if (std::abs(sumQ - tgt) > 1e-8) {
        std::fprintf(stderr, "sumQ != target\n");
        ++fail;
    }
    return fail ? 1 : 0;
}
