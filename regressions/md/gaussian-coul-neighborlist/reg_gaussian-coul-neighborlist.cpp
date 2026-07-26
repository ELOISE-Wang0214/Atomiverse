
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include <cmath>
#include <cstdio>
using namespace ATV;

#include "atomiverse/cpu/potentials/pair_gaussian_coul.h"
int main() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 2.0;
    f.span(kAtom_charge)[0] = 1;
    f.span(kAtom_charge)[1] = -1;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 0;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.box = cpu::Box::cubic(20);
    f.meta.get(kEnergy) = 0;
    cpu::PairGaussianCoul pg(10.0, {0.0, 2.0}, true);
    cpu::Context ctx;
    pg.build({ctx, f});
    for (int i = 0; i < 2; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    f.meta.get(kEnergy) = 0;
    pg.launch({ctx, f, 0, 0});
    double e_n2 = f.meta.get(kEnergy);
    std::printf("gaussian-coul-neighborlist E_n2=%.6e\n", e_n2);
    return std::isfinite(e_n2) ? 0 : 1;
}
