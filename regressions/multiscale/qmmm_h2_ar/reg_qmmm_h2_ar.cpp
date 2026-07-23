#include <cmath>
#include <cstdio>
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/qmmm_mechanical.h"
using namespace ATV;
int main() {
    cpu::Frame f(3);
    f.span(kAtom_x)[0]=0; f.span(kAtom_y)[0]=0; f.span(kAtom_z)[0]=0;
    f.span(kAtom_x)[1]=0.74; f.span(kAtom_y)[1]=0; f.span(kAtom_z)[1]=0;
    f.span(kAtom_x)[2]=5.0; f.span(kAtom_y)[2]=0; f.span(kAtom_z)[2]=0;
    f.span(kAtom_type)[0]=1; f.span(kAtom_type)[1]=1; f.span(kAtom_type)[2]=18;
    f.span(kAtom_tag)[0]=1; f.span(kAtom_tag)[1]=2; f.span(kAtom_tag)[2]=3;
    f.span(kAtom_mass)[0]=1; f.span(kAtom_mass)[1]=1; f.span(kAtom_mass)[2]=40;
    f.simbox=cpu::Box::cubic(30);
    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate)
        .add(cpu::QmmmMechanical{{1,2}, 0.01, 3.4, 10.0});
    eval.build();
    eval.run(f);
    const double E=f.meta.get(kEnergy);
    std::printf("QMMM E=%.6e\n", E);
    if (!std::isfinite(E)) { std::printf("FAIL\n"); return 1; }
    std::printf("PASS\n");
    return 0;
}
