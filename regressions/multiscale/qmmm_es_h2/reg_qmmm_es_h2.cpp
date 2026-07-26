#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"
#include <cmath>
#include <cstdio>
using namespace ATV;
static double energy(double q) {
    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_x)[2] = 2.5;
    f.span(kAtom_y)[2] = 0;
    f.span(kAtom_z)[2] = 0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.span(kAtom_tag)[2] = 3;
    f.box = cpu::Box::cubic(20);
    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate)
        .add(cpu::QmmmElectrostatic{{1, 2}, q, 0.01, 3.4, 8.0});
    eval.build();
    eval.run(f);
    return f.meta.get(kEnergy);
}
int main() {
    const double e0 = energy(0.0), e1 = energy(1.0);
    std::printf("QMMM-ES E(q=0)=%.6f E(q=1)=%.6f eV\n", e0, e1);
    if (!std::isfinite(e0) || !std::isfinite(e1) || std::fabs(e1 - e0) < 1e-4) {
        std::printf("FAIL\n");
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
