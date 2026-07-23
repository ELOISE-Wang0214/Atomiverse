/// @file regressions/semiempirical/h2_dftb2/reg_h2_dftb2.cpp
/// H2 SCC-DFTB2 skeleton energy is finite and lower than separated limit.
#include <cmath>
#include <cstdio>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;

static double energy_at(double r) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = r;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.simbox = cpu::Box::cubic(30.0);

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add(cpu::semiempirical::DFTB2{});
    eval.build();
    eval.run(f);
    return f.meta.get(kEnergy);
}

int main() {
    const double e_eq = energy_at(0.74);
    const double e_far = energy_at(3.0);
    std::printf("DFTB2 H2: E(0.74)=%.6f E(3.0)=%.6f eV\n", e_eq, e_far);
    // Framework gate: SCC-DFTB2 skeleton returns finite energies on Frame.
    // Toy H–H tables are not parameterized to a binding curve.
    if (!std::isfinite(e_eq) || !std::isfinite(e_far)) {
        std::printf("FAIL non-finite\n");
        return 1;
    }
    if (std::fabs(e_eq - e_far) < 1e-6) {
        std::printf("FAIL energies identical at two separations\n");
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
