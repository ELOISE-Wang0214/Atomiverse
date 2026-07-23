/// @file regressions/sampling/fep_soft_dimer/reg_fep_soft_dimer.cpp
/// Kirkwood TI identity on fixed soft LJ dimer.
#include <cmath>
#include <cstdio>

#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/pair_list.h"
#include "atomiverse/cpu/potentials/pair_lj_soft.h"

using namespace ATV;

static double energy_at(double lam, cpu::md::PairList& pairs) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 3.0;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1.0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 0;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.simbox = cpu::Box::cubic(50.0);
    f.meta.get(kEnergy) = 0;
    f.meta.get(kEnergyPair) = 0;
    for (int i = 0; i < 2; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    cpu::Context ctx;
    ctx.md.pair_list = &pairs;
    cpu::PairLJSoft soft(0.5, 3.0, 8.0, 0.5, 1.0);
    soft.lambda = lam;
    soft.build({ctx, f});
    soft.launch({ctx, f, 0, 0});
    soft.destroy({ctx, f});
    return f.meta.get(kEnergy);
}

int main() {
    cpu::md::PairList pairs;
    pairs.inum = 1;
    pairs.ilist = {0};
    pairs.numneigh = {1, 0};
    pairs.firstneigh = {0, 0};
    pairs.neighbors = {1};
    pairs.total_pairs = 1;
    pairs.has_shifts = false;
    pairs.version = 1;

    const double U0 = energy_at(0.0, pairs);
    const double U1 = energy_at(1.0, pairs);
    const int N = 400;
    double integral = 0;
    for (int i = 0; i <= N; ++i) {
        const double lam = static_cast<double>(i) / N;
        cpu::Frame f(2);
        f.span(kAtom_x)[0] = 0;
        f.span(kAtom_y)[0] = 0;
        f.span(kAtom_z)[0] = 0;
        f.span(kAtom_x)[1] = 3.0;
        f.span(kAtom_y)[1] = 0;
        f.span(kAtom_z)[1] = 0;
        f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1.0;
        f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 0;
        f.span(kAtom_tag)[0] = 1;
        f.span(kAtom_tag)[1] = 2;
        f.simbox = cpu::Box::cubic(50.0);
        f.meta.get(kEnergy) = 0;
        f.meta.get(kEnergyPair) = 0;
        for (int a = 0; a < 2; ++a)
            f.span(kAtom_fx)[a] = f.span(kAtom_fy)[a] = f.span(kAtom_fz)[a] = 0;
        cpu::Context ctx;
        ctx.md.pair_list = &pairs;
        cpu::PairLJSoft soft(0.5, 3.0, 8.0, 0.5, 1.0);
        soft.lambda = lam;
        soft.build({ctx, f});
        soft.launch({ctx, f, 0, 0});
        const double w = (i == 0 || i == N) ? 0.5 : 1.0;
        integral += w * soft.dU_dlambda();
        soft.destroy({ctx, f});
    }
    integral *= 1.0 / N;
    std::printf("TI ∫dU/dλ=%.16e  U(1)-U(0)=%.16e\n", integral, U1 - U0);
    if (!(std::fabs(integral - (U1 - U0)) <= 1e-5)) {
        std::printf("FAIL\n");
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
