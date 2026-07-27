/// @file tests/cpu/path/test_hremd.cpp
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/path/hremd.h"
#include "atomiverse/cpu/potentials/pair_lj_soft.h"

using namespace ATV;

namespace {

void build_soft(cpu::md::MD& md, double lambda, void*) {
    cpu::PairLJSoft soft(0.0104, 3.4, 8.5);
    soft.lambda = lambda;
    md.set_evaluator(std::move(soft));
    // No neighbor for soft all-pairs path? PairLJSoft may need neighbor —
    // check: if it requires pair_list, set neighbor.
}

double energy_soft(const cpu::Frame& frame, double lambda, void*) {
    cpu::Frame f = frame; // copy for evaluate
    cpu::Evaluate eval;
    cpu::PairLJSoft soft(0.0104, 3.4, 8.5);
    soft.lambda = lambda;
    eval.set_evaluator(std::move(soft));
    eval.build();
    eval.run(f);
    return f.meta.get(kEnergy);
}

} // namespace

TEST(HRemd, ConstructsAndRunsTwoLambda) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 4.0;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_mass)[i] = 39.948;
        f.span(kAtom_type)[i] = 0;
        f.span(kAtom_tag)[i] = i + 1;
    }
    f.box = cpu::Box::cubic(30.0);

    // PairLJSoft may need NeighborList — try without first; if MD validate fails
    // we'll fix in build.
    cpu::path::HRemd h({0.0, 1.0}, /*T=*/300.0, /*dt=*/1.0, /*steps=*/5,
                       /*seed=*/3);
    h.set_program_builder(
        [](cpu::md::MD& md, double lambda, void*) {
            cpu::PairLJSoft soft(0.0104, 3.4, 8.5);
            soft.lambda = lambda;
            md.set_evaluator(std::move(soft));
            md.set_neighbor(cpu::md::NeighborList<>{8.5, 0.3});
        },
        nullptr);
    h.set_energy_fn(
        [](const cpu::Frame& fr, double lambda, void*) {
            cpu::Frame copy = fr;
            cpu::Evaluate eval;
            cpu::PairLJSoft soft(0.0104, 3.4, 8.5);
            soft.lambda = lambda;
            eval.set_evaluator(std::move(soft));
            eval.set_neighbor(cpu::md::NeighborList<>{8.5, 0.3});
            eval.build();
            eval.run(copy);
            return copy.meta.get(kEnergy);
        },
        nullptr);
    h.run(f, 4);
    EXPECT_GE(h.exchange_attempts(), 2);
    EXPECT_EQ(f.nlocal(), 2);
    (void)build_soft;
    (void)energy_soft;
}
