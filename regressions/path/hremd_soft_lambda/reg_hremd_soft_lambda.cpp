/// @file regressions/path/hremd_soft_lambda/reg_hremd_soft_lambda.cpp
/// Hamiltonian REMD on two-λ soft LJ dimer — exchange stats.
#include <cmath>
#include <cstdio>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/path/hremd.h"
#include "atomiverse/cpu/potentials/pair_lj_soft.h"

using namespace ATV;

int main() {
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
    f.simbox = cpu::Box::cubic(30.0);

    cpu::path::HRemd h({0.2, 0.8}, /*T=*/300.0, /*dt=*/1.0, /*steps=*/15,
                       /*seed=*/11);
    h.set_program_builder(
        [](cpu::md::MD& md, double lambda, void*) {
            cpu::PairLJSoft soft(0.0104, 3.4, 8.5);
            soft.lambda = lambda;
            md.add_potential(std::move(soft));
            md.set_neighbor(cpu::md::NeighborList<>{8.5, 0.3});
        },
        nullptr);
    h.set_energy_fn(
        [](const cpu::Frame& fr, double lambda, void*) {
            cpu::Frame copy = fr;
            cpu::Evaluate eval;
            cpu::PairLJSoft soft(0.0104, 3.4, 8.5);
            soft.lambda = lambda;
            eval.add_potential(std::move(soft));
            eval.set_neighbor(cpu::md::NeighborList<>{8.5, 0.3});
            eval.build();
            eval.run(copy);
            return copy.meta.get(kEnergy);
        },
        nullptr);

    constexpr int n_ex = 30;
    h.run(f, n_ex);
    const int attempts = h.exchange_attempts();
    const int accepts = h.exchange_accepts();
    std::printf("H-REMD attempts=%d accepts=%d ratio=%.3f nlocal=%d\n",
                attempts, accepts, attempts ? double(accepts) / attempts : 0.0,
                f.nlocal());
    if (attempts < n_ex / 2 || f.nlocal() != 2) {
        std::printf("FAIL\n");
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
