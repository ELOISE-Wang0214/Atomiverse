/// @file regressions/md/settle_tip3p/reg_settle_tip3p.cpp
#include <cmath>
#include <cstdio>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/dynamics/fix_settle.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/velocity_verlet.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

static double dist(const cpu::Frame& f, int i, int j) {
    const double dx = f.span(kAtom_x)[i] - f.span(kAtom_x)[j];
    const double dy = f.span(kAtom_y)[i] - f.span(kAtom_y)[j];
    const double dz = f.span(kAtom_z)[i] - f.span(kAtom_z)[j];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

int main() {
    const double d_oh = 0.9572;
    const double d_hh = 1.5139;
    const double half = 0.5 * d_hh;
    const double h = std::sqrt(d_oh * d_oh - half * half);
    const double dt = 2.0;

    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = h;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = -half;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_x)[2] = half;
    f.span(kAtom_y)[2] = 0;
    f.span(kAtom_z)[2] = 0;
    f.span(kAtom_mass)[0] = 15.999;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_mass)[2] = 1.008;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 2;
    f.span(kAtom_type)[2] = 2;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.span(kAtom_tag)[2] = 3;
    f.span(kAtom_vx)[0] = 0.001;
    f.span(kAtom_vy)[1] = -0.002;
    f.span(kAtom_vz)[2] = 0.0015;
    f.box = cpu::Box::cubic(30.0);

    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{dt});
    // Soft LJ on O only to add mild forces without blowing geometry.
    auto lj = cpu::PairLJCut{0.0067, 3.1507, 8.0};
    // Single type — hydrogens share O LJ; dilute box OK for geometry pin.
    sim.set_evaluator(std::move(lj));
    sim.set_neighbor(cpu::md::NeighborList<>{8.0, 0.5});
    cpu::md::FixSettle settle(dt, d_oh, d_hh);
    settle.add_water(1, 2, 3);
    sim.set_observer(std::move(settle));
    sim.build();
    sim.run(f, 500); // 1 ps at 2 fs

    const double e1 = std::fabs(dist(f, 0, 1) - d_oh) / d_oh;
    const double e2 = std::fabs(dist(f, 0, 2) - d_oh) / d_oh;
    const double e3 = std::fabs(dist(f, 1, 2) - d_hh) / d_hh;
    std::printf("SETTLE: rel err OH1=%.3e OH2=%.3e HH=%.3e\n", e1, e2, e3);
    if (!(e1 < 1e-10 && e2 < 1e-10 && e3 < 1e-10)) {
        std::printf("FAIL\n");
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
