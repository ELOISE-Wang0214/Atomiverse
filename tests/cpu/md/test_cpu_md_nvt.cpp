/// @file tests/cpu/md/test_cpu_md_nvt.cpp
/// @brief cpu::md::NVT — the LAMMPS `fix nvt`-style bundled integrator
///        (velocity-Verlet + selectable thermostat), added through the singular
///        set_integrator. Confirms both thermostat choices (CSVR and
///        Nose-Hoover) drive a cheap LJ argon gas to the canonical setpoint
///        temperature, i.e. the bundling composes the existing thermostat
///        physics correctly.

#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/nvt.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

namespace {

cpu::Frame make_fcc_argon(int n_cells, double a) {
    const int n = n_cells * n_cells * n_cells * 4;
    cpu::Frame frame(n);
    const double basis[4][3] = {
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
    int idx = 0;
    for (int ix = 0; ix < n_cells; ++ix)
        for (int iy = 0; iy < n_cells; ++iy)
            for (int iz = 0; iz < n_cells; ++iz)
                for (int b = 0; b < 4; ++b) {
                    frame.atoms.x[idx] = (ix + basis[b][0]) * a;
                    frame.atoms.y[idx] = (iy + basis[b][1]) * a;
                    frame.atoms.z[idx] = (iz + basis[b][2]) * a;
                    frame.atoms.mass[idx] = 39.948;
                    frame.atoms.type[idx] = 18;
                    ++idx;
                }
    frame.box = Mat3::cubic(n_cells * a);
    return frame;
}

double kinetic_ev(const cpu::Frame& f) {
    double s = 0.0;
    for (int i = 0; i < f.atoms.size(); ++i)
        s += f.atoms.mass[i]
             * (f.atoms.vx[i] * f.atoms.vx[i] + f.atoms.vy[i] * f.atoms.vy[i]
                + f.atoms.vz[i] * f.atoms.vz[i]);
    return 0.5 * kMvv2e * s;
}

// Equilibrate a 108-atom LJ argon gas under NVT (the given thermostat) and
// return the mean temperature sampled over the production window.
double nvt_mean_temperature(cpu::md::NVT::Thermostat thermo) {
    const double T = 300.0, dt = 5.0, tau = 100.0;
    const double eps = 0.0104, sigma = 3.4, cutoff = 8.5;

    auto frame = make_fcc_argon(3, 5.26); // 108 atoms
    const int n_dof = 3 * frame.atoms.size() - 3;

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T, /*seed=*/2024});
    // One call bundles VelocityVerlet + the chosen thermostat (LAMMPS fix nvt
    // ergonomic) — no separate set_integrator + set_thermostat.
    sim.set_nvt(cpu::md::NVT{dt, T, tau, thermo, /*seed=*/7});
    sim.set_pair(cpu::PairLJCut{eps, sigma, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff, 1.0});
    sim.build();

    sim.run(frame, 2500); // equilibrate

    const int n_samples = 600, interval = 8;
    double mean_ke = 0.0;
    for (int s = 0; s < n_samples; ++s) {
        sim.run(frame, interval);
        mean_ke += kinetic_ev(frame);
    }
    mean_ke /= n_samples;
    return 2.0 * mean_ke / (n_dof * kBoltzmann);
}

} // namespace

TEST(MdNVT, CsvrReachesSetpoint) {
    const double mean_T = nvt_mean_temperature(cpu::md::NVT::Thermostat::CSVR);
    std::printf("[NVT/CSVR] mean_T = %.2f K (target 300)\n", mean_T);
    EXPECT_NEAR(mean_T, 300.0, 0.08 * 300.0)
        << "NVT(CSVR) did not thermostat to the setpoint";
}

TEST(MdNVT, NoseHooverReachesSetpoint) {
    const double mean_T =
        nvt_mean_temperature(cpu::md::NVT::Thermostat::NoseHoover);
    std::printf("[NVT/NoseHoover] mean_T = %.2f K (target 300)\n", mean_T);
    EXPECT_NEAR(mean_T, 300.0, 0.08 * 300.0)
        << "NVT(NoseHoover) did not thermostat to the setpoint";
}
