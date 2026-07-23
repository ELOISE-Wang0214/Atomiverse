/// @file regressions/md/piston-electrode-pz/reg_piston-electrode-pz.cpp
/// @brief E2E: rigid upper-wall piston holds mean \(P_z\) at target (E6).
///
/// Public API only. Offline MetalWalls-algorithm pin (no live MW binary).
/// See REFERENCE.md. Exits non-zero on failure.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/dynamics/fix_freeze.h"
#include "atomiverse/cpu/dynamics/init.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/nose_hoover.h"
#include "atomiverse/cpu/dynamics/piston_electrode.h"
#include "atomiverse/cpu/dynamics/velocity_verlet.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

int main() {
    // ── Parallel-plate LJ fluid + frozen lower / piston upper (REFERENCE) ─
    constexpr double Lxy = 10.0; // [Å]
    constexpr double Lz = 20.0;
    constexpr double z_lo = 2.0;  // lower wall
    constexpr double z_hi = 8.0;  // upper piston (tight gap → repulsive LJ)
    constexpr double eps = 0.10;  // [eV]
    constexpr double sigma = 3.0; // [Å]
    constexpr double cutoff = 6.0;
    constexpr double skin = 0.5;
    constexpr double dt = 0.5; // [fs]
    constexpr double T = 300.0;
    constexpr double tau_T = 50.0;
    constexpr int n_wall = 4; // 2×2 per wall
    constexpr int n_fluid = 16;
    constexpr int n_warm = 400;
    constexpr int n_prod = 600;
    constexpr int lower_type = 2;
    constexpr int upper_type = 3;
    constexpr int fluid_type = 1;

    const int N = 2 * n_wall + n_fluid;
    cpu::Frame frame(N);
    int idx = 0;

    // Lower wall (frozen)
    for (int ix = 0; ix < 2; ++ix)
        for (int iy = 0; iy < 2; ++iy) {
            frame.span(kAtom_x)[idx] = (ix + 0.5) * (Lxy / 2.0);
            frame.span(kAtom_y)[idx] = (iy + 0.5) * (Lxy / 2.0);
            frame.span(kAtom_z)[idx] = z_lo;
            frame.span(kAtom_type)[idx] = lower_type;
            frame.span(kAtom_mass)[idx] = 196.9665;
            frame.span(kAtom_tag)[idx] = idx + 1;
            frame.span(kAtom_vx)[idx] = 0.0;
            frame.span(kAtom_vy)[idx] = 0.0;
            frame.span(kAtom_vz)[idx] = 0.0;
            ++idx;
        }
    // Upper piston wall
    for (int ix = 0; ix < 2; ++ix)
        for (int iy = 0; iy < 2; ++iy) {
            frame.span(kAtom_x)[idx] = (ix + 0.5) * (Lxy / 2.0);
            frame.span(kAtom_y)[idx] = (iy + 0.5) * (Lxy / 2.0);
            frame.span(kAtom_z)[idx] = z_hi;
            frame.span(kAtom_type)[idx] = upper_type;
            frame.span(kAtom_mass)[idx] = 12.011;
            frame.span(kAtom_tag)[idx] = idx + 1;
            frame.span(kAtom_vx)[idx] = 0.0;
            frame.span(kAtom_vy)[idx] = 0.0;
            frame.span(kAtom_vz)[idx] = 0.0;
            ++idx;
        }
    // Dense fluid slab between walls — two layers near each wall so LJ
    // repulsion on the upper piston is positive (compressive).
    for (int i = 0; i < n_fluid; ++i) {
        const int ix = i % 4;
        const int iy = (i / 4) % 2;
        const int iz = i / 8; // 0 or 1
        frame.span(kAtom_x)[idx] = (ix + 0.5) * (Lxy / 4.0);
        frame.span(kAtom_y)[idx] = (iy + 0.5) * (Lxy / 2.0);
        // Layer 0 just above lower wall; layer 1 just below upper wall.
        frame.span(kAtom_z)[idx] = (iz == 0) ? (z_lo + 2.2) : (z_hi - 2.2);
        frame.span(kAtom_type)[idx] = fluid_type;
        frame.span(kAtom_mass)[idx] = 39.948;
        frame.span(kAtom_tag)[idx] = idx + 1;
        ++idx;
    }
    frame.simbox = cpu::Box::orthorhombic(Lxy, Lxy, Lz);

    // Snapshot lower wall for freeze check
    std::vector<double> lo_z0;
    for (int i = 0; i < N; ++i)
        if (frame.span(kAtom_type)[i] == lower_type)
            lo_z0.push_back(frame.span(kAtom_z)[i]);

    // Probe atom forces once: piston with P_target=0 measures P_atoms.
    cpu::md::PistonElectrode probe(upper_type, /*P=*/0.0,
                                   cpu::md::PistonElectrode::Side::Upper);
    {
        cpu::md::MD probe_sim;
        probe_sim.add_potential(cpu::PairLJCut{eps, sigma, cutoff});
        probe_sim.set_neighbor(cpu::md::NeighborList<>{cutoff, skin});
        probe_sim.set_integrator(cpu::md::VelocityVerlet{dt});
        probe_sim.add_ref(probe);
        probe_sim.enable_timing(false);
        probe_sim.build();
        probe_sim.run(frame, 1);
    }
    const double P_pin = probe.last_P_z();
    std::printf("Probe P_z (atoms/A) = %.6e eV/Å³  F_atoms=%.6e  A=%.3f\n",
                P_pin, probe.last_F_z_atoms(), probe.last_area());
    if (!(std::isfinite(P_pin) && std::fabs(P_pin) > 1e-8)) {
        std::printf("FAIL: probe |P_z| must be nonzero (got %.3e)\n", P_pin);
        return 1;
    }
    if (!(P_pin > 0.0)) {
        std::printf("FAIL: expected compressive probe P_z>0 on upper wall "
                    "(got %.3e); geometry not repulsive enough\n",
                    P_pin);
        return 1;
    }

    // Reset velocities after probe (VV advanced one step)
    for (int i = 0; i < frame.nlocal(); ++i) {
        if (frame.span(kAtom_type)[i] == fluid_type)
            continue;
        frame.span(kAtom_vx)[i] = 0.0;
        frame.span(kAtom_vy)[i] = 0.0;
        frame.span(kAtom_vz)[i] = 0.0;
    }
    // Restore wall z after probe drift (piston moved under unbalanced P=0)
    for (int i = 0; i < frame.nlocal(); ++i) {
        if (frame.span(kAtom_type)[i] == upper_type)
            frame.span(kAtom_z)[i] = z_hi;
        if (frame.span(kAtom_type)[i] == lower_type)
            frame.span(kAtom_z)[i] = z_lo;
    }

    // Production: P_target = probe pin, damped piston, NVT fluid.
    cpu::md::PistonElectrode piston(upper_type, P_pin,
                                    cpu::md::PistonElectrode::Side::Upper);
    piston.tau_damp = 20.0; // [fs] textbook equilibration

    cpu::md::MD sim;
    sim.add_potential(cpu::PairLJCut{eps, sigma, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff, skin});
    sim.set_integrator(cpu::md::VelocityVerlet{dt});
    sim.set_thermostat(cpu::md::NoseHoover{T, tau_T, dt});
    sim.add(cpu::md::InitVelocity{T, /*seed=*/42u});
    sim.add(cpu::md::FixFreeze{lower_type});
    sim.add_ref(piston);
    sim.enable_timing(false);
    sim.build();

    // Force-balance residual at first force eval with P_target=P_pin.
    // (One step: InitVelocity + integrate + force + piston post_force.)
    sim.run(frame, 1);
    const double F_atoms = piston.last_F_z_atoms();
    const double F_net = piston.last_F_net();
    const double res =
        std::fabs(F_net) / (std::fabs(F_atoms) + std::fabs(P_pin * piston.last_area()) + 1e-30);
    std::printf("Step1 force residual |F_net|/(|F_a|+|F_ext|)=%.3e  P_z=%.6e\n",
                res, piston.last_P_z());
    // With damping and velocity after one step, residual may include F_damp;
    // ac-005 is the static identity (unit-tested). E2E gates mean P_z.

    const auto t0 = std::chrono::steady_clock::now();
    // Warmup
    sim.run(frame, n_warm);
    // Production average of piston P_z
    double sum_P = 0.0;
    int n_samp = 0;
    constexpr int chunk = 20;
    for (int s = 0; s < n_prod; s += chunk) {
        sim.run(frame, chunk);
        sum_P += piston.last_P_z();
        ++n_samp;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double wall_s = std::chrono::duration<double>(t1 - t0).count();
    const double mean_P = sum_P / static_cast<double>(n_samp);

    std::printf("mean P_z = %.6e  target = %.6e  rel_err = %.3e  wall=%.4f s\n",
                mean_P, P_pin, std::fabs(mean_P - P_pin) / P_pin, wall_s);

    // ── Gates (REFERENCE.md) ──────────────────────────────────────────
    const double rel_err = std::fabs(mean_P - P_pin) / std::fabs(P_pin);
    if (!(rel_err <= 0.10)) {
        std::printf("FAIL: |<P_z>-P_target|/|P_target| = %.3e > 0.10\n",
                    rel_err);
        return 1;
    }

    // Rigid upper wall: all type-3 atoms share the same z within 1e-9
    double zmin = 1e300, zmax = -1e300;
    double max_vxy = 0.0;
    int n_up = 0;
    for (int i = 0; i < frame.nlocal(); ++i) {
        if (frame.span(kAtom_type)[i] != upper_type)
            continue;
        const double z = frame.span(kAtom_z)[i];
        zmin = std::min(zmin, z);
        zmax = std::max(zmax, z);
        max_vxy = std::max(max_vxy, std::fabs(static_cast<double>(frame.span(kAtom_vx)[i])));
        max_vxy = std::max(max_vxy, std::fabs(static_cast<double>(frame.span(kAtom_vy)[i])));
        ++n_up;
    }
    if (n_up != n_wall) {
        std::printf("FAIL: expected %d upper wall atoms, saw %d\n", n_wall, n_up);
        return 1;
    }
    if (!(zmax - zmin <= 1e-9)) {
        std::printf("FAIL: upper wall not coplanar Δz=%.3e\n", zmax - zmin);
        return 1;
    }
    if (!(max_vxy <= 1e-12)) {
        std::printf("FAIL: upper wall |v_xy| max=%.3e\n", max_vxy);
        return 1;
    }

    // Lower wall frozen
    int ilo = 0;
    double max_lo_dz = 0.0;
    for (int i = 0; i < frame.nlocal(); ++i) {
        if (frame.span(kAtom_type)[i] != lower_type)
            continue;
        max_lo_dz =
            std::max(max_lo_dz, std::fabs(frame.span(kAtom_z)[i] - lo_z0[static_cast<std::size_t>(ilo)]));
        ++ilo;
    }
    if (!(max_lo_dz <= 1e-10)) {
        std::printf("FAIL: lower wall moved Δz=%.3e\n", max_lo_dz);
        return 1;
    }

    // Finite thermo
    const double E = frame.meta.get(kEnergy);
    if (!std::isfinite(E)) {
        std::printf("FAIL: non-finite energy\n");
        return 1;
    }

    std::printf("PASS piston-electrode-pz  mean_P_z=%.6e  target=%.6e  "
                "rel_err=%.3e  wall_s=%.4f\n",
                mean_P, P_pin, rel_err, wall_s);
    return 0;
}
