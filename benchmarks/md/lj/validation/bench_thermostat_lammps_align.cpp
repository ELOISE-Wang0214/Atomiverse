/// @file benchmarks/md/bench_thermostat_lammps_align.cpp
/// @brief lammps-api-align-03 validation benchmark — runtime checks that
///        FixLangevin and NoseHoover match LAMMPS `fix langevin` /
///        `fix nvt` semantics: the `Tstart`→`Tstop` linear temperature
///        ramp and the damp / tau_T parameter being a *time constant*
///        (not a rate). Ported from the gtest unit-test suite; this is a
///        third-party-comparison validation benchmark, not a ctest gate.
///
/// Acceptance: ac-004 (ramp interpolates linearly across the run span),
/// ac-005 (halving damp/tau_T roughly doubles the relaxation rate),
/// ac-006 (NVT equilibration holds the LAMMPS-equivalent mean T).
///
/// Both thermostats expose, per the spec:
///   * FixLangevin ctor `(T_start, T_stop, damp, dt, seed)` — LAMMPS
///     `fix langevin Tstart Tstop damp seed` argument order;
///   * NoseHoover ctor `(T_start, tau_T, dt)` (+ `set_T_stop` /
///     `set_run_span` for the ramp) — LAMMPS `fix nvt temp`;
///   * `set_run_span(nfirst, nlast)` on both, mirroring LAMMPS
///     `update->beginstep` / `endstep` captured at run setup.

#include "bench_check.h"

#include "atomiverse/core/constants.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/fix_langevin.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/velocity_verlet.h"

#include <cmath>
#include <iostream>

using namespace ATV;

namespace {

/// @brief Build an FCC argon box (`n_cells^3 · 4` atoms) with no
/// initial velocity — InitVelocity seeds them in the driver.
cpu::Frame make_fcc_argon(int n_cells, double a, double mass) {
    constexpr int atoms_per_cell = 4;
    const int n_atoms = n_cells * n_cells * n_cells * atoms_per_cell;
    cpu::Frame frame(n_atoms);
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
                    frame.atoms.mass[idx] = mass;
                    frame.atoms.type[idx] = 0;
                    ++idx;
                }
    frame.box = Mat3::cubic(n_cells * a);
    return frame;
}

/// @brief Instantaneous temperature with 3(N-1) DOF (COM removed),
/// matching the LAMMPS `temp` compute.
double instantaneous_T(const cpu::Frame& f) {
    const int n = f.atoms.size();
    double KE = 0.0;
    for (int i = 0; i < n; ++i) {
        const double v2 = f.atoms.vx[i] * f.atoms.vx[i]
                          + f.atoms.vy[i] * f.atoms.vy[i]
                          + f.atoms.vz[i] * f.atoms.vz[i];
        KE += 0.5 * f.atoms.mass[i] * v2;
    }
    KE *= kMvv2e;
    return 2.0 * KE / (3.0 * (n - 1) * kBoltzmann);
}

} // namespace

// ── ac-002: FixLangevin LAMMPS-ordered ctor + default no-ramp ─────────
//
// The ctor is `(T_start, T_stop, damp, dt, seed)`. With `T_stop ==
// T_start` the run is a constant-temperature thermostat — no behavioural
// change from the pre-spec single-temperature callers.
BENCH_CASE(ThermostatLammpsAlign_FixLangevinConstantTemperatureHoldsSetpoint) {
    constexpr double T = 300.0;
    constexpr double damp = 50.0; // fs
    constexpr double dt = 1.0;    // fs

    auto frame = make_fcc_argon(4, 4.2, 16.0);

    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::InitVelocity{T, /*seed=*/777u});
    // LAMMPS argument order: Tstart Tstop damp seed.
    sim.add(cpu::md::FixLangevin{T, T, damp, dt, /*seed=*/904297u});
    sim.build();
    sim.run(frame, 3000);

    BM_NEAR(instantaneous_T(frame), T, 0.15 * T)
        << "constant-T FixLangevin failed to hold the setpoint";
}

// ── ac-004: FixLangevin Tstart→Tstop linear ramp ─────────────────────
//
// With `T_stop != T_start` and a run span set via `set_run_span`, the
// live thermostat setpoint each step is `T_start + delta·(T_stop −
// T_start)`, `delta = (step − nfirst)/(nlast − nfirst) ∈ [0,1]`, exactly
// as LAMMPS `fix_langevin.cpp::compute_target()`. The measured T must
// sweep with the ramp.
BENCH_CASE(ThermostatLammpsAlign_FixLangevinTemperatureRampTracksSetpoint) {
    constexpr double T_start = 150.0;
    constexpr double T_stop = 450.0;
    constexpr double damp = 25.0; // fs (stiff coupling so T tracks)
    constexpr double dt = 1.0;    // fs
    constexpr int total_steps = 12000;

    auto frame = make_fcc_argon(4, 4.2, 16.0);

    cpu::md::FixLangevin lang{T_start, T_stop, damp, dt, /*seed=*/904297u};
    lang.set_run_span(0, total_steps);

    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::InitVelocity{T_start, /*seed=*/777u});
    sim.add(lang); // lvalue → non-owning storage
    sim.build();

    constexpr int n_windows = 6;
    constexpr int window = total_steps / n_windows;
    int step = 0;
    for (int w = 0; w < n_windows; ++w) {
        sim.run(frame, window);
        step += window;
        const double delta = static_cast<double>(step) / total_steps;
        const double setpoint = T_start + delta * (T_stop - T_start);
        const double T_meas = instantaneous_T(frame);
        std::cerr << "langevin ramp step " << step
                  << ": setpoint = " << setpoint << " K, measured = " << T_meas
                  << " K\n";
        BM_NEAR(T_meas, setpoint, 60.0)
            << "FixLangevin T does not track the linear ramp at fraction "
            << delta;
    }
    // A no-ramp run would stay near T_start; the ramp must have heated it.
    BM_GT(instantaneous_T(frame), 300.0)
        << "FixLangevin ramp did not heat the system toward T_stop";
}

// ── ac-005: FixLangevin damp behaves as a time constant ──────────────
//
// `damp` is a friction *time*; γ = 1/damp. Halving `damp` doubles γ and
// so roughly doubles the temperature-relaxation rate. Starting cold
// (T0 ≈ 0) and thermostatting to `T_target`, the approach is
// exponential `T(t) = T_target·(1 − e^{−2t/damp})`. After a fixed short
// time the run with the smaller `damp` must be measurably hotter.
BENCH_CASE(ThermostatLammpsAlign_FixLangevinDampIsTimeConstant) {
    constexpr double T_target = 400.0;
    constexpr double dt = 1.0;      // fs
    constexpr int relax_steps = 40; // short window, mid-relaxation

    auto run_relax = [&](double damp) {
        auto frame = make_fcc_argon(4, 4.2, 16.0);
        // Start essentially at rest (no InitVelocity) so the only heating
        // is the Langevin random force; the rise rate then isolates damp.
        cpu::md::MD sim;
        sim.set_integrator(cpu::md::VelocityVerlet{dt});
        sim.add(cpu::md::FixLangevin{T_target, T_target, damp, dt,
                                     /*seed=*/55u});
        sim.build();
        sim.run(frame, relax_steps);
        return instantaneous_T(frame);
    };

    const double T_fast = run_relax(20.0); // small damp → fast heating
    const double T_slow = run_relax(40.0); // large damp → slow heating

    std::cerr << "langevin damp=20 → T=" << T_fast << " K, "
              << "damp=40 → T=" << T_slow << " K\n";

    // Faster coupling heats more in the same wall-time window.
    BM_GT(T_fast, T_slow)
        << "halving damp did not speed up the relaxation — damp is not "
           "behaving as a time constant";
    // Exponential rise: ratio of (1−e^{−2t/damp}) for damp 20 vs 40 at
    // t = 40 fs is (1−e^{−4})/(1−e^{−2}) ≈ 1.13 — require at least a
    // clear 5% separation to rule out noise.
    BM_GT(T_fast / std::max(T_slow, 1e-9), 1.05)
        << "relaxation-rate ratio too small for a halved time constant";
}

// ── ac-003 + ac-004: NoseHoover Tstart→Tstop linear ramp ─────────────
//
// NoseHoover exposes `tau_T` (renamed from `tau`) and a `T_stop` ramp.
// With a run span set, the live setpoint sweeps linearly, mirroring
// LAMMPS `fix_nh.cpp::compute_temp_target()`.
BENCH_CASE(ThermostatLammpsAlign_NoseHooverTemperatureRampTracksSetpoint) {
    constexpr double T_start = 200.0;
    constexpr double T_stop = 500.0;
    constexpr double tau_T = 30.0; // fs
    constexpr double dt = 1.0;     // fs
    constexpr int total_steps = 12000;
    // FCC argon: a = 5.79 Å, mass 39.948 amu, LJ ε/σ/cutoff.
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    auto frame = make_fcc_argon(4, 5.79, 39.948);

    cpu::md::NoseHoover nh{T_start, tau_T, dt};
    nh.set_T_stop(T_stop);
    nh.set_run_span(0, total_steps);

    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::InitVelocity{T_start, /*seed=*/4242u});
    sim.add(nh); // lvalue → non-owning storage
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    sim.set_pair(cpu::PairLJCut{kEpsilon, kSigma, kCutoff});
    sim.build();

    constexpr int n_windows = 6;
    constexpr int window = total_steps / n_windows;
    int step = 0;
    for (int w = 0; w < n_windows; ++w) {
        sim.run(frame, window);
        step += window;
        const double delta = static_cast<double>(step) / total_steps;
        const double setpoint = T_start + delta * (T_stop - T_start);
        const double T_meas = instantaneous_T(frame);
        std::cerr << "nvt ramp step " << step << ": setpoint = " << setpoint
                  << " K, measured = " << T_meas << " K\n";
        BM_NEAR(T_meas, setpoint, 70.0)
            << "NoseHoover T does not track the linear ramp at fraction "
            << delta;
    }
    BM_GT(instantaneous_T(frame), 350.0)
        << "NoseHoover ramp did not heat the system toward T_stop";
}

// ── ac-005: NoseHoover tau_T behaves as a time constant ──────────────
//
// `tau_T` (LAMMPS `Tdamp`) is the thermostat coupling *time*. Halving it
// stiffens the coupling and roughly doubles the temperature-relaxation
// rate. Start hot (2·T_target) and thermostat down: after a fixed
// window the stiffer (smaller tau_T) run must be measurably colder.
BENCH_CASE(ThermostatLammpsAlign_NoseHooverTauIsTimeConstant) {
    constexpr double T_target = 200.0;
    constexpr double T_init = 400.0;
    constexpr double dt = 1.0; // fs
    constexpr int relax_steps = 200;
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    auto run_relax = [&](double tau_T) {
        auto frame = make_fcc_argon(4, 5.79, 39.948);
        cpu::md::MD sim;
        sim.set_integrator(cpu::md::VelocityVerlet{dt});
        sim.add(cpu::md::InitVelocity{T_init, /*seed=*/909u});
        sim.add(cpu::md::NoseHoover{T_target, tau_T, dt});
        sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
        sim.set_pair(cpu::PairLJCut{kEpsilon, kSigma, kCutoff});
        sim.build();
        sim.run(frame, relax_steps);
        return instantaneous_T(frame);
    };

    const double T_stiff = run_relax(20.0); // small tau_T → fast cooling
    const double T_soft = run_relax(60.0);  // large tau_T → slow cooling

    std::cerr << "nvt tau_T=20 → T=" << T_stiff << " K, "
              << "tau_T=60 → T=" << T_soft << " K\n";

    // The stiffer coupling pulls the temperature toward the setpoint
    // faster — the smaller-tau run is closer to T_target.
    BM_LT(std::abs(T_stiff - T_target), std::abs(T_soft - T_target))
        << "halving tau_T did not speed up the NVT relaxation — tau_T is "
           "not behaving as a time constant";
}

// ── ac-006: NVT equilibration holds the LAMMPS-equivalent mean T ─────
//
// A constant-T NoseHoover run reaches and holds the target mean
// temperature, consistent with LAMMPS `fix nvt` at the same setpoint.
BENCH_CASE(ThermostatLammpsAlign_NoseHooverHoldsMeanTemperature) {
    constexpr double T_target = 300.0;
    constexpr double tau_T = 100.0; // fs
    constexpr double dt = 1.0;      // fs
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    auto frame = make_fcc_argon(4, 5.79, 39.948);

    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/2024u});
    sim.add(cpu::md::NoseHoover{T_target, tau_T, dt});
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    sim.set_pair(cpu::PairLJCut{kEpsilon, kSigma, kCutoff});
    sim.build();
    sim.run(frame, 2000); // burn-in

    // Block-average the temperature over the sampling window.
    double T_sum = 0.0;
    int samples = 0;
    for (int c = 0; c < 40; ++c) {
        sim.run(frame, 50);
        T_sum += instantaneous_T(frame);
        ++samples;
    }
    const double T_mean = T_sum / samples;
    std::cerr << "nvt equilibration <T> = " << T_mean << " K\n";
    BM_NEAR(T_mean, T_target, 0.10 * T_target)
        << "NoseHoover NVT mean temperature off the LAMMPS setpoint";
}

BENCH_MAIN("MD thermostats — LAMMPS API alignment");
