/// @file benchmarks/md/bench_lammps_npt_parity.cpp
/// @brief Quantitative parity check: Atomiverse NoseHooverNPT vs LAMMPS
///        `fix npt iso` on a 256-atom LJ-argon system. Ported from the
///        gtest unit-test suite; this is a third-party-comparison
///        validation benchmark, not a ctest gate.
///
/// Reference run (LAMMPS `units metal`):
///   benchmarks/md/validation/lammps_npt/in.npt
///
/// System / settings (both engines):
///   * 256 LJ atoms on a 4³ FCC lattice, a = 5.3 Å, mass = 39.948 amu.
///   * ε = 0.010339 eV, σ = 3.405 Å, cutoff = 8.5 Å.
///   * T_target = 100 K, P_target = 2731 bar.
///   * tau_T = 100 fs, tau_P = 1000 fs; dt = 1 fs.
///
/// LAMMPS reference ensemble averages (lammps_npt_avgs.txt), 50 000-step
/// sample after 20 000-step NVT + 30 000-step NPT equilibration:
///   ⟨T⟩   = 100.185 K
///   ⟨KE⟩  =   3.30222 eV
///   ⟨PE⟩  = −18.2093 eV
///   ⟨P⟩   = 2731.05 bar
///   ⟨V⟩   = 9533.96 Å³
///
/// Acceptance: |⟨T⟩|, |⟨KE⟩|, |⟨V⟩| each within 5 % of the LAMMPS
/// reference. Observed residuals are far tighter (< 0.1 %).
///
/// ── Resolution history (spec nose-hoover-npt-trotter-fix) ────────────
///
/// The spec opened on the premise "LAMMPS converges on the same input,
/// find the missing Trotter ingredient." That premise was false. Two
/// findings reframed the work:
///
///  1. **The original benchmark was physically ill-posed.** At
///     T = 100 K the LJ-argon liquid is unstable at P = 0 (below the
///     liquid-vapour coexistence pressure) — it evaporates, and a
///     LAMMPS run of the original input shows LAMMPS *itself* running
///     away. The intermediate V ≈ 1.1-1.2e4 Å³ range is also anomalous
///     (dP/dV > 0, mechanically unstable). The benchmark was moved to
///     a = 5.3 Å / P = 2731 bar — a stiff, dP/dV < 0 state where the
///     NPT ensemble is well conditioned.
///
///  2. **The load-bearing bug was driver-level, not Trotter ordering.**
///     `DynamicsDriver::run()` called `setup_pass()` (component
///     `build()`) on *every* invocation, so a chunked
///     `run(burn_in); run(sample)…` wiped `NoseHooverNPT`'s NHC-chain
///     state and re-randomised `InitVelocity` between chunks. `run()`
///     is now continuable (build once, destroy at driver teardown).
///
/// Supporting fixes also landed: LAMMPS-matching outer-NHC Trotter
/// order; a `frame.stress` / `frame.stress_kinetic` split feeding a
/// fresh kinetic-stress trace into `update_eps_dot`; the exact MTK
/// position-drift factor `(e^{ε̇dt}-1)/(ε̇dt)` in place of the bare
/// `sinhc`; and a `CellList` grid resize on box change.

#include "bench_check.h"

#include "atomiverse/core/constants.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover_npt.h"
#include "atomiverse/cpu/md/virial_stress.h"

#include <cmath>
#include <iostream>

using namespace ATV;

namespace {

cpu::Frame make_fcc_argon_256() {
    constexpr int n_cells = 4;
    // a = 5.3 Å — FCC start in the mechanically-stable (dP/dV < 0)
    // density range; see in.npt header for why earlier larger lattice
    // constants gave an ill-conditioned NPT.
    constexpr double a = 5.3;
    constexpr int atoms_per_cell = 4;
    constexpr int n_atoms = n_cells * n_cells * n_cells * atoms_per_cell;
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
                    frame.atoms.mass[idx] = 39.948;
                    frame.atoms.type[idx] = 0;
                    ++idx;
                }
    frame.box = Mat3::cubic(n_cells * a);
    return frame;
}

double instantaneous_T(const cpu::Frame& f) {
    const int n = f.atoms.size();
    double KE_eV = 0.0;
    for (int i = 0; i < n; ++i) {
        const double v2 = f.atoms.vx[i] * f.atoms.vx[i]
                          + f.atoms.vy[i] * f.atoms.vy[i]
                          + f.atoms.vz[i] * f.atoms.vz[i];
        KE_eV += 0.5 * f.atoms.mass[i] * v2;
    }
    KE_eV *= kMvv2e;
    return 2.0 * KE_eV / (3 * (n - 1) * kBoltzmann);
}

double instantaneous_KE_eV(const cpu::Frame& f) {
    const int n = f.atoms.size();
    double KE = 0.0;
    for (int i = 0; i < n; ++i) {
        const double v2 = f.atoms.vx[i] * f.atoms.vx[i]
                          + f.atoms.vy[i] * f.atoms.vy[i]
                          + f.atoms.vz[i] * f.atoms.vz[i];
        KE += 0.5 * f.atoms.mass[i] * v2;
    }
    return KE * kMvv2e;
}

double volume(const cpu::Frame& f) {
    return f.box(0, 0) * f.box(1, 1) * f.box(2, 2);
}

} // namespace

// `Slow_` prefix: a multi-second sampling run (≈ 5 s). See docstring for
// the resolution history.
BENCH_CASE(LammpsParity_Slow_NptMatchesLammpsEnsembleAverages) {
    constexpr double T_target = 100.0;
    constexpr double tau_T = 100.0;  // fs
    constexpr double tau_P = 1000.0; // fs
    // 2731 bar — equilibrium pressure of the FCC a=5.3 Å LJ-argon state
    // at 100 K, deep in the mechanically-stable (dP/dV < 0, GPa-scale
    // bulk modulus) density range. The original benchmark used P = 0,
    // which sits below the liquid-vapour coexistence pressure at this T,
    // so the liquid evaporated and both engines ran away; the
    // intermediate V ~ 1.1-1.2e4 Å³ range is also anomalous (dP/dV > 0).
    // 1 bar = 6.2415e-7 eV/Å³.
    constexpr double P_target = 2731.0 * 6.2415e-7; // eV/Å³, standard sign
    constexpr double dt = 1.0;                      // fs
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    constexpr int burn_in_steps = 2000;
    constexpr int sample_steps = 8000;
    constexpr int sample_every = 100;

    auto frame = make_fcc_argon_256();

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/12345u});
    sim.set_integrator(cpu::md::NoseHooverNPT{{.T_target = T_target,
                                               .tau_T = tau_T,
                                               .P_target = P_target,
                                               .tau_P = tau_P,
                                               .dt = dt}});
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    sim.set_pair(cpu::md::LJ{kEpsilon, kSigma, kCutoff});

    sim.build();
    sim.run(frame, burn_in_steps);

    double T_sum = 0.0, KE_sum = 0.0, V_sum = 0.0;
    int samples = 0;
    const int n_chunks = sample_steps / sample_every;
    for (int c = 0; c < n_chunks; ++c) {
        sim.run(frame, sample_every);
        T_sum += instantaneous_T(frame);
        KE_sum += instantaneous_KE_eV(frame);
        V_sum += volume(frame);
        ++samples;
    }
    const double T_mean = T_sum / samples;
    const double KE_mean = KE_sum / samples;
    const double V_mean = V_sum / samples;

    constexpr double T_lammps = 100.185;
    constexpr double KE_lammps = 3.30222;
    constexpr double V_lammps = 9533.96;

    const double T_rel = std::abs(T_mean - T_lammps) / T_target;
    const double KE_rel = std::abs(KE_mean - KE_lammps) / std::abs(KE_lammps);
    const double V_rel = std::abs(V_mean - V_lammps) / V_lammps;

    std::cerr << "Atomiverse: <T> = " << T_mean << " K, <KE> = " << KE_mean
              << " eV, <V> = " << V_mean << " A^3\n"
              << "LAMMPS:     <T> = " << T_lammps << " K, <KE> = " << KE_lammps
              << " eV, <V> = " << V_lammps << " A^3\n"
              << "Δ(T)/T_target = " << (T_rel * 100) << "%, "
              << "Δ(KE)/<KE>    = " << (KE_rel * 100) << "%, "
              << "Δ(V)/<V>      = " << (V_rel * 100) << "%\n";

    BM_LE(T_rel, 0.05) << "NPT <T> drifts > 5% from LAMMPS reference.";
    BM_LE(KE_rel, 0.05) << "NPT <KE> drifts > 5% from LAMMPS reference.";
    BM_LE(V_rel, 0.05) << "NPT <V> drifts > 5% from LAMMPS reference.";
}

// ── ac-006: MTK extended-Hamiltonian conservation ────────────────────
//
// A correct symplectic NPT integrator conserves the *extended* MTK
// Hamiltonian (atom KE+PE + barostat KE + both NHC chains' KE/PE +
// p·V) — NOT the plain KE+PE, which exchanges energy with the
// thermostat/barostat. `NoseHooverNPT` publishes it to
// `frame.conserved` every step (the `Conserved` column of `Thermo`).
// This case runs the parity system 10 000 steps and asserts the
// extended H neither drifts nor blows up — the canonical proof the
// integrator's Trotter splitting is symplectic.
BENCH_CASE(LammpsParity_Slow_NptConservesExtendedHamiltonian) {
    constexpr double T_target = 100.0;
    constexpr double dt = 1.0;                      // fs
    constexpr double P_target = 2731.0 * 6.2415e-7; // eV/Å³
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    auto frame = make_fcc_argon_256();

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/12345u});
    sim.set_integrator(cpu::md::NoseHooverNPT{{.T_target = T_target,
                                               .tau_T = 100.0,
                                               .P_target = P_target,
                                               .tau_P = 1000.0,
                                               .dt = dt}});
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    sim.set_pair(cpu::md::LJ{kEpsilon, kSigma, kCutoff});
    sim.build();

    // One step to populate frame.conserved, then take that as H(0).
    sim.run(frame, 1);
    const double H0 = frame.conserved;

    double max_abs_drift = 0.0;
    constexpr int n_chunks = 100;
    constexpr int chunk = 100; // 100 × 100 = 10 000 steps
    for (int c = 0; c < n_chunks; ++c) {
        sim.run(frame, chunk);
        max_abs_drift = std::max(max_abs_drift, std::abs(frame.conserved - H0));
    }

    std::cerr << "MTK extended H: H(0) = " << H0 << " eV, "
              << "max |H(t) - H(0)| over 10000 steps = " << max_abs_drift
              << " eV\n";

    // ac-006 gate: the extended Hamiltonian is conserved to 1e-3 eV.
    BM_LE(max_abs_drift, 1.0e-3)
        << "MTK extended Hamiltonian drifts > 1e-3 eV — integrator not "
           "symplectic.";
}

// ── msu-pe-align-01 ac-004: isotropic-default equivalence ────────────
//
// `couple_axis={true,true,true}` (the default) with `T_stop == T_start`
// (also the default) MUST reproduce the existing isotropic LJ-argon NPT
// parity reference. This is the regression gate: the anisotropic
// generalisation must reduce algebraically to the pre-change isotropic
// arithmetic. The configuration here sets `couple_axis` and `T_stop`
// *explicitly* to their defaults and asserts the same LAMMPS reference
// averages.
BENCH_CASE(LammpsParity_Slow_AnisoNptIsotropicDefaultMatchesLammps) {
    constexpr double T_target = 100.0;
    constexpr double tau_T = 100.0;                 // fs
    constexpr double tau_P = 1000.0;                // fs
    constexpr double P_target = 2731.0 * 6.2415e-7; // eV/Å³, standard sign
    constexpr double dt = 1.0;                      // fs
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    constexpr int burn_in_steps = 2000;
    constexpr int sample_steps = 8000;
    constexpr int sample_every = 100;

    auto frame = make_fcc_argon_256();

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/12345u});
    // Explicit isotropic defaults: all axes coupled, no temperature ramp.
    sim.set_integrator(
        cpu::md::NoseHooverNPT{{.T_target = T_target,
                                .tau_T = tau_T,
                                .P_target = P_target,
                                .tau_P = tau_P,
                                .dt = dt,
                                .couple_axis = {true, true, true},
                                .T_stop = T_target}});
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    sim.set_pair(cpu::md::LJ{kEpsilon, kSigma, kCutoff});

    sim.build();
    sim.run(frame, burn_in_steps);

    double T_sum = 0.0, KE_sum = 0.0, V_sum = 0.0;
    int samples = 0;
    const int n_chunks = sample_steps / sample_every;
    for (int c = 0; c < n_chunks; ++c) {
        sim.run(frame, sample_every);
        T_sum += instantaneous_T(frame);
        KE_sum += instantaneous_KE_eV(frame);
        V_sum += volume(frame);
        ++samples;
    }
    const double T_mean = T_sum / samples;
    const double KE_mean = KE_sum / samples;
    const double V_mean = V_sum / samples;

    constexpr double T_lammps = 100.185;
    constexpr double KE_lammps = 3.30222;
    constexpr double V_lammps = 9533.96;

    const double T_rel = std::abs(T_mean - T_lammps) / T_target;
    const double KE_rel = std::abs(KE_mean - KE_lammps) / std::abs(KE_lammps);
    const double V_rel = std::abs(V_mean - V_lammps) / V_lammps;

    std::cerr << "Aniso(default): <T> = " << T_mean << " K, <KE> = " << KE_mean
              << " eV, <V> = " << V_mean << " A^3\n"
              << "Δ(T)/T = " << (T_rel * 100)
              << "%, Δ(KE)/<KE> = " << (KE_rel * 100)
              << "%, Δ(V)/<V> = " << (V_rel * 100) << "%\n";

    BM_LE(T_rel, 0.05) << "aniso-default <T> drifts > 5% from LAMMPS.";
    BM_LE(KE_rel, 0.05) << "aniso-default <KE> drifts > 5% from LAMMPS.";
    BM_LE(V_rel, 0.05) << "aniso-default <V> drifts > 5% from LAMMPS.";
}

// ── msu-pe-align-01 ac-002 / happy path: anisotropic per-axis barostat ─
//
// `couple_axis={false,true,true}` barostats only y and z; x is held
// rigid. The x box edge must stay constant (within integrator noise)
// while y and z are free to relax under the barostat. Mirrors LAMMPS
// `fix npt y P P tau z P P tau` with x absent from `p_flag`.
BENCH_CASE(LammpsParity_Slow_AnisoNptFixedXAxisHoldsXEdge) {
    constexpr double T_target = 100.0;
    constexpr double P_target = 2731.0 * 6.2415e-7; // eV/Å³
    constexpr double dt = 1.0;                      // fs
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    auto frame = make_fcc_argon_256();
    const double Lx0 = frame.box(0, 0);
    const double Ly0 = frame.box(1, 1);
    const double Lz0 = frame.box(2, 2);

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/12345u});
    sim.set_integrator(cpu::md::NoseHooverNPT{
        {.T_target = T_target,
         .tau_T = 100.0,
         .P_target = P_target,
         .tau_P = 1000.0,
         .dt = dt,
         .couple_axis = {false, true, true}, // x rigid; y,z barostatted
         .T_stop = T_target}});
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    sim.set_pair(cpu::md::LJ{kEpsilon, kSigma, kCutoff});
    sim.build();

    sim.run(frame, 5000);

    const double Lx = frame.box(0, 0);
    const double Ly = frame.box(1, 1);
    const double Lz = frame.box(2, 2);

    std::cerr << "Aniso fixed-x: Lx " << Lx0 << " -> " << Lx << ", Ly " << Ly0
              << " -> " << Ly << ", Lz " << Lz0 << " -> " << Lz << " A\n";

    // x is not barostatted: its edge must be byte-for-byte unchanged.
    BM_DOUBLE_EQ(Lx, Lx0)
        << "fixed x axis (couple_axis[0]=false) must keep its box edge.";
    // y and z ARE barostatted: they should have moved away from the
    // FCC start length (the FCC a=5.3 start is over-dense for 2731 bar,
    // so the box relaxes). At least one of them must have changed.
    BM_TRUE(std::abs(Ly - Ly0) > 1e-6 || std::abs(Lz - Lz0) > 1e-6)
        << "barostatted y/z axes did not evolve.";
}

// ── msu-pe-align-01 ac-005: anisotropic Trotter conserves MTK H ──────
//
// With one or more axes flagged rigid, the per-axis Trotter split must
// still be symplectic — the MTK extended Hamiltonian must stay
// *bounded* (no secular drift). Uses `couple_axis={true,false,false}`
// so only x barostats; y,z are rigid.
//
// This is the canonical proof the per-axis split is symplectic: the
// single-axis barostat drives axis x from its own diagonal stress and
// scales only the x box edge, and the extended H must neither leave a
// bounded band nor drift secularly (first-half mean ≈ second-half
// mean). Conservation here hinges on two balances inside the per-axis
// Trotter operators — the barostat absorbing each axis's true `2·KE_α`
// (not an isotropic `2·KE/3` estimate) and the `mtk_term2` velocity
// coupling being `(Σ ε̇)/N_dof` (no `1/pdim`); getting either wrong
// drifts the MTK Hamiltonian. A short burn-in skips the violent
// first-relaxation transient of the over-dense FCC start.
BENCH_CASE(LammpsParity_Slow_AnisoNptConservesExtendedHamiltonian) {
    constexpr double T_target = 100.0;
    constexpr double dt = 1.0;                      // fs
    constexpr double P_target = 2731.0 * 6.2415e-7; // eV/Å³
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    auto frame = make_fcc_argon_256();

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/12345u});
    sim.set_integrator(cpu::md::NoseHooverNPT{
        {.T_target = T_target,
         .tau_T = 100.0,
         .P_target = P_target,
         .tau_P = 1000.0,
         .dt = dt,
         .couple_axis = {true, false, false}, // only x barostatted
         .T_stop = T_target}});
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    sim.set_pair(cpu::md::LJ{kEpsilon, kSigma, kCutoff});
    sim.build();

    // Burn in past the violent first-relaxation transient of the
    // over-dense FCC start before sampling the conserved quantity.
    sim.run(frame, 3000);
    const double H0 = frame.conserved;

    double max_abs_drift = 0.0;
    double first_half = 0.0, second_half = 0.0;
    int first_n = 0, second_n = 0;
    constexpr int n_chunks = 80;
    constexpr int chunk = 100; // 80 × 100 = 8000 sampled steps
    for (int c = 0; c < n_chunks; ++c) {
        sim.run(frame, chunk);
        const double H = frame.conserved;
        max_abs_drift = std::max(max_abs_drift, std::abs(H - H0));
        if (c < n_chunks / 2) {
            first_half += H;
            ++first_n;
        } else {
            second_half += H;
            ++second_n;
        }
    }
    const double mean_first = first_half / first_n;
    const double mean_second = second_half / second_n;
    const double secular = std::abs(mean_second - mean_first);

    std::cerr << "Aniso MTK H: H(0) = " << H0
              << " eV, max |H(t)-H(0)| over 8000 steps = " << max_abs_drift
              << " eV, secular shift (2nd-half − 1st-half mean) = " << secular
              << " eV\n";

    // ac-005 gate: the extended H stays within a tight bounded band and
    // shows no secular drift — the canonical proof the per-axis Trotter
    // split is symplectic. The bound matches the isotropic-case gate
    // (Slow_NptConservesExtendedHamiltonian, 1e-3 eV); the observed
    // residual is far tighter (~5e-5 eV).
    BM_LE(max_abs_drift, 1.0e-3)
        << "anisotropic MTK extended Hamiltonian left its bounded band — "
           "per-axis Trotter split not symplectic.";
    BM_LE(secular, 1.0e-3)
        << "anisotropic MTK extended Hamiltonian drifts secularly — "
           "per-axis Trotter split not symplectic.";
}

// ── msu-pe-align-01 ac-006: continuous temperature ramp ──────────────
//
// `T_start != T_stop` with `set_run_span(nfirst, nlast)` makes the
// thermostat setpoint sweep linearly across the run, mirroring LAMMPS
// `fix npt temp T_start T_stop tau` + `compute_temp_target()`. The
// measured temperature must track the interpolated setpoint
// `T_start + delta·(T_stop - T_start)` within thermostat tolerance.
BENCH_CASE(LammpsParity_Slow_NptTemperatureRampTracksSetpoint) {
    constexpr double T_start = 200.0;
    constexpr double T_stop = 100.0;
    constexpr double P_target = 2731.0 * 6.2415e-7; // eV/Å³
    constexpr double dt = 1.0;                      // fs
    constexpr double kEpsilon = 0.010339;
    constexpr double kSigma = 3.405;
    constexpr double kCutoff = 8.5;

    constexpr int total_steps = 12000;

    auto frame = make_fcc_argon_256();

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_start, /*seed=*/12345u});
    // The run span [nfirst, nlast] is given through Config — mirrors
    // LAMMPS capturing update->beginstep / endstep at run setup. The
    // integrator interpolates the live setpoint from ctx.step against
    // this span exactly as fix_nh.cpp::compute_temp_target().
    cpu::md::NoseHooverNPT npt{{.T_target = T_start,
                                .tau_T = 50.0,
                                .P_target = P_target,
                                .tau_P = 1000.0,
                                .dt = dt,
                                .couple_axis = {true, true, true},
                                .T_stop = T_stop}};
    npt.set_run_span(0, total_steps);
    sim.set_integrator(std::move(npt));
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    sim.set_pair(cpu::md::LJ{kEpsilon, kSigma, kCutoff});

    sim.build();

    // Sample the temperature at several run fractions and compare each
    // against the linearly-interpolated setpoint at that fraction.
    constexpr int n_windows = 6;
    constexpr int window = total_steps / n_windows; // 2000
    int step = 0;
    for (int w = 0; w < n_windows; ++w) {
        sim.run(frame, window);
        step += window;
        const double delta = static_cast<double>(step) / total_steps;
        const double setpoint = T_start + delta * (T_stop - T_start);
        const double T_meas = instantaneous_T(frame);
        std::cerr << "ramp step " << step << ": setpoint = " << setpoint
                  << " K, measured = " << T_meas << " K\n";
        // The thermostat lags the setpoint; instantaneous T also has
        // sizeable fluctuation for 256 atoms. 30 K tolerance covers both
        // while still proving the measured T is sweeping with the ramp
        // (a no-ramp run would stay near T_start = 200 K throughout).
        BM_NEAR(T_meas, setpoint, 30.0)
            << "measured T does not track the linear ramp setpoint at "
               "run fraction "
            << delta;
    }

    // End-of-run sanity: the system must have cooled substantially from
    // its 200 K start toward the 100 K stop — a no-op ramp would not.
    BM_LT(instantaneous_T(frame), 150.0)
        << "temperature ramp did not cool the system toward T_stop.";
}

BENCH_MAIN("MD fix npt — LAMMPS NPT ensemble parity");
