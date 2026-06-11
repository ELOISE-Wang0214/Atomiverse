/// @file benchmarks/md/bench_setup_lammps_align.cpp
/// @brief LAMMPS-API-alignment validation benchmark — simulation-setup
///        components. Ported from the gtest unit-test suite; this is a
///        third-party-comparison validation benchmark, not a ctest gate.
///
/// Verifies that `cpu::md::NeighborList` and `cpu::md::InitVelocity` match
/// the semantics of LAMMPS `neighbor` / `neigh_modify` and `velocity create`
/// (`lammps/src/neighbor.cpp`, `lammps/src/velocity.cpp`):
///
///  - `NeighborList` exposes a `check` flag (LAMMPS `neigh_modify check`).
///    `check no` (`dist_check == 0` in `Neighbor::decide()`) makes the
///    rebuild fire on *every* `every`/`delay`-eligible step regardless of
///    displacement — so a slow-moving system rebuilds more often than with
///    `check yes`.
///  - `InitVelocity` rescales the drawn velocities so the *measured*
///    kinetic temperature equals the requested `T` exactly after COM
///    removal — the defining behaviour of `Velocity::create()`
///    (`rescale(t, t_desired)` at velocity.cpp:368).
///  - `InitVelocity` supports both `dist` modes (`uniform` / `gaussian`).

#include "bench_check.h"

#include <cmath>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/md/velocity_verlet.h"

using namespace ATV;

namespace {

/// 3×3×3 FCC argon crystal at a = 5.26 Å → 108 atoms.
cpu::Frame make_fcc_argon() {
    constexpr double a = 5.26;
    constexpr int n_cell = 3;
    constexpr double box_len = n_cell * a;

    const double basis[4][3] = {
        {0.0, 0.0, 0.0},
        {0.0, 0.5 * a, 0.5 * a},
        {0.5 * a, 0.0, 0.5 * a},
        {0.5 * a, 0.5 * a, 0.0},
    };

    std::vector<double> x, y, z;
    for (int ix = 0; ix < n_cell; ++ix)
        for (int iy = 0; iy < n_cell; ++iy)
            for (int iz = 0; iz < n_cell; ++iz)
                for (auto& b : basis) {
                    x.push_back(ix * a + b[0]);
                    y.push_back(iy * a + b[1]);
                    z.push_back(iz * a + b[2]);
                }

    const int n = static_cast<int>(x.size());
    cpu::Frame frame(n);
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = x[i];
        frame.atoms.y[i] = y[i];
        frame.atoms.z[i] = z[i];
        frame.atoms.mass[i] = 39.948; // Ar amu
        frame.atoms.type[i] = 0;
    }
    frame.box = Mat3::cubic(box_len);
    return frame;
}

/// Measured kinetic temperature with `3N − 3` DOF (COM removed), matching
/// LAMMPS `compute temp` with `extra_dof = dimension`.
double measured_temperature(const cpu::Frame& f) {
    const int n = f.atoms.size();
    double ke2 = 0.0; // 2 * KE / kMvv2e  ==  sum(m v^2)
    for (int i = 0; i < n; ++i) {
        const double v2 = f.atoms.vx[i] * f.atoms.vx[i]
                          + f.atoms.vy[i] * f.atoms.vy[i]
                          + f.atoms.vz[i] * f.atoms.vz[i];
        ke2 += f.atoms.mass[i] * v2;
    }
    const double ke = 0.5 * kMvv2e * ke2;
    const int ndof = 3 * n - 3;
    return 2.0 * ke / (ndof * kBoltzmann);
}

} // namespace

// ── NeighborList: check no vs check yes ───────────────────────────────

BENCH_CASE(SetupLammpsAlign_CheckNoRebuildsEveryEligibleStep) {
    // A slow-moving system with a large skin: under `check yes` the
    // displacement gate (max_disp*2 > skin) essentially never fires, so
    // only the initial build() counts. Under `check no` the rebuild fires
    // on every `every`/`delay`-eligible step regardless of displacement.
    constexpr double cutoff = 6.8;
    constexpr double skin =
        10.0; // huge: skin/2 = 5 Å, atoms never move that far
    constexpr int every = 1;
    constexpr int delay = 0;
    constexpr int n_steps = 50;

    // --- check yes (default) ---
    int rebuilds_yes = 0;
    {
        auto frame = make_fcc_argon();
        cpu::md::NeighborList<> nl{cutoff, skin, every, delay, /*check=*/true};
        cpu::md::MD sim;
        sim.set_integrator(cpu::md::VelocityVerlet{0.5});
        sim.set_pair(cpu::PairLJCut{0.0104, 3.4, cutoff});
        sim.set_neighbor(nl);
        sim.build();
        sim.run(frame, n_steps);
        rebuilds_yes = nl.n_rebuilds();
    }

    // --- check no ---
    int rebuilds_no = 0;
    {
        auto frame = make_fcc_argon();
        cpu::md::NeighborList<> nl{cutoff, skin, every, delay, /*check=*/false};
        cpu::md::MD sim;
        sim.set_integrator(cpu::md::VelocityVerlet{0.5});
        sim.set_pair(cpu::PairLJCut{0.0104, 3.4, cutoff});
        sim.set_neighbor(nl);
        sim.build();
        sim.run(frame, n_steps);
        rebuilds_no = nl.n_rebuilds();
    }

    // check yes: enormous skin + zero velocities → only the initial build.
    BM_EQ(rebuilds_yes, 1)
        << "check yes with enormous skin should rebuild only at build()";

    // check no: rebuild fires on every every/delay-eligible step.
    // every=1, delay=0 → eligible at every step → 1 (initial) + n_steps.
    BM_EQ(rebuilds_no, 1 + n_steps)
        << "check no with every=1 delay=0 must rebuild every step";

    // The binding criterion: check no rebuilds strictly more often.
    BM_GT(rebuilds_no, rebuilds_yes)
        << "check no must rebuild more than check yes for a slow system";
}

BENCH_CASE(SetupLammpsAlign_CheckNoHonoursEveryAndDelayGates) {
    // check no still respects the every/delay gates: `decide()` increments
    // ago then tests `ago >= delay && ago % every == 0` BEFORE the
    // dist_check==0 short-circuit. With every=5, delay=0 the rebuild fires
    // only on steps 5, 10, 15, ...
    constexpr double cutoff = 6.8;
    constexpr double skin = 10.0;
    constexpr int every = 5;
    constexpr int delay = 0;
    constexpr int n_steps = 50;

    auto frame = make_fcc_argon();
    cpu::md::NeighborList<> nl{cutoff, skin, every, delay, /*check=*/false};
    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{0.5});
    sim.set_pair(cpu::PairLJCut{0.0104, 3.4, cutoff});
    sim.set_neighbor(nl);
    sim.build();
    sim.run(frame, n_steps);

    // 1 (initial) + n_steps/every rebuilds.
    BM_EQ(nl.n_rebuilds(), 1 + n_steps / every)
        << "check no must still honour the every gate";
}

// ── InitVelocity: exact temperature ───────────────────────────────────

/// Build an `MD` driver with a Force-capable pair + neighbor list plus
/// the given `InitVelocity`, build it, and return. `MD::build()` requires
/// a Force-capable component; the LJ pair satisfies that without affecting
/// the build-time velocity initialization.
void build_with_init(cpu::md::MD& sim, cpu::md::InitVelocity init) {
    constexpr double cutoff = 6.8;
    sim.set_integrator(cpu::md::VelocityVerlet{0.5});
    sim.set_pair(cpu::PairLJCut{0.0104, 3.4, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff});
    sim.add(init);
    sim.build();
}

BENCH_CASE(SetupLammpsAlign_InitVelocityTemperatureExactGaussian) {
    constexpr double T = 300.0;
    auto frame = make_fcc_argon();

    cpu::md::MD sim;
    build_with_init(sim,
                    cpu::md::InitVelocity{T, /*seed=*/12345u,
                                          cpu::md::VelocityDist::gaussian});
    (void)frame;
    sim.run(frame, 0);

    const double T_meas = measured_temperature(frame);
    BM_NEAR(T_meas, T, 1e-6)
        << "Gaussian InitVelocity must rescale measured T to exactly " << T
        << " K (LAMMPS velocity create), got " << T_meas;
}

BENCH_CASE(SetupLammpsAlign_InitVelocityTemperatureExactUniform) {
    constexpr double T = 150.0;
    auto frame = make_fcc_argon();

    cpu::md::MD sim;
    build_with_init(sim, cpu::md::InitVelocity{T, /*seed=*/777u,
                                               cpu::md::VelocityDist::uniform});
    sim.run(frame, 0);

    const double T_meas = measured_temperature(frame);
    BM_NEAR(T_meas, T, 1e-6)
        << "Uniform InitVelocity must rescale measured T to exactly " << T
        << " K, got " << T_meas;
}

BENCH_CASE(SetupLammpsAlign_InitVelocityDefaultDistGaussian) {
    // Atomiverse keeps Gaussian as its default `dist` (documented
    // divergence from LAMMPS `dist uniform`): the two-arg ctor must still
    // produce an exact temperature.
    constexpr double T = 250.0;
    auto frame = make_fcc_argon();

    cpu::md::MD sim;
    build_with_init(sim, cpu::md::InitVelocity{T, /*seed=*/42u});
    sim.run(frame, 0);

    const double T_meas = measured_temperature(frame);
    BM_NEAR(T_meas, T, 1e-6)
        << "default-dist InitVelocity must rescale to exact T";
}

// ── InitVelocity: COM momentum removal ────────────────────────────────

BENCH_CASE(SetupLammpsAlign_InitVelocityZeroesLinearMomentum) {
    constexpr double T = 300.0;
    auto frame = make_fcc_argon();

    cpu::md::MD sim;
    build_with_init(sim,
                    cpu::md::InitVelocity{T, /*seed=*/9001u,
                                          cpu::md::VelocityDist::gaussian});
    sim.run(frame, 0);

    const int n = frame.atoms.size();
    double p[3] = {0, 0, 0}, mtot = 0;
    for (int i = 0; i < n; ++i) {
        const double m = frame.atoms.mass[i];
        p[0] += m * frame.atoms.vx[i];
        p[1] += m * frame.atoms.vy[i];
        p[2] += m * frame.atoms.vz[i];
        mtot += m;
    }
    // Net COM velocity must be zero (LAMMPS `mom yes` default).
    BM_NEAR(p[0] / mtot, 0.0, 1e-10);
    BM_NEAR(p[1] / mtot, 0.0, 1e-10);
    BM_NEAR(p[2] / mtot, 0.0, 1e-10);
}

// ── InitVelocity: dist modes produce distinct draws ───────────────────

BENCH_CASE(SetupLammpsAlign_InitVelocityDistModesDiffer) {
    // Same seed + same T but different `dist`: the per-atom velocities
    // must differ (different underlying draws), while both still hit the
    // exact target temperature.
    constexpr double T = 300.0;
    auto frame_g = make_fcc_argon();
    auto frame_u = make_fcc_argon();

    {
        cpu::md::MD sim;
        build_with_init(sim,
                        cpu::md::InitVelocity{T, /*seed=*/555u,
                                              cpu::md::VelocityDist::gaussian});
        sim.run(frame_g, 0);
    }
    {
        cpu::md::MD sim;
        build_with_init(sim,
                        cpu::md::InitVelocity{T, /*seed=*/555u,
                                              cpu::md::VelocityDist::uniform});
        sim.run(frame_u, 0);
    }

    bool any_diff = false;
    for (int i = 0; i < frame_g.atoms.size(); ++i) {
        if (std::abs(frame_g.atoms.vx[i] - frame_u.atoms.vx[i]) > 1e-9) {
            any_diff = true;
            break;
        }
    }
    BM_TRUE(any_diff)
        << "gaussian and uniform dist must draw distinct velocities";

    BM_NEAR(measured_temperature(frame_g), T, 1e-6);
    BM_NEAR(measured_temperature(frame_u), T, 1e-6);
}

BENCH_MAIN("MD simulation-setup — LAMMPS API alignment");
