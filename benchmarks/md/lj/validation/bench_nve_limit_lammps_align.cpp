/// @file benchmarks/md/bench_nve_limit_lammps_align.cpp
/// @brief LAMMPS API-alignment validation benchmark for the NVE family
///        (spec lammps-api-align-02-nve): `cpu::md::VelocityVerlet` ↔
///        LAMMPS `fix nve` and `cpu::md::NVELimit` ↔ LAMMPS
///        `fix nve/limit`. Ported from the gtest unit-test suite; this is
///        a third-party-comparison validation benchmark, not a ctest gate.
///
/// Grounded in the verbatim LAMMPS sources `src/fix_nve.cpp` and
/// `src/fix_nve_limit.cpp`:
///
///   - `fix nve/limit` (`init()`): `vlimitsq = (xlimit/dtv)^2` — `xlimit`
///     is a *maximum distance moved per step*, not a velocity.
///   - `initial_integrate()` / `final_integrate()`: after `v += dtfm*f`,
///     `vsq = v·v` (Euclidean); `if (vsq > vlimitsq) scale =
///     sqrt(vlimitsq/vsq); v *= scale` — the cap rescales `v` so that
///     `|v|·dt <= xlimit`. The cap is applied in BOTH integrate halves.
///   - `fix nve` is plain velocity-Verlet, parameter-free.
///
/// These checks assert the Atomiverse semantics match: `xmax` is a
/// per-step distance cap on the Euclidean speed, engaging in both the
/// initial and final integrate halves; a large `xmax` reduces `NVELimit`
/// to plain `VelocityVerlet`; and `VelocityVerlet` conserves NVE energy.

#include "bench_check.h"

#include <cmath>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/nve_limit.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/potentials/bond_harmonic.h"

using namespace ATV;

namespace {

constexpr int kInitial = static_cast<int>(ATV::md::Phase::initial_integrate);
constexpr int kFinal = static_cast<int>(ATV::md::Phase::final_integrate);

double speed(const cpu::Frame& f, int i) {
    return std::sqrt(f.atoms.vx[i] * f.atoms.vx[i]
                     + f.atoms.vy[i] * f.atoms.vy[i]
                     + f.atoms.vz[i] * f.atoms.vz[i]);
}

// Single free atom at the box centre, unit mass, zero velocity.
cpu::Frame make_single_atom() {
    cpu::Frame frame(1);
    frame.atoms.x[0] = 10.0;
    frame.atoms.y[0] = 10.0;
    frame.atoms.z[0] = 10.0;
    frame.atoms.vx[0] = 0.0;
    frame.atoms.vy[0] = 0.0;
    frame.atoms.vz[0] = 0.0;
    frame.atoms.mass[0] = 1.0;
    frame.atoms.type[0] = 0;
    frame.box = Mat3::cubic(1000.0); // huge box → no PBC wrap interference
    return frame;
}

} // namespace

// ── ac-002: xmax is a per-step DISTANCE cap on the Euclidean speed ────
//
// LAMMPS `fix nve/limit` caps so |v|·dt <= xmax. Inject a force spike in
// the INITIAL half that, unclamped, would accelerate the atom past the
// cap; verify the post-kick Euclidean speed equals exactly xmax/dt and
// the displacement of that step equals exactly xmax.

BENCH_CASE(NveLimitLammpsAlign_InitialHalfDistanceCap) {
    auto frame = make_single_atom();
    const double dt = 0.5;    // fs
    const double xmax = 0.05; // Å  → vmax = xmax/dt = 0.1 Å/fs

    // Force spike (in eV/Å). dtfm = 0.5*dt*kFtm2v/m. With a large enough
    // force the post-kick speed far exceeds vmax, so the cap MUST engage.
    const double F = 1.0e4;
    frame.atoms.fx[0] = 0.6 * F;
    frame.atoms.fy[0] = 0.8 * F; // |f| direction is (0.6, 0.8, 0)
    frame.atoms.fz[0] = 0.0;

    const double x0 = frame.atoms.x[0];
    const double y0 = frame.atoms.y[0];
    const double z0 = frame.atoms.z[0];

    cpu::md::NVELimit limit(dt, xmax);
    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, kInitial, 0};
    limit.launch(lctx);

    const double vmax = xmax / dt;
    BM_NEAR(speed(frame, 0), vmax, 1e-12)
        << "post-kick Euclidean speed must be capped to xmax/dt";

    // Displacement of this single step must be exactly xmax (distance cap).
    const double disp =
        std::sqrt((frame.atoms.x[0] - x0) * (frame.atoms.x[0] - x0)
                  + (frame.atoms.y[0] - y0) * (frame.atoms.y[0] - y0)
                  + (frame.atoms.z[0] - z0) * (frame.atoms.z[0] - z0));
    BM_NEAR(disp, xmax, 1e-12)
        << "per-step displacement must be bounded by xmax (a distance)";

    // Cap preserves direction: v stays parallel to the force spike.
    BM_NEAR(frame.atoms.vx[0] / speed(frame, 0), 0.6, 1e-12);
    BM_NEAR(frame.atoms.vy[0] / speed(frame, 0), 0.8, 1e-12);
}

// ── ac-003: the cap also engages in the FINAL integrate half ─────────
//
// final_integrate() does v += dtfm*f only (no x update). Inject the
// spike here; verify the post-kick speed is still capped to xmax/dt.

BENCH_CASE(NveLimitLammpsAlign_FinalHalfDistanceCap) {
    auto frame = make_single_atom();
    const double dt = 0.5;
    const double xmax = 0.05;

    const double F = 1.0e4;
    frame.atoms.fx[0] = 0.0;
    frame.atoms.fy[0] = 0.0;
    frame.atoms.fz[0] = F; // pure +z spike

    cpu::md::NVELimit limit(dt, xmax);
    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, kFinal, 0};
    limit.launch(lctx);

    const double vmax = xmax / dt;
    BM_NEAR(speed(frame, 0), vmax, 1e-12)
        << "final_integrate half must also cap |v| to xmax/dt";
    // final_integrate must NOT move the atom.
    BM_NEAR(frame.atoms.x[0], 10.0, 1e-12);
    BM_NEAR(frame.atoms.y[0], 10.0, 1e-12);
    BM_NEAR(frame.atoms.z[0], 10.0, 1e-12);
}

// ── ac-002: the cap uses the Euclidean (L2) norm, not per-component ──
//
// A velocity equal in all three components must be rescaled as a vector:
// the resulting |v|_2 equals xmax/dt, while each component is xmax/dt/√3.

BENCH_CASE(NveLimitLammpsAlign_EuclideanNormNotPerComponent) {
    auto frame = make_single_atom();
    const double dt = 0.5;
    const double xmax = 0.05;
    const double vmax = xmax / dt;

    // Pre-existing isotropic velocity, each component already > vmax so a
    // per-component cap would give |v| = vmax*√3 ≠ vmax.
    frame.atoms.vx[0] = 1.0;
    frame.atoms.vy[0] = 1.0;
    frame.atoms.vz[0] = 1.0;
    frame.atoms.fx[0] = 0.0;
    frame.atoms.fy[0] = 0.0;
    frame.atoms.fz[0] = 0.0;

    cpu::md::NVELimit limit(dt, xmax);
    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, kFinal, 0};
    limit.launch(lctx);

    BM_NEAR(speed(frame, 0), vmax, 1e-12)
        << "cap must use Euclidean speed vsq = v·v";
    BM_NEAR(frame.atoms.vx[0], vmax / std::sqrt(3.0), 1e-12);
    BM_NEAR(frame.atoms.vy[0], vmax / std::sqrt(3.0), 1e-12);
    BM_NEAR(frame.atoms.vz[0], vmax / std::sqrt(3.0), 1e-12);
}

// ── ac-002: below-threshold velocities pass through unchanged ─────────

BENCH_CASE(NveLimitLammpsAlign_NoCapBelowThreshold) {
    auto frame = make_single_atom();
    const double dt = 0.5;
    const double xmax = 0.05; // vmax = 0.1 Å/fs

    frame.atoms.vx[0] = 0.03;
    frame.atoms.vy[0] = 0.04; // |v| = 0.05 < vmax = 0.1 → untouched
    frame.atoms.vz[0] = 0.0;

    cpu::md::NVELimit limit(dt, xmax);
    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, kFinal, 0};
    limit.launch(lctx);

    BM_NEAR(frame.atoms.vx[0], 0.03, 1e-14);
    BM_NEAR(frame.atoms.vy[0], 0.04, 1e-14);
    BM_NEAR(frame.atoms.vz[0], 0.0, 1e-14);
}

// ── ac-004: large xmax → NVELimit reduces to plain VelocityVerlet ────
//
// Drive a single harmonic-bonded pair with both integrators. With xmax
// far larger than any per-step displacement, the cap never engages and
// the two trajectories must be bit-for-bit identical.

BENCH_CASE(NveLimitLammpsAlign_LargeXmaxReducesToVelocityVerlet) {
    auto build_pair = []() {
        cpu::Frame f(2);
        f.atoms.x[0] = 9.0;
        f.atoms.y[0] = 10.0;
        f.atoms.z[0] = 10.0;
        f.atoms.x[1] = 11.4;
        f.atoms.y[1] = 10.0;
        f.atoms.z[1] = 10.0; // stretched bond
        for (int i = 0; i < 2; ++i) {
            f.atoms.mass[i] = 12.0;
            f.atoms.type[i] = 0;
            f.atoms.vx[i] = f.atoms.vy[i] = f.atoms.vz[i] = 0.0;
        }
        f.bonds.atomi = xt::zeros<int>({1});
        f.bonds.atomj = xt::zeros<int>({1});
        f.bonds.type = xt::zeros<int>({1});
        f.bonds.atomi[0] = 0;
        f.bonds.atomj[0] = 1;
        f.box = Mat3::cubic(1000.0);
        return f;
    };

    const double dt = 0.5;
    const double K = 20.0;  // eV/Å²
    const double r0 = 1.53; // Å
    const int steps = 200;

    cpu::Frame frame_vv = build_pair();
    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{dt});
        sim.on(ATV::md::Phase::force).add(cpu::BondHarmonic{{K}, {r0}});
        sim.build();
        sim.run(frame_vv, steps);
    }

    cpu::Frame frame_lim = build_pair();
    {
        cpu::md::MD sim;
        // xmax = 100 Å — far beyond any per-step displacement here.
        sim.add(cpu::md::NVELimit{dt, 100.0});
        sim.on(ATV::md::Phase::force).add(cpu::BondHarmonic{{K}, {r0}});
        sim.build();
        sim.run(frame_lim, steps);
    }

    for (int i = 0; i < 2; ++i) {
        BM_NEAR(frame_lim.atoms.x[i], frame_vv.atoms.x[i], 1e-12)
            << "atom " << i << " x diverged: NVELimit with large xmax must "
            << "match plain VelocityVerlet";
        BM_NEAR(frame_lim.atoms.vx[i], frame_vv.atoms.vx[i], 1e-12)
            << "atom " << i << " vx diverged";
    }
}

// ── ac-004: VelocityVerlet conserves NVE total energy ────────────────
//
// A harmonic-bonded pair given a small stretch oscillates; total energy
// (potential + kinetic) must stay constant within symplectic drift.

BENCH_CASE(NveLimitLammpsAlign_VelocityVerletConservesEnergy) {
    cpu::Frame frame(2);
    frame.atoms.x[0] = 9.3;
    frame.atoms.y[0] = 10.0;
    frame.atoms.z[0] = 10.0;
    frame.atoms.x[1] = 11.0;
    frame.atoms.y[1] = 10.0;
    frame.atoms.z[1] = 10.0;
    for (int i = 0; i < 2; ++i) {
        frame.atoms.mass[i] = 12.0;
        frame.atoms.type[i] = 0;
        frame.atoms.vx[i] = frame.atoms.vy[i] = frame.atoms.vz[i] = 0.0;
    }
    frame.bonds.atomi = xt::zeros<int>({1});
    frame.bonds.atomj = xt::zeros<int>({1});
    frame.bonds.type = xt::zeros<int>({1});
    frame.bonds.atomi[0] = 0;
    frame.bonds.atomj[0] = 1;
    frame.box = Mat3::cubic(1000.0);

    // dt = 0.1 fs keeps the velocity-Verlet symplectic energy error
    // (which scales as O(dt^2)) comfortably inside the 1e-4 drift bound
    // for this large-amplitude bond oscillation.
    const double dt = 0.1;
    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.on(ATV::md::Phase::force).add(cpu::BondHarmonic{{20.0}, {1.53}});
    sim.build();

    sim.run(frame, 1);
    const double e0 = frame.scalars.energy + frame.scalars.kinetic;
    if (!(std::abs(e0) > 0.0)) {
        BM_TRUE(false) << "initial total energy must be nonzero";
        return;
    }

    sim.run(frame, 1999);
    const double e1 = frame.scalars.energy + frame.scalars.kinetic;

    const double drift = std::abs(e1 - e0) / std::abs(e0);
    BM_LT(drift, 1e-4) << "VelocityVerlet NVE drift " << drift
                       << " exceeds 1e-4 over "
                       << "2000 steps (E0=" << e0 << ", E1=" << e1 << ")";
}

// ── ac-003: far-from-equilibrium run — per-step displacement bounded ─
//
// Two atoms placed almost on top of each other → a colossal repulsive
// bond restoring force. Without the cap they would fly apart in one
// step; NVELimit must bound every per-step displacement by xmax.

BENCH_CASE(NveLimitLammpsAlign_FarFromEquilibriumDisplacementBounded) {
    cpu::Frame frame(2);
    frame.atoms.x[0] = 10.00;
    frame.atoms.y[0] = 10.0;
    frame.atoms.z[0] = 10.0;
    frame.atoms.x[1] = 10.02;
    frame.atoms.y[1] = 10.0;
    frame.atoms.z[1] = 10.0;
    for (int i = 0; i < 2; ++i) {
        frame.atoms.mass[i] = 1.0;
        frame.atoms.type[i] = 0;
        frame.atoms.vx[i] = frame.atoms.vy[i] = frame.atoms.vz[i] = 0.0;
    }
    frame.bonds.atomi = xt::zeros<int>({1});
    frame.bonds.atomj = xt::zeros<int>({1});
    frame.bonds.type = xt::zeros<int>({1});
    frame.bonds.atomi[0] = 0;
    frame.bonds.atomj[0] = 1;
    frame.box = Mat3::cubic(1000.0);

    const double dt = 0.5;
    const double xmax = 0.1; // Å per step
    const int steps = 100;

    cpu::md::MD sim;
    sim.add(cpu::md::NVELimit{dt, xmax});
    // r0 far from the 0.02 Å start → enormous restoring force.
    sim.on(ATV::md::Phase::force).add(cpu::BondHarmonic{{500.0}, {1.53}});
    sim.build();

    double prev_x0 = frame.atoms.x[0], prev_y0 = frame.atoms.y[0],
           prev_z0 = frame.atoms.z[0];
    double prev_x1 = frame.atoms.x[1], prev_y1 = frame.atoms.y[1],
           prev_z1 = frame.atoms.z[1];

    for (int s = 0; s < steps; ++s) {
        sim.run(frame, 1);
        const double d0 = std::sqrt(
            (frame.atoms.x[0] - prev_x0) * (frame.atoms.x[0] - prev_x0)
            + (frame.atoms.y[0] - prev_y0) * (frame.atoms.y[0] - prev_y0)
            + (frame.atoms.z[0] - prev_z0) * (frame.atoms.z[0] - prev_z0));
        const double d1 = std::sqrt(
            (frame.atoms.x[1] - prev_x1) * (frame.atoms.x[1] - prev_x1)
            + (frame.atoms.y[1] - prev_y1) * (frame.atoms.y[1] - prev_y1)
            + (frame.atoms.z[1] - prev_z1) * (frame.atoms.z[1] - prev_z1));
        // Per-step displacement = |v_capped|·dt ≤ xmax. A small tolerance
        // covers floating-point round-off in the rescale.
        BM_LE(d0, xmax + 1e-9) << "atom 0 step " << s << " moved " << d0;
        BM_LE(d1, xmax + 1e-9) << "atom 1 step " << s << " moved " << d1;
        if (std::isnan(frame.atoms.x[0])) {
            BM_TRUE(false) << "atom 0 blew up at step " << s;
            return;
        }

        prev_x0 = frame.atoms.x[0];
        prev_y0 = frame.atoms.y[0];
        prev_z0 = frame.atoms.z[0];
        prev_x1 = frame.atoms.x[1];
        prev_y1 = frame.atoms.y[1];
        prev_z1 = frame.atoms.z[1];
    }
}

BENCH_MAIN("MD fix nve / nve/limit — LAMMPS API alignment");
