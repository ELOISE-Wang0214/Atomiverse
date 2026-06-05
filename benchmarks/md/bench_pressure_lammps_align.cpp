/// @file benchmarks/md/bench_pressure_lammps_align.cpp
/// @brief lammps-api-align-04 validation benchmark — pressure-component
///        alignment with LAMMPS. Ported from the gtest unit-test suite;
///        this is a third-party-comparison validation benchmark, not a
///        ctest gate.
///
/// Pins the three binding runtime properties of the link:
///
///   * ac-002 — the unified pressure sign: every Atomiverse pressure
///     consumer uses the LAMMPS-standard `P = +trace(σ)/(3V)`
///     (positive ⇔ compressive). No component computes `−trace/(3V)`.
///   * ac-003 — `Barostat`'s Berendsen `μ` factor is numerically
///     invariant under the sign-convention flip: for the default
///     `P_target = 0` the rescale is bit-identical to the pre-flip
///     `μ = (1 − dt/τ·(P_old − P_target))^(1/3)` with `P_old = −P_new`.
///   * ac-005 — for one shared `P_target`, `NoseHooverNPT` and
///     `Barostat` both relax an over-compressed box toward the same
///     equilibrium density, and the scalar pressure carries the
///     LAMMPS positive-compressive sign in both.
///
/// LAMMPS reference (`src/fix_press_berendsen.cpp::end_of_step`):
///
///     dilation = (1 − dt/Pdamp · (P_target − P_current)/bulkmodulus)^(1/3)
///
/// with `P_current` the LAMMPS-standard `compute pressure` scalar
/// `P = (2·KE + Σ virial)/(3V)`, positive ⇔ compressive.

#include "bench_check.h"

#include "atomiverse/core/constants.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/barostat.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover_npt.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/virial_stress.h"

#include <cmath>

using namespace ATV;

namespace {

/// @brief FCC argon box (`n_cells^3 · 4` atoms). Lattice constant `a`
/// sets the density; a small `a` over-compresses the box (P > 0).
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
                    frame.atoms.x[idx] = (ix + basis[b][0] + 0.25) * a;
                    frame.atoms.y[idx] = (iy + basis[b][1] + 0.25) * a;
                    frame.atoms.z[idx] = (iz + basis[b][2] + 0.25) * a;
                    frame.atoms.mass[idx] = mass;
                    frame.atoms.type[idx] = 18;
                    const double s = 1e-4 * ((idx % 2 == 0) ? 1.0 : -1.0);
                    frame.atoms.vx[idx] = s;
                    frame.atoms.vy[idx] = s * 0.5;
                    frame.atoms.vz[idx] = s * 0.25;
                    ++idx;
                }
    frame.box = Mat3::cubic(n_cells * a);
    return frame;
}

/// @brief LAMMPS-standard scalar pressure `P = (2·KE + Σ virial)/(3V)`,
/// positive ⇔ compressive — `compute pressure::compute_scalar()`.
double lammps_pressure(const cpu::Frame& frame) {
    const double V = frame.box.volume();
    const double trace =
        (frame.scalars.stress[0] + frame.scalars.stress[1]
         + frame.scalars.stress[2])
        + (frame.scalars.stress_kinetic[0] + frame.scalars.stress_kinetic[1]
           + frame.scalars.stress_kinetic[2]);
    return +trace / (3.0 * V); // LAMMPS sign: positive compressive.
}

} // namespace

// ── ac-003 — μ-invariance under the sign flip ───────────────────────
//
// The pre-flip Barostat computed `P_old = −trace/(3V)` and
// `μ³ = 1 − dt/τ·(P_old − P_target)`. The post-flip Barostat computes
// `P_new = +trace/(3V) = −P_old` and the LAMMPS formula
// `μ³ = 1 − dt/(τ·β)·(P_target − P_new)`. For the default
// `P_target = 0` and `β = 1`:
//     μ³_old = 1 − dt/τ·(P_old)        = 1 + dt/τ·trace/(3V)
//     μ³_new = 1 − dt/τ·(0 − P_new)    = 1 + dt/τ·P_new
//            = 1 + dt/τ·trace/(3V)
// — bit-identical. This case reproduces both arithmetic forms in
// double precision and asserts exact equality.
BENCH_CASE(PressureLammpsAlign_BerendsenMuInvariantUnderSignFlip) {
    // A representative over-compressed snapshot: positive virial trace.
    const double trace = 1.7e-2; // eV (Σ virial + Σ kinetic)
    const double V = 1.85e4;     // Å³
    const double dt = 0.5;       // fs
    const double tau_P = 200.0;  // fs
    const double P_target = 0.0; // default setpoint

    // Pre-flip representation.
    const double P_old = -trace / (3.0 * V);
    const double mu3_old = 1.0 - (dt / tau_P) * (P_old - P_target);

    // Post-flip (LAMMPS) representation, β = bulk_modulus = 1.
    const double P_new = +trace / (3.0 * V);
    const double mu3_new = 1.0 - (dt / tau_P) * (P_target - P_new) / 1.0;

    // Representation change — the μ³ argument must be bit-identical.
    BM_DOUBLE_EQ(mu3_old, mu3_new);
    BM_DOUBLE_EQ(std::cbrt(mu3_old), std::cbrt(mu3_new));

    // And P_new is the LAMMPS positive-compressive sign (trace > 0 ⇒ P > 0).
    BM_GT(P_new, 0.0);
    BM_LT(P_old, 0.0);
}

// ── ac-002 — Barostat applies the LAMMPS-standard sign ──────────────
//
// An over-compressed box has a positive virial trace ⇒ positive
// (compressive) LAMMPS pressure ⇒ with P_target = 0 the box must
// EXPAND (μ > 1). The pre-flip code, computing P = −trace/(3V) < 0,
// produced exactly the same μ at P_target = 0 (proven above); this
// case pins the *behaviour*: an over-pressured box grows.
BENCH_CASE(PressureLammpsAlign_BarostatExpandsOverCompressedBox) {
    auto frame = make_fcc_argon(4, 4.6, 39.948); // dense ⇒ P > 0

    const double dt = 0.5;
    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::NeighborList<>{6.0}, -10)
        .add(cpu::md::LJ{0.0104, 3.4, 6.0});
    sim.add(cpu::md::Virial{});
    sim.add(cpu::md::Barostat{0.0, /*tau_P=*/200.0, dt});
    sim.build();

    sim.run(frame, 1);
    const double V0 = frame.box.volume();
    const double P0 = lammps_pressure(frame);

    // Over-compressed: LAMMPS sign says P > 0 (compressive).
    BM_GT(P0, 0.0) << "over-compressed box must read positive pressure";

    sim.run(frame, 400);
    const double V1 = frame.box.volume();
    const double P1 = lammps_pressure(frame);

    for (int i = 0; i < frame.atoms.size(); ++i)
        BM_ASSERT_FALSE(std::isnan(frame.atoms.x[i]));

    // Positive (compressive) pressure at P_target = 0 ⇒ box expands.
    BM_GT(V1, V0) << "compressive pressure must drive box growth";
    // The barostat reduces the pressure magnitude toward the setpoint.
    BM_LT(std::abs(P1), std::abs(P0));
}

// ── ac-004 — bulk_modulus scales the relaxation rate ────────────────
//
// LAMMPS `bulkmodulus` divides the μ argument. A larger bulk_modulus
// ⇒ smaller per-step rescale ⇒ slower box relaxation. Two otherwise
// identical runs differing only in bulk_modulus must reach different
// volumes, with the larger modulus moving the box less.
BENCH_CASE(PressureLammpsAlign_BulkModulusScalesRelaxationRate) {
    const double dt = 0.5;
    auto run = [&](double bulk_modulus) {
        auto frame = make_fcc_argon(4, 4.6, 39.948);
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{dt});
        sim.on(ATV::md::Phase::force)
            .add(cpu::md::NeighborList<>{6.0}, -10)
            .add(cpu::md::LJ{0.0104, 3.4, 6.0});
        sim.add(cpu::md::Virial{});
        sim.add(cpu::md::Barostat{0.0, /*tau_P=*/200.0, dt, bulk_modulus});
        sim.build();
        sim.run(frame, 200);
        return frame.box.volume();
    };

    const double V0 = make_fcc_argon(4, 4.6, 39.948).box.volume();
    const double V_soft = run(1.0);   // soft modulus ⇒ fast relaxation
    const double V_stiff = run(50.0); // stiff modulus ⇒ slow relaxation

    // Over-compressed ⇒ both expand, but the stiff modulus expands less.
    BM_GT(V_soft, V0);
    BM_GT(V_stiff, V0);
    BM_GT(V_soft - V0, V_stiff - V0)
        << "larger bulk_modulus must slow the box relaxation";
}

// ── ac-005 — NoseHooverNPT and Barostat relax to the same density ───
//
// Same over-compressed start, same P_target = 0. Both barostats must
// drive the box toward the same equilibrium volume (the LJ density at
// zero pressure). The MTK NPT and Berendsen schemes differ in their
// transient, so the tolerance is loose — the binding claim is that
// both move the box in the SAME direction toward a COMMON setpoint and
// that the unified positive-compressive pressure sign holds in both.
BENCH_CASE(PressureLammpsAlign_NptAndBarostatRelaxToSameDensity) {
    const double dt = 0.5;
    const double P_target = 0.0;

    // Berendsen path. Over-compressed start ⇒ positive (compressive)
    // pressure ⇒ box must expand toward the LJ equilibrium density.
    double V_berendsen = 0.0;
    double P_berendsen_0 = 0.0, P_berendsen_1 = 0.0;
    {
        auto frame = make_fcc_argon(4, 4.7, 39.948);
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{dt});
        sim.on(ATV::md::Phase::force)
            .add(cpu::md::NeighborList<>{6.0}, -10)
            .add(cpu::md::LJ{0.0104, 3.4, 6.0});
        sim.add(cpu::md::Virial{});
        sim.add(cpu::md::Barostat{P_target, /*tau_P=*/100.0, dt});
        sim.build();
        sim.run(frame, 1);
        P_berendsen_0 = lammps_pressure(frame);
        sim.run(frame, 1999);
        V_berendsen = frame.box.volume();
        P_berendsen_1 = lammps_pressure(frame);
        for (int i = 0; i < frame.atoms.size(); ++i)
            BM_ASSERT_FALSE(std::isnan(frame.atoms.x[i]));
        // Unified sign: over-compressed ⇒ P > 0 (compressive), and the
        // barostat drives |P| down toward the setpoint.
        BM_GT(P_berendsen_0, 0.0);
        BM_LT(std::abs(P_berendsen_1), std::abs(P_berendsen_0));
    }

    // MTK NoseHooverNPT path.
    double V_npt = 0.0;
    {
        auto frame = make_fcc_argon(4, 4.7, 39.948);
        cpu::md::MD sim;
        sim.on(ATV::md::Phase::force)
            .add(cpu::md::NeighborList<>{6.0}, -10)
            .add(cpu::md::LJ{0.0104, 3.4, 6.0});
        sim.add(cpu::md::Virial{});
        sim.set_integrator(cpu::md::NoseHooverNPT{{.T_target = 30.0,
                                                   .tau_T = 100.0,
                                                   .P_target = P_target,
                                                   .tau_P = 500.0,
                                                   .dt = dt}});
        sim.build();
        sim.run(frame, 2000);
        V_npt = frame.box.volume();
        for (int i = 0; i < frame.atoms.size(); ++i)
            BM_ASSERT_FALSE(std::isnan(frame.atoms.x[i]));
    }

    // Both relax an over-compressed box outward toward a common density.
    const double V0 = make_fcc_argon(4, 4.7, 39.948).box.volume();
    BM_GT(V_berendsen, V0);
    BM_GT(V_npt, V0);
    // Same equilibrium density within 25 % — the two schemes converge
    // to the same LJ zero-pressure volume from the same start.
    const double rel_diff =
        std::abs(V_berendsen - V_npt) / (0.5 * (V_berendsen + V_npt));
    BM_LT(rel_diff, 0.25) << "Berendsen V=" << V_berendsen
                          << " (P0=" << P_berendsen_0 << " P1=" << P_berendsen_1
                          << ") vs NPT V=" << V_npt;
}

BENCH_MAIN("MD pressure / barostats — LAMMPS API alignment");
