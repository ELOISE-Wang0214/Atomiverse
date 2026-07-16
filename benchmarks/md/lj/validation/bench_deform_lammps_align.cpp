/// @file benchmarks/md/bench_deform_lammps_align.cpp
/// @brief LAMMPS `fix deform` API-alignment validation benchmark for
///        cpu::md::Deform. Ported from the gtest unit-test suite; this is
///        a third-party-comparison validation benchmark, not a ctest gate.
///
/// Audited against lammps/src/fix_deform.cpp (ERATE / TRATE strain modes
/// and the L0-at-fix-issue snapshot). Covers:
///   * Mode::erate engineering-strain trajectory L(t) = L0*(1 + R*t)
///   * Mode::trate true-strain trajectory      L(t) = L0*exp(R*t)
///   * the reference-frame regression: strain referenced to the box edge
///     at Deform activation, NOT the simulation-initial box -- the exact
///     failure mode the MSU 10000-atom polyethylene run hit.
///
/// `Deform` is driven directly through its component lifecycle
/// (`build()` + repeated `launch()`) rather than through a full
/// `cpu::md::MD` driver. This isolates the box-deformation semantics
/// under test from the unrelated `cpu::md::MD` requirement that a
/// Force-capable evaluator be present.

#include "bench_check.h"

#include <cmath>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/deform.h"

using namespace ATV;

namespace {

/// Two atoms in a cubic box. Only Deform moves the box / coordinates.
cpu::Frame make_two_atom_box(double box_len) {
    cpu::Frame frame(2);
    frame.span(ATV::kAtom_x)[0] = 0.25 * box_len;
    frame.span(ATV::kAtom_y)[0] = 0.25 * box_len;
    frame.span(ATV::kAtom_z)[0] = 0.25 * box_len;
    frame.span(ATV::kAtom_x)[1] = 0.75 * box_len;
    frame.span(ATV::kAtom_y)[1] = 0.75 * box_len;
    frame.span(ATV::kAtom_z)[1] = 0.75 * box_len;
    for (int i = 0; i < 2; ++i) {
        frame.span(ATV::kAtom_mass)[i] = 12.011;
        frame.span(ATV::kAtom_type)[i] = 6;
    }
    frame.simbox = ATV::cpu::Box::cubic(box_len);
    return frame;
}

/// Run a Deform component through `build()` + `N` `launch()` calls,
/// mirroring how `cpu::md::MD` would drive it at `initial_integrate`.
void run_deform(cpu::md::Deform& deform, cpu::Frame& frame, int n_steps) {
    cpu::Context runtime;
    cpu::BuildContext bctx{runtime, frame};
    deform.build(bctx);
    for (int step = 0; step < n_steps; ++step) {
        cpu::LaunchContext lctx{runtime, frame, /*stage=*/0, step};
        deform.launch(lctx);
    }
}

} // namespace

// ── ac-002 / ac-004: Mode::erate reproduces L(t) = L0*(1 + R*t) ─────────

BENCH_CASE(DeformLammpsAlign_ErateLinearEngineeringStrainTrajectory) {
    const double L0 = 40.0;   // Angstrom
    const double dt = 1.0;    // fs
    const double rate = 1e-3; // 1/fs
    const int N = 50;         // steps

    auto frame = make_two_atom_box(L0);
    cpu::md::Deform deform{cpu::md::Deform::Mode::erate, rate, dt, /*axis=*/0};
    run_deform(deform, frame, N);

    // LAMMPS fix deform ERATE: L = L0*(1 + R*t_elapsed), t_elapsed = N*dt.
    const double expected = L0 * (1.0 + rate * N * dt);
    BM_NEAR(frame.simbox.cell()(0, 0), expected, 1e-9)
        << "erate box edge diverged from L0*(1+R*t)";

    // y and z untouched.
    BM_DOUBLE_EQ(frame.simbox.cell()(1, 1), L0);
    BM_DOUBLE_EQ(frame.simbox.cell()(2, 2), L0);
}

// ── ac-002 / ac-003: Mode::trate reproduces L(t) = L0*exp(R*t) ──────────

BENCH_CASE(DeformLammpsAlign_TrateExponentialTrueStrainTrajectory) {
    const double L0 = 40.0;
    const double dt = 1.0;
    const double rate = 1e-3;
    const int N = 50;

    auto frame = make_two_atom_box(L0);
    cpu::md::Deform deform{cpu::md::Deform::Mode::trate, rate, dt, /*axis=*/0};
    run_deform(deform, frame, N);

    // LAMMPS fix deform TRATE: L = L0*exp(R*t_elapsed).
    const double expected = L0 * std::exp(rate * N * dt);
    BM_NEAR(frame.simbox.cell()(0, 0), expected, 1e-9)
        << "trate box edge diverged from L0*exp(R*t)";

    // erate and trate must differ beyond first order.
    const double erate_value = L0 * (1.0 + rate * N * dt);
    BM_GT(std::abs(frame.simbox.cell()(0, 0) - erate_value), 1e-6)
        << "trate must be distinct from erate beyond first order";
}

// ── erate affine remap: atom coords scale with the box edge ─────────────

BENCH_CASE(DeformLammpsAlign_ErateAffinelyRemapsAtomCoordinates) {
    const double L0 = 40.0;
    const double dt = 1.0;
    const double rate = 2e-3;
    const int N = 20;

    auto frame = make_two_atom_box(L0);
    const double x0_atom1 = frame.span(ATV::kAtom_x)[1];

    cpu::md::Deform deform{cpu::md::Deform::Mode::erate, rate, dt, /*axis=*/0};
    run_deform(deform, frame, N);

    const double box_factor = (1.0 + rate * N * dt);
    BM_NEAR(frame.span(ATV::kAtom_x)[1], x0_atom1 * box_factor, 1e-9)
        << "erate did not affinely remap atom x coordinate";
}

// ── ac-003: reference-frame regression -- the MSU failure mode ──────────
//
// The MSU stress-strain protocol compressed the box from 80 -> 73 Angstrom
// under NPT BEFORE the stage-5 Deform fix was issued. LAMMPS fix deform
// snapshots L0 at fix-issue time (lo_start/hi_start in fix_deform.cpp::init).
// The old Deform implementation per-step compounded from whatever edge
// existed, with no explicit snapshot, so referencing strain to the 80 A
// simulation-initial box gave all-negative strain. After this fix the
// strain must reference the 73 A activation-time length: a positive ramp.

BENCH_CASE(DeformLammpsAlign_ErateReferencesActivationLengthNotSimInitial) {
    const double L_sim_initial = 80.0; // box at frame construction
    const double L_activation = 73.0;  // box edge AFTER NPT compression
    const double dt = 1.0;
    const double rate = 1e-5; // matches lammps.in erate 1e-5
    const int N = 1000;

    // Frame built at 80 A, then "pre-compressed" to 73 A before Deform runs.
    auto frame = make_two_atom_box(L_sim_initial);
    frame.simbox = ATV::cpu::Box::cubic(L_activation);

    cpu::md::Deform deform{cpu::md::Deform::Mode::erate, rate, dt, /*axis=*/0};
    run_deform(deform, frame, N);

    // Engineering strain referenced to the ACTIVATION length (73 A).
    const double L_final = frame.simbox.cell()(0, 0);
    const double strain_vs_activation = (L_final - L_activation) / L_activation;
    const double strain_vs_sim_initial =
        (L_final - L_sim_initial) / L_sim_initial;

    // The LAMMPS-correct strain is positive and equals R*N*dt.
    BM_GT(strain_vs_activation, 0.0)
        << "strain must be a positive ramp referenced to activation length";
    BM_NEAR(strain_vs_activation, rate * N * dt, 1e-9)
        << "engineering strain must equal R*N*dt referenced to L_activation";

    // The buggy sim-initial reference would have reported NEGATIVE strain.
    BM_LT(strain_vs_sim_initial, 0.0)
        << "sanity: sim-initial reference is negative -- the old bug";

    // Box edge itself: L_activation*(1 + R*N*dt).
    BM_NEAR(L_final, L_activation * (1.0 + rate * N * dt), 1e-9);
}

// ── trate reference-frame regression: L0 snapshot at activation ─────────

BENCH_CASE(DeformLammpsAlign_TrateReferencesActivationLength) {
    const double L_activation = 73.0;
    const double dt = 1.0;
    const double rate = 1e-4;
    const int N = 100;

    auto frame = make_two_atom_box(80.0);
    frame.simbox = ATV::cpu::Box::cubic(L_activation);

    cpu::md::Deform deform{cpu::md::Deform::Mode::trate, rate, dt, /*axis=*/0};
    run_deform(deform, frame, N);

    BM_NEAR(frame.simbox.cell()(0, 0), L_activation * std::exp(rate * N * dt),
            1e-9)
        << "trate must reference the activation-time box edge";
}

// ── build() re-arms the snapshot -- a rebuilt driver re-references L0 ────

BENCH_CASE(DeformLammpsAlign_RebuildReSnapshotsActivationLength) {
    const double dt = 1.0;
    const double rate = 1e-3;
    const int N = 30;

    cpu::md::Deform deform{cpu::md::Deform::Mode::erate, rate, dt, /*axis=*/0};

    // First segment from L = 40.
    auto frame_a = make_two_atom_box(40.0);
    run_deform(deform, frame_a, N);
    BM_NEAR(frame_a.simbox.cell()(0, 0), 40.0 * (1.0 + rate * N * dt), 1e-9);

    // Second segment: same component, different box edge. build() must
    // re-snapshot L0 so strain references the new activation length.
    auto frame_b = make_two_atom_box(55.0);
    run_deform(deform, frame_b, N);
    BM_NEAR(frame_b.simbox.cell()(0, 0), 55.0 * (1.0 + rate * N * dt), 1e-9)
        << "build() must re-arm the L0 snapshot for a new segment";
}

// ── set_rate(0) freezes the box; component stays usable ─────────────────

BENCH_CASE(DeformLammpsAlign_SetRateZeroFreezesBox) {
    const double L0 = 40.0;
    const double dt = 1.0;

    auto frame = make_two_atom_box(L0);

    cpu::md::Deform deform{cpu::md::Deform::Mode::erate, 1e-3, dt, /*axis=*/0};
    deform.set_rate(0.0);
    BM_DOUBLE_EQ(deform.rate(), 0.0);
    BM_EQ(deform.mode(), cpu::md::Deform::Mode::erate);

    run_deform(deform, frame, 100);

    BM_DOUBLE_EQ(frame.simbox.cell()(0, 0), L0)
        << "set_rate(0) must freeze the box";
    BM_DOUBLE_EQ(frame.span(ATV::kAtom_x)[1], 0.75 * L0)
        << "set_rate(0) must leave atom coordinates untouched";
}

// ── negative rate compresses (sign convention matches fix deform) ───────

BENCH_CASE(DeformLammpsAlign_NegativeRateCompresses) {
    const double L0 = 40.0;
    const double dt = 1.0;
    const double rate = -5e-4;
    const int N = 50;

    auto frame = make_two_atom_box(L0);
    cpu::md::Deform deform{cpu::md::Deform::Mode::erate, rate, dt, /*axis=*/0};
    run_deform(deform, frame, N);

    BM_LT(frame.simbox.cell()(0, 0), L0)
        << "negative erate must compress the box";
    BM_NEAR(frame.simbox.cell()(0, 0), L0 * (1.0 + rate * N * dt), 1e-9);
}

// ── axis selection: y axis (axis=1) deforms, x and z untouched ──────────

BENCH_CASE(DeformLammpsAlign_AxisSelectionDeformsRequestedAxisOnly) {
    const double L0 = 40.0;
    const double dt = 1.0;
    const double rate = 1e-3;
    const int N = 30;

    auto frame = make_two_atom_box(L0);
    cpu::md::Deform deform{cpu::md::Deform::Mode::erate, rate, dt, /*axis=*/1};
    run_deform(deform, frame, N);

    BM_NEAR(frame.simbox.cell()(1, 1), L0 * (1.0 + rate * N * dt), 1e-9);
    BM_DOUBLE_EQ(frame.simbox.cell()(0, 0), L0);
    BM_DOUBLE_EQ(frame.simbox.cell()(2, 2), L0);
}

BENCH_MAIN("MD fix deform — LAMMPS API alignment");
