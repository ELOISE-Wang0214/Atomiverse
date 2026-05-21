/// @file tests/cpu/vib/test_tr_projection.cpp
/// @brief Unit tests for translation/rotation projection, rotor-type
///        classification, genuine mode count, and imaginary-mode
///        handling (spec phase4-vibrations-02).
///
/// These tests drive `cpu::vib::MassWeighting` + `cpu::vib::EigenSolver`
/// directly on a hand-built Hessian via the shared `HessianWorkspace`,
/// so they exercise the TR projector without an SCF run.
///
/// Coverage map:
///   ac-003 NonlinearWaterYields3Modes   — projector removes 6 TR vectors
///   ac-003 LinearMoleculeYields3Nminus5 — is_linear, 3N-5 modes
///   ac-004 ImaginaryModeCountedAndNegative

#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/vib/eigen_solver.h"
#include "atomiverse/cpu/vib/hessian_workspace.h"
#include "atomiverse/cpu/vib/mass_weighting.h"
#include "atomiverse/core/mat3.h"

using namespace ATV;

namespace {

/// @brief Run MassWeighting then EigenSolver on a workspace whose `r0`,
///        `hessian`, and atom types are pre-populated.
///
/// @param frame Frame carrying atom types (atomic numbers).
/// @param ws    Workspace with `r0` and `hessian` filled.
void run_diagonalize(cpu::Frame& frame, cpu::vib::HessianWorkspace& ws) {
    cpu::Context runtime;
    const int dim = 3 * ws.n_atoms;
    cpu::LaunchContext ctx{runtime, frame, /*stage=*/0, /*step=*/2 * dim};
    cpu::vib::MassWeighting{&ws}.launch(ctx);
    cpu::vib::EigenSolver{&ws}.launch(ctx);
}

/// @brief Optimized-ish H2O geometry (nonlinear), Angstrom.
cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.x = {0.0, 0.7572, -0.7572};
    f.atoms.y = {0.0, 0.5865,  0.5865};
    f.atoms.z = {0.0, 0.0,     0.0};
    f.atoms.type = {8, 1, 1};
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief Linear CO2 geometry: O=C=O along the x-axis, Angstrom.
cpu::Frame make_co2() {
    cpu::Frame f(3);
    f.atoms.x = {0.0, 1.16, -1.16}; // C at origin, two O symmetric.
    f.atoms.y = {0.0, 0.0,  0.0};
    f.atoms.z = {0.0, 0.0,  0.0};
    f.atoms.type = {6, 8, 8};
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief Diatomic HF along the x-axis, Angstrom.
cpu::Frame make_hf() {
    cpu::Frame f(2);
    f.atoms.x = {0.0, 0.917};
    f.atoms.y = {0.0, 0.0};
    f.atoms.z = {0.0, 0.0};
    f.atoms.type = {1, 9}; // H, F
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief Fill the workspace `r0` from the frame coordinates (Angstrom).
void fill_r0(const cpu::Frame& f, cpu::vib::HessianWorkspace& ws) {
    const int n = f.atoms.size();
    for (int i = 0; i < n; ++i) {
        ws.r0[3 * i + 0] = f.atoms.x(i);
        ws.r0[3 * i + 1] = f.atoms.y(i);
        ws.r0[3 * i + 2] = f.atoms.z(i);
    }
    ws.r0_ready = true;
}

} // namespace

// ── ac-003: nonlinear H2O yields exactly 3N-6 = 3 finite modes ──────────

TEST(VibTRProjection, NonlinearWaterYields3Modes) {
    auto frame = make_h2o();
    const int n = 3, dim = 9;

    cpu::vib::HessianWorkspace ws;
    ws.allocate(n);
    fill_r0(frame, ws);

    // A physically reasonable symmetric positive-definite Cartesian
    // Hessian: a diagonally-dominant SPD matrix. Its content does not
    // affect the projected mode COUNT — that is fixed by the TR
    // subspace rank — but a finite SPD matrix keeps the modes finite.
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j)
            ws.hessian(i, j) = 0.1 * std::cos(0.7 * i + 1.3 * j);
        ws.hessian(i, i) += 10.0; // diagonal dominance -> SPD.
    }
    // Symmetrize (the driver's Symmetrize would have done this).
    for (int i = 0; i < dim; ++i)
        for (int j = i + 1; j < dim; ++j) {
            const double a = 0.5 * (ws.hessian(i, j) + ws.hessian(j, i));
            ws.hessian(i, j) = ws.hessian(j, i) = a;
        }

    run_diagonalize(frame, ws);

    EXPECT_FALSE(ws.is_linear) << "H2O must be classified nonlinear";
    EXPECT_EQ(ws.frequencies.size(), 3u)
        << "TR projection must leave exactly 3N-6 = 3 modes";
    EXPECT_EQ(ws.normal_modes.shape()[0], 9u);
    EXPECT_EQ(ws.normal_modes.shape()[1], 3u);
}

// ── ac-003: linear CO2 yields exactly 3N-5 = 4 finite modes ─────────────

TEST(VibTRProjection, LinearMoleculeYields3Nminus5) {
    auto frame = make_co2();
    const int n = 3, dim = 9;

    cpu::vib::HessianWorkspace ws;
    ws.allocate(n);
    fill_r0(frame, ws);

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j)
            ws.hessian(i, j) = 0.1 * std::cos(0.5 * i + 1.1 * j);
        ws.hessian(i, i) += 10.0;
    }
    for (int i = 0; i < dim; ++i)
        for (int j = i + 1; j < dim; ++j) {
            const double a = 0.5 * (ws.hessian(i, j) + ws.hessian(j, i));
            ws.hessian(i, j) = ws.hessian(j, i) = a;
        }

    run_diagonalize(frame, ws);

    EXPECT_TRUE(ws.is_linear) << "CO2 must be classified linear";
    EXPECT_EQ(ws.frequencies.size(), 4u)
        << "linear molecule must leave exactly 3N-5 = 4 modes";
}

// ── ac-004: a negative-curvature Hessian yields one imaginary mode ──────

TEST(VibTRProjection, ImaginaryModeCountedAndNegative) {
    // Diatomic HF on the x-axis. The standard spring Hessian for the
    // stretch coordinate: H(0,0)=H(3,3)=k, H(0,3)=H(3,0)=-k. With k<0
    // the single stretch mode is imaginary; the 5 TR modes are removed.
    auto frame = make_hf();
    const int n = 2, dim = 6;

    cpu::vib::HessianWorkspace ws;
    ws.allocate(n);
    fill_r0(frame, ws);

    const double k = -3.0; // negative curvature -> imaginary stretch.
    ws.hessian(0, 0) = k;
    ws.hessian(3, 3) = k;
    ws.hessian(0, 3) = -k;
    ws.hessian(3, 0) = -k;

    run_diagonalize(frame, ws);

    // Diatomic: 3N-5 = 1 genuine vibrational mode.
    ASSERT_EQ(ws.frequencies.size(), 1u);
    EXPECT_EQ(ws.n_imaginary, 1)
        << "exactly one mode must be flagged imaginary";
    EXPECT_LT(ws.frequencies(0), 0.0)
        << "the imaginary mode must be reported as a negative wavenumber";
}
