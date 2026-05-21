/// @file tests/cpu/vib/test_mass_weighting.cpp
/// @brief Unit tests for the nist constants, IUPAC-2013 masses, and the
///        cpu::vib::MassWeighting transform (spec phase4-vibrations-02).
///
/// Coverage map:
///   ac-001 ConstantsIntegrity / MassTableValues
///   ac-002 MassWeightingTransform / MassLookupThrows

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/vib/atomic_masses.h"
#include "atomiverse/cpu/vib/hessian_workspace.h"
#include "atomiverse/cpu/vib/mass_weighting.h"
#include "atomiverse/cpu/vib/nist_constants.h"
#include "atomiverse/core/mat3.h"

using namespace ATV;

// ── ac-001: nist constants compose to the expected cm^-1 factor ─────────

TEST(VibConstants, ConstantsIntegrity) {
    // au2hz = sqrt(HARTREE2J/(ATOMIC_MASS*BOHR_SI^2))/(2*pi) ~ 1.541e14 Hz
    // per unit freq_au; the freq_au->cm^-1 factor au2hz/c*1e-2 ~ 5140.487.
    // (The spec's Domain-basis parenthetical "~2.194746e5" is a doc slip:
    // 2.194746e5 is nist.HARTREE2WAVENUMBER, a different quantity. The
    // spec FORMULA, mirrored verbatim here, yields PySCF's 5140.487138 —
    // confirmed against a live `pyscf.data.nist` computation.)
    EXPECT_NEAR(cpu::vib::nist::au2hz(), 1.5410792745663e14, 1e3)
        << "au2hz drifted: " << cpu::vib::nist::au2hz();
    const double au2cm = cpu::vib::nist::au2wavenumber();
    EXPECT_NEAR(au2cm, 5140.487138493, 1e-6)
        << "au2wavenumber drifted: " << au2cm;

    // HARTREE2J must be composed from HBAR/E_MASS/BOHR_SI, ~4.359744e-18 J.
    EXPECT_NEAR(cpu::vib::nist::HARTREE2J, 4.359744650e-18, 1e-24);

    // HARTREE2EV must land on the CODATA value ~27.2114 eV.
    EXPECT_NEAR(cpu::vib::nist::HARTREE2EV, 27.21138, 1e-3);
}

// ── ac-001: IUPAC-2013 atomic masses for the test elements ─────────────

TEST(VibConstants, MassTableValues) {
    EXPECT_DOUBLE_EQ(cpu::vib::mass_amu(1), 1.008);   // H
    EXPECT_DOUBLE_EQ(cpu::vib::mass_amu(6), 12.011);  // C
    EXPECT_DOUBLE_EQ(cpu::vib::mass_amu(7), 14.007);  // N
    EXPECT_DOUBLE_EQ(cpu::vib::mass_amu(8), 15.999);  // O
}

// ── ac-002: mass_amu throws on an unsupported atomic number ────────────

TEST(VibConstants, MassLookupThrows) {
    EXPECT_THROW(cpu::vib::mass_amu(0), std::out_of_range);
    EXPECT_THROW(cpu::vib::mass_amu(-3), std::out_of_range);
    EXPECT_THROW(cpu::vib::mass_amu(99), std::out_of_range);
}

// ── ac-002: MassWeighting divides the Hessian by sqrt(m_i m_j) ──────────

TEST(VibMassWeighting, MassWeightingTransform) {
    // 2-atom system: atom 0 = O (15.999), atom 1 = H (1.008).
    cpu::Frame frame(2);
    frame.atoms.x = {0.0, 0.96};
    frame.atoms.y = {0.0, 0.0};
    frame.atoms.z = {0.0, 0.0};
    frame.atoms.type = {8, 1};
    frame.box = Mat3::cubic(50.0);

    const int n   = 2;
    const int dim = 6;

    cpu::vib::HessianWorkspace ws;
    ws.allocate(n);

    // Hand-built 6x6 Hessian: H(i,j) = (i+1)*(j+1) (symmetric, nonzero).
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            ws.hessian(i, j) = static_cast<double>((i + 1) * (j + 1));

    // Keep a copy of the raw Hessian for the post-condition check.
    auto raw = ws.hessian;

    cpu::Context runtime;
    cpu::vib::MassWeighting mw{&ws};
    // Finalize pass: ctx.step == 6N == 12.
    cpu::LaunchContext ctx{runtime, frame,
                           /*stage=*/0, /*step=*/2 * dim};
    mw.launch(ctx);

    ASSERT_TRUE(ws.mass_weighted);

    // Expected: H^mw(ic,jc) = H(ic,jc) * unit_conv / sqrt(m[ic/3] m[jc/3]).
    const double m0 = cpu::vib::mass_amu(8); // O
    const double m1 = cpu::vib::mass_amu(1); // H
    const double mass[2] = {m0, m1};
    const double conv = cpu::vib::nist::EV_PER_ANG2_TO_HARTREE_PER_BOHR2;

    for (int ic = 0; ic < dim; ++ic) {
        for (int jc = 0; jc < dim; ++jc) {
            const double expected =
                raw(ic, jc) * conv /
                std::sqrt(mass[ic / 3] * mass[jc / 3]);
            EXPECT_NEAR(ws.hessian_mw(ic, jc), expected, 1e-14)
                << "mass-weighted entry (" << ic << "," << jc << ")";
        }
    }

    // The raw Hessian must be left untouched (link-1 units preserved).
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            EXPECT_DOUBLE_EQ(ws.hessian(i, j), raw(i, j));
}
