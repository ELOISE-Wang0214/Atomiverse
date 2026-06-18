/// @file tests/cpu/gpw/test_cpu_gpw_evaluator.cpp
/// @brief GPWEvaluator plugs into the Evaluate driver and writes frame energy.

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/driver/single_point.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"
#include "atomiverse/facade/component_role.h"

using namespace ATV::cpu::gpw;

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;
constexpr double kBohr2Ang = 1.0 / 1.8897259885789;
} // namespace

TEST(GpwEvaluator, H2SinglePointThroughDriver) {
    // H2 in an 8x4x4 Bohr box at 1 Bohr separation (matches the H2/CP2K test),
    // expressed in Angstrom for the Frame.
    ATV::cpu::Frame frame(2);
    frame.box = ATV::cpu::Box::orthorhombic(8.0 * kBohr2Ang, 4.0 * kBohr2Ang, 4.0 * kBohr2Ang);
    frame.atoms.x(0) = 3.5 * kBohr2Ang;
    frame.atoms.y(0) = 2.0 * kBohr2Ang;
    frame.atoms.z(0) = 2.0 * kBohr2Ang;
    frame.atoms.x(1) = 4.5 * kBohr2Ang;
    frame.atoms.y(1) = 2.0 * kBohr2Ang;
    frame.atoms.z(1) = 2.0 * kBohr2Ang;
    frame.atoms.type(0) = 1; // H
    frame.atoms.type(1) = 1;

    GpwElement H;
    H.z = 1;
    H.basis = parse_gth_basis(kData + "/BASIS_SET", "H", "DZV-GTH-PADE");
    H.pot = parse_gth_potential(kData + "/POTENTIAL", "H", "GTH-PADE-q1");

    GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = 140.0;
    cfg.xc = XcKind::PADE;
    cfg.nonperiodic = false;
    cfg.mix = 0.3;
    cfg.elements = {H};

    GPWEvaluator eval(cfg);
    ATV::cpu::single_point(frame, eval);

    // The evaluator computed n_occ=1 from 2 valence electrons, ran the SCF, and
    // wrote eV onto the frame. Compare to CP2K H2 (-1.06345575 Ha) in eV.
    const double e_ev = frame.scalars.energy;
    const double ref_ev = -1.06345574964568 * 27.211386245988;
    EXPECT_TRUE(eval.last_result().converged);
    EXPECT_NEAR(e_ev, ref_ev, 0.1) // ~3.7 mHa window (same as the direct H2 test)
        << "frame.energy=" << e_ev << " eV  ref=" << ref_ev << " eV";
    // Frame energy must equal the solver's Hartree result times the conversion.
    EXPECT_NEAR(e_ev, eval.last_energy_hartree() * 27.211386245988, 1e-9);
}
