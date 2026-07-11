/// @file tests/cpu/gpw/test_cpu_gpw_geoopt.cpp
/// @brief GPW geometry optimization: H2 relaxes toward its equilibrium bond
///        through the Minimize driver + LBFGS, driven by GPWEvaluator forces.
///
/// This is the getting-started tutorial's "Next Step" (geometry optimization),
/// realized with GPW: the same Driver/Component vocabulary, with a GPWEvaluator
/// advertising Capability::Force (finite-difference). Starting from a compressed
/// 1.0 Bohr bond, the optimizer should lengthen it toward the H2 equilibrium
/// (~1.4 Bohr).

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/minimizer/minimize.h"

using namespace ATV::cpu::gpw;

namespace {
constexpr double kBohr2Ang = 1.0 / 1.8897259885789;

double bond_bohr(const ATV::cpu::Frame& f) {
    const double dx = (f.atoms.x(1) - f.atoms.x(0)) / kBohr2Ang;
    const double dy = (f.atoms.y(1) - f.atoms.y(0)) / kBohr2Ang;
    const double dz = (f.atoms.z(1) - f.atoms.z(0)) / kBohr2Ang;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
} // namespace

TEST(GpwGeoOpt, H2RelaxesTowardEquilibrium) {
    ATV::cpu::Frame frame(2);
    // 14x12x12 Bohr, not the old 8x4x4: at 4 Bohr the molecule interacted with its
    // own periodic images hard enough to push the equilibrium bond out to 1.73 Bohr,
    // so the test was asserting a physical bond length in a box that had no physical
    // equilibrium to find. Enlarged until the isolated-molecule minimum is recovered.
    frame.box = ATV::cpu::Box::orthorhombic(14.0 * kBohr2Ang, 12.0 * kBohr2Ang, 12.0 * kBohr2Ang);
    frame.atoms.x(0) = 6.5 * kBohr2Ang; // 1.0 Bohr apart along x (compressed)
    frame.atoms.x(1) = 7.5 * kBohr2Ang;
    frame.atoms.y(0) = frame.atoms.y(1) = 6.0 * kBohr2Ang;
    frame.atoms.z(0) = frame.atoms.z(1) = 6.0 * kBohr2Ang;
    frame.atoms.type(0) = 1;
    frame.atoms.type(1) = 1;

    Kind H;
    H.element = ATV::Element::H;
    H.basis = GaussianBasis(ATV::Element::H, "DZV-GTH-PADE");
    H.pot = GthPotential(ATV::Element::H, "GTH-PADE-q1");

    GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = 140.0;
    cfg.mix = 0.3;
    cfg.compute_forces = true;
    cfg.fd_step = 2e-3;
    cfg.kinds = {H};

    const double bond0 = bond_bohr(frame);
    ASSERT_NEAR(bond0, 1.0, 1e-6);

    ATV::cpu::minimizer::Minimize opt;
    opt.add(ATV::cpu::minimizer::LBFGS{5, 0.01, 0.05, 10, 0.03, 0.1, 0.05, 1e-3});
    opt.on(ATV::minimizer::Phase::force).add(GPWEvaluator{cfg});
    opt.build();
    // 30 steps, not 8: the assertion below is about WHERE the optimizer settles, so
    // it has to be allowed to settle. (The old test stopped at 8 steps and capped
    // the bond at 1.6 Bohr. That cap was calibrated against forces that were wrong
    // — the analytic gradient was taken of a different functional than the SCF
    // converged, see internal/local_pp.h — so it encoded a mid-flight position of a
    // trajectory driven by a bad gradient, not a physical bound.)
    opt.run(frame, 30);

    const double bond1 = bond_bohr(frame);
    // The compressed bond must lengthen and settle at the H2 equilibrium, which for
    // GTH-PADE/DZV sits near 1.4-1.5 Bohr (gas-phase H2 is 1.40).
    EXPECT_GT(bond1, bond0 + 0.02) << "bond went " << bond0 << " -> " << bond1 << " Bohr";
    EXPECT_GT(bond1, 1.25) << "undershoot: bond=" << bond1 << " Bohr";
    EXPECT_LT(bond1, 1.70) << "overshoot: bond=" << bond1 << " Bohr";
}
