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
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/minimizer/minimize.h"

using namespace ATV::cpu::gpw;

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;
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
    frame.box = ATV::cpu::Box::orthorhombic(8.0 * kBohr2Ang, 4.0 * kBohr2Ang, 4.0 * kBohr2Ang);
    frame.atoms.x(0) = 3.5 * kBohr2Ang; // 1.0 Bohr apart along x (compressed)
    frame.atoms.x(1) = 4.5 * kBohr2Ang;
    frame.atoms.y(0) = frame.atoms.y(1) = 2.0 * kBohr2Ang;
    frame.atoms.z(0) = frame.atoms.z(1) = 2.0 * kBohr2Ang;
    frame.atoms.type(0) = 1;
    frame.atoms.type(1) = 1;

    GpwElement H;
    H.z = 1;
    H.basis = parse_gth_basis(kData + "/BASIS_SET", "H", "DZV-GTH-PADE");
    H.pot = parse_gth_potential(kData + "/POTENTIAL", "H", "GTH-PADE-q1");

    GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = 140.0;
    cfg.mix = 0.3;
    cfg.compute_forces = true;
    cfg.fd_step = 2e-3;
    cfg.elements = {H};

    const double bond0 = bond_bohr(frame);
    ASSERT_NEAR(bond0, 1.0, 1e-6);

    ATV::cpu::minimizer::Minimize opt;
    opt.add(ATV::cpu::minimizer::LBFGS{5, 0.01, 0.05, 10, 0.03, 0.1, 0.05, 1e-3});
    opt.on(ATV::minimizer::Phase::force).add(GPWEvaluator{cfg});
    opt.build();
    opt.run(frame, 8);

    const double bond1 = bond_bohr(frame);
    // The compressed bond should lengthen toward the ~1.4 Bohr equilibrium.
    EXPECT_GT(bond1, bond0 + 0.02) << "bond went " << bond0 << " -> " << bond1 << " Bohr";
    EXPECT_LT(bond1, 1.6) << "overshoot: bond=" << bond1 << " Bohr";
}
