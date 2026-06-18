/// @file tests/cpu/gpw/test_cpu_gpw_aimd.cpp
/// @brief GPW ab-initio MD: a short NVE trajectory of H2 driven by GPWEvaluator
///        forces through the MD driver — the tutorial's MD "Next Step".

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"

using namespace ATV::cpu::gpw;

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;
constexpr double kBohr2Ang = 1.0 / 1.8897259885789;

double bond_ang(const ATV::cpu::Frame& f) {
    const double dx = f.atoms.x(1) - f.atoms.x(0);
    const double dy = f.atoms.y(1) - f.atoms.y(0);
    const double dz = f.atoms.z(1) - f.atoms.z(0);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
} // namespace

TEST(GpwAIMD, H2NveTrajectoryRuns) {
    ATV::cpu::Frame frame(2);
    frame.box = ATV::cpu::Box::orthorhombic(8.0 * kBohr2Ang, 4.0 * kBohr2Ang, 4.0 * kBohr2Ang);
    frame.atoms.x(0) = 3.6 * kBohr2Ang; // ~1.4 Bohr near equilibrium
    frame.atoms.x(1) = 5.0 * kBohr2Ang;
    frame.atoms.y(0) = frame.atoms.y(1) = 2.0 * kBohr2Ang;
    frame.atoms.z(0) = frame.atoms.z(1) = 2.0 * kBohr2Ang;
    frame.atoms.type(0) = frame.atoms.type(1) = 1;
    frame.atoms.mass(0) = frame.atoms.mass(1) = 1.008; // H [amu]

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

    const double bond0 = bond_ang(frame);

    ATV::cpu::md::MD sim;
    sim.add(ATV::cpu::md::InitVelocity{300.0, /*seed=*/7u}); // 300 K
    sim.add(ATV::cpu::md::VelocityVerlet{0.0005});           // 0.5 fs
    sim.on(ATV::md::Phase::force).add(GPWEvaluator{cfg});
    sim.build();
    sim.run(frame, 4);

    // The trajectory ran: energy is finite and the bond has moved (atoms have
    // nonzero velocity at 300 K, integrated by velocity-Verlet on GPW forces).
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.scalars.energy)));
    EXPECT_GT(std::abs(bond_ang(frame) - bond0), 1e-6) << "frame did not evolve";
}
