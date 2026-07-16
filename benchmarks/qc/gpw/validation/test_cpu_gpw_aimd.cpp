/// @file tests/cpu/gpw/test_cpu_gpw_aimd.cpp
/// @brief GPW ab-initio MD: a short NVE trajectory of H2 driven by GPWEvaluator
///        forces through the MD driver — the tutorial's MD "Next Step".

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/md/fix_csvr.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"

using namespace ATV::cpu::gpw;

namespace {

double bond_ang(const ATV::cpu::Frame& f) {
    const double dx = f.span(ATV::kAtom_x)[1] - f.span(ATV::kAtom_x)[0];
    const double dy = f.span(ATV::kAtom_y)[1] - f.span(ATV::kAtom_y)[0];
    const double dz = f.span(ATV::kAtom_z)[1] - f.span(ATV::kAtom_z)[0];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
} // namespace

TEST(GpwAIMD, H2NveTrajectoryRuns) {
    ATV::cpu::Frame frame(2);
    frame.simbox = ATV::cpu::Box::orthorhombic(
        8.0 * ATV::units::kBohrInAngstrom, 4.0 * ATV::units::kBohrInAngstrom,
        4.0 * ATV::units::kBohrInAngstrom);
    frame.span(ATV::kAtom_x)[0] =
        3.6 * ATV::units::kBohrInAngstrom; // ~1.4 Bohr near equilibrium
    frame.span(ATV::kAtom_x)[1] = 5.0 * ATV::units::kBohrInAngstrom;
    frame.span(ATV::kAtom_y)[0] = frame.span(ATV::kAtom_y)[1] =
        2.0 * ATV::units::kBohrInAngstrom;
    frame.span(ATV::kAtom_z)[0] = frame.span(ATV::kAtom_z)[1] =
        2.0 * ATV::units::kBohrInAngstrom;
    frame.span(ATV::kAtom_type)[0] = frame.span(ATV::kAtom_type)[1] = 1;
    frame.span(ATV::kAtom_mass)[0] = frame.span(ATV::kAtom_mass)[1] =
        1.008; // H [amu]

    const Kind H{molrs::Element::H,
                 GaussianBasis(molrs::Element::H, "DZV-GTH-PADE"),
                 GthPotential(molrs::Element::H, "GTH-PADE-q1")};

    GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = 140.0;
    cfg.mix = 0.3;
    cfg.compute_forces = true;
    cfg.fd_step = 2e-3;
    cfg.kinds = {H};

    const double bond0 = bond_ang(frame);

    ATV::cpu::md::MD sim;
    sim.add(ATV::cpu::md::InitVelocity{300.0, /*seed=*/7u}); // 300 K
    sim.add(ATV::cpu::md::VelocityVerlet{0.0005});           // 0.5 fs
    sim.on(ATV::md::Phase::force).add(GPWEvaluator{cfg});
    sim.build();
    sim.run(frame, 4);

    // The trajectory ran: energy is finite and the bond has moved (atoms have
    // nonzero velocity at 300 K, integrated by velocity-Verlet on GPW forces).
    EXPECT_TRUE(
        std::isfinite(static_cast<double>(frame.meta.get(ATV::kEnergy))));
    EXPECT_GT(std::abs(bond_ang(frame) - bond0), 1e-6)
        << "frame did not evolve";
}

// GPW-CSVR NVT: the same textbook H2, thermostatted by FixCSVR (Bussi). The
// canonical-distribution and H-tilde drift statistics are proven on classical
// Ar in benchmarks/md/lj/validation/test_md_csvr.cpp; this is the
// integration gate — GPW forces + CSVR run together, the trajectory stays
// physical and the Bussi conserved quantity is finite.
TEST(GpwAIMD, H2CsvrNvtTrajectoryRuns) {
    ATV::cpu::Frame frame(2);
    frame.simbox = ATV::cpu::Box::orthorhombic(
        8.0 * ATV::units::kBohrInAngstrom, 4.0 * ATV::units::kBohrInAngstrom,
        4.0 * ATV::units::kBohrInAngstrom);
    frame.span(ATV::kAtom_x)[0] = 3.6 * ATV::units::kBohrInAngstrom;
    frame.span(ATV::kAtom_x)[1] = 5.0 * ATV::units::kBohrInAngstrom;
    frame.span(ATV::kAtom_y)[0] = frame.span(ATV::kAtom_y)[1] =
        2.0 * ATV::units::kBohrInAngstrom;
    frame.span(ATV::kAtom_z)[0] = frame.span(ATV::kAtom_z)[1] =
        2.0 * ATV::units::kBohrInAngstrom;
    frame.span(ATV::kAtom_type)[0] = frame.span(ATV::kAtom_type)[1] = 1;
    frame.span(ATV::kAtom_mass)[0] = frame.span(ATV::kAtom_mass)[1] = 1.008;

    const Kind H{molrs::Element::H,
                 GaussianBasis(molrs::Element::H, "DZV-GTH-PADE"),
                 GthPotential(molrs::Element::H, "GTH-PADE-q1")};

    GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = 140.0;
    cfg.mix = 0.3;
    cfg.compute_forces = true;
    cfg.fd_step = 2e-3;
    cfg.kinds = {H};

    const double dt = 0.0005; // 0.5 fs
    ATV::cpu::md::MD sim;
    sim.add(ATV::cpu::md::InitVelocity{330.0, /*seed=*/7u});
    sim.add(ATV::cpu::md::VelocityVerlet{dt});
    sim.set_thermostat(ATV::cpu::md::FixCSVR{330.0, /*tau [fs]=*/20.0, dt, 7u});
    sim.on(ATV::md::Phase::force).add(GPWEvaluator{cfg});
    sim.build();
    sim.run(frame, 10);

    EXPECT_TRUE(
        std::isfinite(static_cast<double>(frame.meta.get(ATV::kEnergy))));
    EXPECT_TRUE(
        std::isfinite(static_cast<double>(frame.meta.get(ATV::kConserved))))
        << "Bussi H-tilde bookkeeping broke";
    const double bond = bond_ang(frame);
    EXPECT_GT(bond, 0.4);
    EXPECT_LT(bond, 1.5) << "H2 dissociated in 10 NVT steps — unphysical";
}
