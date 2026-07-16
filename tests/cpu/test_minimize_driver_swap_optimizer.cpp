/// @file tests/cpu/test_minimize_driver_swap_optimizer.cpp
/// @brief Phase 2.6 — same `Minimize` driver accepts either `LBFGS` or `CG`
///        as the Optimizer component, and both converge to the same H₂ /
///        STO-3G minimum (within a loose tolerance).
///
/// Acceptance: ac-007 in `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/minimizer/cg.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

namespace {

cpu::Frame make_h2(double bond_length) {
    cpu::Frame f(2);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = bond_length;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_z)[1] = 0.0;
    f.span(ATV::kAtom_mass)[0] = 1.008;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.span(ATV::kAtom_type)[0] = 1;
    f.span(ATV::kAtom_type)[1] = 1;
    f.simbox = ATV::cpu::Box::cubic(100.0);
    return f;
}

double bond_length(const cpu::Frame& f) {
    const double dx = f.span(ATV::kAtom_x)[1] - f.span(ATV::kAtom_x)[0];
    const double dy = f.span(ATV::kAtom_y)[1] - f.span(ATV::kAtom_y)[0];
    const double dz = f.span(ATV::kAtom_z)[1] - f.span(ATV::kAtom_z)[0];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double energy_at(const cpu::Frame& frame_template) {
    auto integrals =
        cpu::scf::build_integrals(frame_template, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);
    ATV::cpu::scf::SCFSolver result_solver{config};
    result_solver.build(integrals);
    const auto& result = result_solver.run(integrals);
    return result.energy;
}

} // namespace

TEST(MinimizeOptimizerSwap, LBFGSAndCGReachSameH2Minimum) {
    // Optimize twice from the same starting geometry, identical configs,
    // identical max_iter — once with LBFGS, once with CG.
    constexpr double start_bond = 1.2; // Å (compressed/stretched)
    constexpr int max_iter = 200;

    // ── Run 1: LBFGS ───────────────────────────────────────────────
    auto frame_lbfgs = make_h2(start_bond);
    {
        auto integrals =
            cpu::scf::build_integrals(frame_lbfgs, cpu::scf::BasisSet::STO3G);
        auto cfg = cpu::scf::rhf_config(integrals);

        cpu::minimizer::Minimize opt;
        opt.add(cpu::minimizer::LBFGS{});
        opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCF{cfg});
        opt.build();
        opt.run(frame_lbfgs, max_iter);
    }
    const double bond_lbfgs = bond_length(frame_lbfgs);
    const double energy_lbfgs = energy_at(frame_lbfgs);

    // ── Run 2: CG ──────────────────────────────────────────────────
    auto frame_cg = make_h2(start_bond);
    {
        auto integrals =
            cpu::scf::build_integrals(frame_cg, cpu::scf::BasisSet::STO3G);
        auto cfg = cpu::scf::rhf_config(integrals);

        cpu::minimizer::Minimize opt;
        opt.add(cpu::minimizer::CG{});
        opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCF{cfg});
        opt.build();
        opt.run(frame_cg, max_iter);
    }
    const double bond_cg = bond_length(frame_cg);
    const double energy_cg = energy_at(frame_cg);

    // ── Sanity: both contracted from the start. ────────────────────
    EXPECT_LT(bond_lbfgs, start_bond) << "LBFGS should contract from 1.2 Å";
    EXPECT_LT(bond_cg, start_bond) << "CG should contract from 1.2 Å";
    EXPECT_GT(bond_lbfgs, 0.5);
    EXPECT_GT(bond_cg, 0.5);

    // ── Headline: same minimum within a structurally loose tolerance.
    // The 0.05 Å bond gap and 1e-3 Ha (~ 0.027 eV) energy gap are sloppy
    // enough that a basic CG without line search still passes; the test
    // is *structural* (Minimize accepts both as Optimizer components),
    // not a numerical-quality test.
    EXPECT_NEAR(bond_lbfgs, bond_cg, 0.05)
        << "LBFGS and CG should land at the same H₂ minimum bond length";
    EXPECT_NEAR(energy_lbfgs, energy_cg, 1e-3 * units::kHartreeInEv)
        << "LBFGS and CG should land at the same H₂ minimum energy "
        << "(within 1e-3 Ha)";
}
