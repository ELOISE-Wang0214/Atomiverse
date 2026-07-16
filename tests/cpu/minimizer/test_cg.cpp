/// @file tests/cpu/minimizer/test_cg.cpp
/// @brief Tests for CPU Polak-Ribière+ CG with the strong-Wolfe line search.

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/minimizer/cg.h"
#include "atomiverse/cpu/minimizer/minimize.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_types.h"

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
} // namespace

// ac-008: CG with the default strong-Wolfe line search (c2=0.1) and PR+
// relaxes H2/STO-3G to the RHF minimum bond length.
TEST(CGTest, H2RHFStrongWolfeReachesMinimum) {
    auto frame = make_h2(1.2);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::minimizer::CG cg{/*alpha0=*/0.05};
    cpu::minimizer::Minimize opt;
    opt.set_optimizer_ref(cg);
    opt.add_handler(cpu::ForceConverged{/*ftol=*/0.01, /*rms_ftol=*/0.005,
                                        /*max_dtol=*/0.01, /*rms_dtol=*/0.005,
                                        /*etol=*/1e-7});
    opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCF{config});
    opt.build();
    opt.run(frame, 100);

    const double bond =
        frame.span(ATV::kAtom_x)[1] - frame.span(ATV::kAtom_x)[0];
    // H2/STO-3G RHF minimum bond ≈ 0.7122 Å.
    EXPECT_NEAR(static_cast<double>(bond), 0.7122, 0.02);
}

// Compressed H2 expands under CG.
TEST(CGTest, CompressedH2Expands) {
    auto frame = make_h2(0.5);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::minimizer::Minimize opt;
    opt.set_optimizer(cpu::minimizer::CG{});
    opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCF{config});
    opt.build();
    opt.run(frame, 60);

    const double bond =
        frame.span(ATV::kAtom_x)[1] - frame.span(ATV::kAtom_x)[0];
    EXPECT_GT(static_cast<double>(bond), 0.5) << "compressed H2 should expand";
}
