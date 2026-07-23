/// @file tests/cpu/minimize/test_cell_relax.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/driver/criteria.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/minimize/cell_relax.h"
#include "atomiverse/cpu/minimize/minimize.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

namespace {

// Two LJ atoms in a cubic box; short L → high pressure.
cpu::Frame compressed_dimer(double L) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_x)[1] = L * 0.5;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_mass)[i] = 39.948;
        f.span(kAtom_type)[i] = 0;
        f.span(kAtom_tag)[i] = i + 1;
    }
    f.simbox = cpu::Box::cubic(L);
    return f;
}

} // namespace

TEST(IsotropicCellRelax, NameAndKind) {
    cpu::minimizer::IsotropicCellRelax cell;
    EXPECT_EQ(cell.name(), "isotropic_cell_relax");
    EXPECT_NE(dynamic_cast<cpu::Optimizer*>(&cell), nullptr);
}

TEST(IsotropicCellRelax, ExpandsCompressedBox) {
    const double eps = 0.0104, sigma = 3.4, rc = 2.5 * sigma;
    // Box too small → high virial pressure → should expand toward P=0.
    auto frame = compressed_dimer(6.0);
    const double L0 = frame.simbox.length(0);

    cpu::minimizer::Minimize opt;
    cpu::minimizer::IsotropicCellRelax cell(0.0, 0.1);
    cell.max_scale = 0.1;
    opt.set_optimizer_ref(cell);
    // Replace default ForceConverged with pressure criterion.
    // Minimize always adds ForceConverged in install_defaults — add pressure
    // with Stop severity preference: user can clear? Check if we can only use
    // PressureConverged by adding it — ForceConverged may stop early on forces.
    opt.add_criterion(cpu::PressureConverged{0.0, 1e-3});
    auto lj = cpu::PairLJCut{eps, sigma, rc};
    lj.set_compute_stress(true);
    opt.add_potential(std::move(lj));
    opt.set_neighbor(cpu::md::NeighborList<>{rc, 0.3});
    opt.build();
    opt.run(frame, 80);

    const double L1 = frame.simbox.length(0);
    EXPECT_GT(L1, L0); // expanded
    EXPECT_TRUE(std::isfinite(L1));
}

TEST(PressureConverged, Name) {
    cpu::PressureConverged pc(0.0, 1e-4);
    EXPECT_EQ(pc.name(), "pressure_converged");
}
