/// @file tests/cpu/minimize/test_cell_relax.cpp
#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/potential.h"
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

namespace {

/// Constant intensive hydrostatic stress writer (PW/GPW kStress convention).
class ConstHydrostaticStress : public cpu::Potential {
  public:
    explicit ConstHydrostaticStress(double s) : s_(s) {}
    std::string name() const override { return "const_hydro_stress"; }
    void launch(const cpu::LaunchContext& ctx) override {
        if (!ctx.at(ATV::minimizer::Phase::force))
            return;
        ctx.frame.meta.get(kEnergy) = 0.0;
        auto& st = ctx.frame.meta.get(kStress);
        st[0] = st[1] = st[2] = s_;
        st[3] = st[4] = st[5] = 0.0;
    }

  private:
    double s_;
};

} // namespace

TEST(IsotropicCellRelax, IntensiveStressExpandsCompressed) {
    // Mechanical intensive σ > 0 (compressed) → s = 1 + alpha*P/B > 1.
    cpu::Frame frame(1);
    const double L0 = 10.0;
    frame.span(kAtom_x)[0] = frame.span(kAtom_y)[0] = frame.span(kAtom_z)[0] =
        0.25 * L0;
    frame.span(kAtom_mass)[0] = 1.0;
    frame.span(kAtom_type)[0] = 0;
    frame.span(kAtom_tag)[0] = 1;
    frame.simbox = cpu::Box::cubic(L0);

    cpu::minimizer::Minimize opt;
    cpu::minimizer::IsotropicCellRelax cell(0.0, /*alpha*/ 0.5);
    cell.bulk_mod = 1.0;
    cell.max_scale = 0.05;
    cell.stress_is_virial = false;
    opt.set_optimizer_ref(cell);
    opt.on(minimizer::Phase::force).add(ConstHydrostaticStress{0.02});
    opt.enable_timing(false);
    opt.build();
    opt.run(frame, 2); // seed + one scale step

    const double expected = L0 * (1.0 + 0.5 * 0.02 / 1.0);
    EXPECT_NEAR(frame.simbox.length(0), expected, 1e-12 * L0);
}

TEST(IsotropicCellRelax, FullMatrixScaleOnFccPrimitive) {
    // Diamond-like fcc cell has zero diagonal; isotropic scale must grow
    // off-diagonal entries (pre-product scale_isotropic was a no-op).
    cpu::Frame frame(1);
    frame.span(kAtom_x)[0] = frame.span(kAtom_y)[0] = frame.span(kAtom_z)[0] =
        0.0;
    frame.span(kAtom_mass)[0] = 1.0;
    frame.span(kAtom_type)[0] = 0;
    frame.span(kAtom_tag)[0] = 1;
    const double h = 2.715;
    Mat3 cell;
    cell(0, 0) = 0;
    cell(0, 1) = h;
    cell(0, 2) = h;
    cell(1, 0) = h;
    cell(1, 1) = 0;
    cell(1, 2) = h;
    cell(2, 0) = h;
    cell(2, 1) = h;
    cell(2, 2) = 0;
    frame.simbox.set_cell(cell);
    const double vol0 = frame.simbox.volume();

    cpu::minimizer::Minimize opt;
    cpu::minimizer::IsotropicCellRelax relax(0.0, /*alpha*/ 1.0);
    relax.bulk_mod = 1.0;
    relax.max_scale = 0.05;
    relax.stress_is_virial = false;
    opt.set_optimizer_ref(relax);
    opt.on(minimizer::Phase::force).add(ConstHydrostaticStress{0.05});
    opt.enable_timing(false);
    opt.build();
    opt.run(frame, 2);

    EXPECT_NEAR(frame.simbox.cell()(0, 1), h * 1.05, 1e-12);
    EXPECT_NEAR(frame.simbox.volume(), vol0 * (1.05 * 1.05 * 1.05), 1e-9);
}
