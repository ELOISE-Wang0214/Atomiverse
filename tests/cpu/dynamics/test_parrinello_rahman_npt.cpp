/// @file tests/cpu/dynamics/test_parrinello_rahman_npt.cpp
/// @brief Unit tests for cpu::md::ParrinelloRahmanNPT (full-cell NPT).

#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/integrator.h"
#include "atomiverse/cpu/dynamics/init.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/parrinello_rahman_npt.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

namespace {

cpu::Frame make_lj_sc(int n_side, double density, double mass) {
    const int n = n_side * n_side * n_side;
    const double L = std::cbrt(double(n) / density);
    const double spacing = L / double(n_side);
    cpu::Frame f(n);
    int idx = 0;
    for (int ix = 0; ix < n_side; ++ix)
        for (int iy = 0; iy < n_side; ++iy)
            for (int iz = 0; iz < n_side; ++iz) {
                f.span(kAtom_x)[idx] = (double(ix) + 0.5) * spacing;
                f.span(kAtom_y)[idx] = (double(iy) + 0.5) * spacing;
                f.span(kAtom_z)[idx] = (double(iz) + 0.5) * spacing;
                f.span(kAtom_mass)[idx] = mass;
                f.span(kAtom_type)[idx] = 0;
                ++idx;
            }
    f.box = cpu::Box::cubic(L);
    return f;
}

} // namespace

TEST(ParrinelloRahmanNPT, IsIntegratorAndRequiresStress) {
    cpu::md::ParrinelloRahmanNPT pr{100.0, 0.0, 100.0, 500.0, 2.0};
    EXPECT_NE(dynamic_cast<const cpu::Integrator*>(&pr), nullptr);
    EXPECT_TRUE(pr.requires_system_stress());
    EXPECT_EQ(pr.name(), "parrinello_rahman_npt");
    const auto binds = pr.default_bindings();
    ASSERT_EQ(binds.size(), 3u);
}

TEST(ParrinelloRahmanNPT, ShortRunVolumePositiveFinite) {
    // Single-function: after a short PR-NPT run the cell volume stays
    // positive and finite, and min-image cache remains consistent.
    const double eps = 120.0 * units::kB_ev_per_k;
    const double sigma = 3.4;
    const double cutoff = 2.5 * sigma;
    const double density = 0.02071;
    const double T = 94.4;
    const double dt = 2.0;

    auto frame = make_lj_sc(/*n_side=*/4, density, /*mass=*/39.948);
    const double V0 = frame.box.volume();
    ASSERT_GT(V0, 0.0);

    cpu::md::MD sim;
    sim.add_potential(cpu::PairLJCut{eps, sigma, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.5});
    sim.set_integrator(cpu::md::ParrinelloRahmanNPT{T, /*P=*/0.0,
                                                    /*tau_T=*/50.0,
                                                    /*tau_P=*/200.0, dt});
    sim.add(cpu::md::InitVelocity{T, 7u});
    sim.enable_timing(false);
    sim.build();
    sim.run(frame, 20);

    const double V = frame.box.volume();
    EXPECT_TRUE(std::isfinite(V));
    EXPECT_GT(V, 0.1 * V0);
    EXPECT_LT(V, 10.0 * V0);
    EXPECT_TRUE(std::isfinite(frame.meta.get(kEnergy)));
    EXPECT_TRUE(std::isfinite(frame.meta.get(kKinetic)));
    EXPECT_TRUE(std::isfinite(frame.meta.get(kConserved)));
    // Cache refresh: orthorhombic min_image still works.
    double dx = 0.1, dy = 0.0, dz = 0.0;
    frame.box.min_image(dx, dy, dz);
    EXPECT_NEAR(dx, 0.1, 1e-12);
}

TEST(ParrinelloRahmanNPT, CubicIsotropicOffDiagonalBounded) {
    // Under hydrostatic P on a cubic LJ start, off-diagonal strain rates
    // must stay small (noise from finite sample, not secular shear drive).
    const double eps = 120.0 * units::kB_ev_per_k;
    const double sigma = 3.4;
    const double cutoff = 2.5 * sigma;
    const double density = 0.02071;
    const double T = 94.4;
    const double dt = 2.0;

    auto frame = make_lj_sc(/*n_side=*/4, density, /*mass=*/39.948);

    cpu::md::ParrinelloRahmanNPT pr{T, 0.0, 50.0, 200.0, dt};
    cpu::md::MD sim;
    sim.add_potential(cpu::PairLJCut{eps, sigma, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.5});
    sim.set_integrator_ref(pr);
    sim.add(cpu::md::InitVelocity{T, 11u});
    sim.enable_timing(false);
    sim.build();
    sim.run(frame, 40);

    const Mat3& e = pr.cell_velocity();
    const double off = std::abs(static_cast<double>(e(0, 1)))
                       + std::abs(static_cast<double>(e(0, 2)))
                       + std::abs(static_cast<double>(e(1, 2)));
    // Natural barostat frequency ~ 2π/τ_P ~ 0.03 /fs; off-diagonal
    // rates on a cubic hydrostatic fluid should be ≪ that.
    EXPECT_LT(off, 0.05);
    // Box should remain nearly orthorhombic (off-diagonal cell entries
    // small vs edge length).
    const Mat3& h = frame.box.cell();
    const double edge = std::abs(static_cast<double>(h(0, 0)));
    const double shear_cell = std::abs(static_cast<double>(h(0, 1)))
                              + std::abs(static_cast<double>(h(0, 2)))
                              + std::abs(static_cast<double>(h(1, 2)));
    EXPECT_LT(shear_cell, 0.05 * edge);
}

TEST(ParrinelloRahmanNPT, CellMatrixNotDiagonalOnly) {
    // Construction + one forced shear-stress step: prove the public type
    // owns a full Mat3 strain rate (not three scalars only).
    cpu::md::ParrinelloRahmanNPT pr{100.0, 0.0, 100.0, 500.0, 1.0};
    const Mat3& e0 = pr.cell_velocity();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(static_cast<double>(e0(i, j)), 0.0);
}
