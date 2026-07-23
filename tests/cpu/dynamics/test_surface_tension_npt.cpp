/// @file tests/cpu/dynamics/test_surface_tension_npt.cpp
/// @brief Unit tests for SurfaceTensionNPT / NoseHooverNPT couple-xy path.

#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/integrator.h"
#include "atomiverse/cpu/dynamics/init.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/surface_tension_npt.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

namespace {

cpu::Frame make_lj_fcc(int n_cells, double a_lat, double mass) {
    const int n = n_cells * n_cells * n_cells * 4;
    cpu::Frame f(n);
    const double basis[4][3] = {
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
    int idx = 0;
    for (int ix = 0; ix < n_cells; ++ix)
        for (int iy = 0; iy < n_cells; ++iy)
            for (int iz = 0; iz < n_cells; ++iz)
                for (int b = 0; b < 4; ++b) {
                    f.span(kAtom_x)[idx] = (ix + basis[b][0]) * a_lat;
                    f.span(kAtom_y)[idx] = (iy + basis[b][1]) * a_lat;
                    f.span(kAtom_z)[idx] = (iz + basis[b][2]) * a_lat;
                    f.span(kAtom_mass)[idx] = mass;
                    f.span(kAtom_type)[idx] = 18;
                    ++idx;
                }
    f.simbox = cpu::Box::cubic(n_cells * a_lat);
    return f;
}

/// Instantaneous planar γ [eV/Å²] from frame stress + kinetic (diagnostic).
double gamma_inst(const cpu::Frame& f) {
    const double V = f.simbox.volume();
    if (!(V > 0.0))
        return 0.0;
    const auto& s = f.meta.get(kStress);
    const auto& sk = f.meta.get(kStressKinetic);
    // P_αα = (virial_α + kinetic_α) / V — kinetic is Σ m v_α² = 2 KE_α
    const double Pxx = (s[0] + sk[0]) / V;
    const double Pyy = (s[1] + sk[1]) / V;
    const double Pzz = (s[2] + sk[2]) / V;
    const double Lz = f.simbox.length(2);
    return 0.5 * Lz * (Pzz - 0.5 * (Pxx + Pyy));
}

} // namespace

TEST(SurfaceTensionNPT, IsIntegratorApiNoConfig) {
    cpu::md::SurfaceTensionNPT st{/*T=*/100.0, /*tau_T=*/50.0,
                                  /*P_n=*/0.0, /*gamma=*/0.001,
                                  /*tau_P=*/500.0, /*dt=*/1.0};
    EXPECT_NE(dynamic_cast<const cpu::Integrator*>(&st), nullptr);
    EXPECT_TRUE(st.requires_system_stress());
    EXPECT_EQ(st.name(), "surface_tension_npt");
    EXPECT_TRUE(st.couple_xy());
    EXPECT_DOUBLE_EQ(st.surface_tension(), 0.001);
    EXPECT_EQ(st.n_surfaces(), 2);
}

TEST(SurfaceTensionNPT, CoupleXyKeepsLxEqualsLy) {
    const double eps = 0.0104;
    const double sigma = 3.4;
    const double cutoff = 8.5;
    const double a_lat = 5.1;
    auto frame = make_lj_fcc(/*n_cells=*/3, a_lat, 39.948);
    const double L0 = frame.simbox.length(0);

    cpu::md::MD sim;
    sim.add_potential(cpu::PairLJCut{eps, sigma, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.5});
    sim.set_integrator(cpu::md::SurfaceTensionNPT{
        /*T=*/100.0, /*tau_T=*/50.0, /*P_n=*/0.0, /*gamma=*/0.001,
        /*tau_P=*/500.0, /*dt=*/1.0, /*n_surfaces=*/2});
    sim.add(cpu::md::InitVelocity{100.0, 777u});
    sim.enable_timing(false);
    sim.build();
    sim.run(frame, 80);

    const double Lx = frame.simbox.length(0);
    const double Ly = frame.simbox.length(1);
    const double Lz = frame.simbox.length(2);
    EXPECT_TRUE(std::isfinite(Lx) && std::isfinite(Ly) && std::isfinite(Lz));
    EXPECT_NEAR(Lx, Ly, 1e-10);
    // Box should move from compressed start (not frozen).
    EXPECT_GT(std::fabs(Lx - L0) + std::fabs(Lz - L0), 1e-6);
    EXPECT_GT(frame.simbox.volume(), 0.0);

    const double g = gamma_inst(frame);
    EXPECT_TRUE(std::isfinite(g));
}

TEST(NoseHooverNPT, SetSurfaceTensionEnablesCoupleXy) {
    cpu::md::NoseHooverNPT npt{100.0, 50.0, 0.0, 500.0, 1.0};
    EXPECT_FALSE(npt.couple_xy());
    npt.set_surface_tension(0.002, 2);
    EXPECT_TRUE(npt.couple_xy());
    EXPECT_DOUBLE_EQ(npt.surface_tension(), 0.002);
}
