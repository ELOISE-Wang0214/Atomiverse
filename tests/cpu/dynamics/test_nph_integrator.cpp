/// @file tests/cpu/dynamics/test_nph_integrator.cpp
/// @brief Unit tests for cpu::md::NphIntegrator (thermostat-free NPH).

#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/integrator.h"
#include "atomiverse/cpu/dynamics/init.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/nph_integrator.h"
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

/// @brief Public PE+KE+PV piece of the NPH enthalpy (barostat KE lives inside
///        the integrator's kConserved).
double enthalpy_from_frame(const cpu::Frame& f, double P_target) {
    return f.meta.get(kKinetic) + f.meta.get(kEnergy)
           + P_target * f.box.volume();
}

} // namespace

TEST(NphIntegrator, IsIntegratorRequiresStressNoFactory) {
    cpu::md::NphIntegrator nph{/*P=*/0.0, /*tau_P=*/200.0, /*dt=*/2.0};
    EXPECT_NE(dynamic_cast<const cpu::Integrator*>(&nph), nullptr);
    EXPECT_TRUE(nph.requires_system_stress());
    EXPECT_EQ(nph.name(), "nph_integrator");
}

TEST(NphIntegrator, EnthalpyDefinitionSelfConsistent) {
    // After one short NPH segment, frame.meta[kConserved] (extended H with
    // barostat KE) is finite and the PE+KE+PV piece is self-consistent.
    const double eps = 120.0 * units::kB_ev_per_k;
    const double sigma = 3.4;
    const double cutoff = 2.5 * sigma;
    const double density = 0.02071;
    const double T_init = 94.4;
    const double P = 0.0;
    const double dt = 2.0;

    auto frame = make_lj_sc(4, density, 39.948);
    const double V0 = frame.box.volume();

    cpu::md::MD sim;
    sim.set_evaluator(cpu::PairLJCut{eps, sigma, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.5});
    sim.set_integrator(cpu::md::NphIntegrator{P, /*tau_P=*/200.0, dt});
    sim.set_observer(cpu::md::InitVelocity{T_init, 13u});
    sim.enable_timing(false);
    sim.build();
    sim.run(frame, 1);

    const double H_pub = enthalpy_from_frame(frame, P);
    EXPECT_TRUE(std::isfinite(H_pub));
    EXPECT_TRUE(std::isfinite(frame.meta.get(kConserved)));
    // With P=0, PE+KE equals the enthalpy piece; conserved adds ½Wε̇² ≥ that.
    EXPECT_GE(frame.meta.get(kConserved) + 1e-12, H_pub - 1e-6);

    sim.run(frame, 15);
    const double V1 = frame.box.volume();
    EXPECT_GT(V1, 0.0);
    EXPECT_TRUE(std::isfinite(V1));
    EXPECT_NE(V1, V0);
    EXPECT_TRUE(std::isfinite(frame.meta.get(kConserved)));
    EXPECT_TRUE(std::isfinite(enthalpy_from_frame(frame, P)));
}

TEST(NphIntegrator, VolumeStaysPositive) {
    const double eps = 120.0 * units::kB_ev_per_k;
    const double sigma = 3.4;
    const double cutoff = 2.5 * sigma;

    auto frame = make_lj_sc(4, 0.02071, 39.948);
    const double V0 = frame.box.volume();

    cpu::md::MD sim;
    sim.set_evaluator(cpu::PairLJCut{eps, sigma, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.5});
    sim.set_integrator(cpu::md::NphIntegrator{0.0, 200.0, 2.0});
    sim.set_observer(cpu::md::InitVelocity{94.4, 17u});
    sim.enable_timing(false);
    sim.build();
    sim.run(frame, 30);

    EXPECT_GT(frame.box.volume(), 0.05 * V0);
    EXPECT_LT(frame.box.volume(), 20.0 * V0);
}
