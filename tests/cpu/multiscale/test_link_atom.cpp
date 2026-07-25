/// @file tests/cpu/multiscale/test_link_atom.cpp
/// @brief Unit tests for H link-atom placement and force projection.
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/link_atom.h"
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"

using namespace ATV;

TEST(LinkAtom, PlacementFormula) {
    double hx, hy, hz;
    // Q at origin, M at (2,0,0), g=0.5 → H at (1,0,0)
    cpu::place_link_hydrogen(0, 0, 0, 2, 0, 0, 0.5, hx, hy, hz);
    EXPECT_NEAR(hx, 1.0, 1e-14);
    EXPECT_NEAR(hy, 0.0, 1e-14);
    EXPECT_NEAR(hz, 0.0, 1e-14);
    // g=0.709 → H at 0.709 * r_QM along the bond
    cpu::place_link_hydrogen(1, 2, 3, 4, 2, 3, 0.709, hx, hy, hz);
    EXPECT_NEAR(hx, 1.0 + 0.709 * 3.0, 1e-12);
    EXPECT_NEAR(hy, 2.0, 1e-12);
    EXPECT_NEAR(hz, 3.0, 1e-12);
}

TEST(LinkAtom, ForceProjectionPartition) {
    double fqx = 0, fqy = 0, fqz = 0;
    double fmx = 0, fmy = 0, fmz = 0;
    const double g = 0.7;
    const double fhx = 10.0, fhy = -4.0, fhz = 2.0;
    cpu::project_link_force(g, fhx, fhy, fhz, fqx, fqy, fqz, fmx, fmy, fmz);
    EXPECT_NEAR(fqx, (1.0 - g) * fhx, 1e-14);
    EXPECT_NEAR(fqy, (1.0 - g) * fhy, 1e-14);
    EXPECT_NEAR(fqz, (1.0 - g) * fhz, 1e-14);
    EXPECT_NEAR(fmx, g * fhx, 1e-14);
    EXPECT_NEAR(fmy, g * fhy, 1e-14);
    EXPECT_NEAR(fmz, g * fhz, 1e-14);
    // Total projected force equals F_H (translation invariance of the map)
    EXPECT_NEAR(fqx + fmx, fhx, 1e-14);
    EXPECT_NEAR(fqy + fmy, fhy, 1e-14);
    EXPECT_NEAR(fqz + fmz, fhz, 1e-14);
}

TEST(LinkAtom, RejectsIllegalLinks) {
    EXPECT_THROW((cpu::QmmmElectrostatic{
                     {1, 2}, 0.0, 0.01, 3.4, 10.0, {{1, 1}} /* same tags */}),
                 std::invalid_argument);
    EXPECT_THROW(
        (cpu::QmmmElectrostatic{
            {1, 2}, 0.0, 0.01, 3.4, 10.0, {{1, 2}} /* mm_tag is QM */}),
        std::invalid_argument);
    EXPECT_THROW((cpu::QmmmElectrostatic{{1}, 0.0, 0.01, 3.4, 10.0, {{9, 3}}}),
                 std::invalid_argument);
}

TEST(LinkAtom, QmmmPlacesLinkHydrogen) {
    // Closed-shell cut: H(QM host) + link H → H2 (2 e⁻) for RHF.
    // MM partner C + free MM site.
    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 1.5;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_x)[2] = 4.0;
    f.span(kAtom_y)[2] = 0;
    f.span(kAtom_z)[2] = 0;
    f.span(kAtom_type)[0] = 1; // H (QM host) + link H → closed shell
    f.span(kAtom_type)[1] = 6; // C (MM link partner)
    f.span(kAtom_type)[2] = 1; // H-like MM site
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.span(kAtom_tag)[2] = 3;
    f.simbox = cpu::Box::cubic(20);
    f.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    cpu::QmmmElectrostatic qm({1}, /*mm_q=*/0.0, 0.01, 3.4, 10.0, {{1, 2}},
                              0.709);
    qm.build({ctx, f});
    cpu::LaunchContext lc{ctx, f, 0, 0};
    lc.forces_requested = false;
    qm.launch(lc);
    EXPECT_EQ(qm.last_n_links(), 1);
    EXPECT_TRUE(std::isfinite(f.meta.get(kEnergy)));
    EXPECT_TRUE(std::isfinite(qm.last_qm_energy()));
}
