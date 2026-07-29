/// @file tests/cpu/multiscale/test_link_atom.cpp
/// @brief Unit tests for H link-atom placement and force projection.
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/column.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/link_atom.h"
#include "atomiverse/cpu/multiscale/qm_mask.h"
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

TEST(LinkAtom, PlacementFormula) {
    double hx, hy, hz;
    cpu::place_link_hydrogen(0, 0, 0, 2, 0, 0, 0.5, hx, hy, hz);
    EXPECT_NEAR(hx, 1.0, 1e-14);
    EXPECT_NEAR(hy, 0.0, 1e-14);
    EXPECT_NEAR(hz, 0.0, 1e-14);
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
    EXPECT_NEAR(fqx + fmx, fhx, 1e-14);
    EXPECT_NEAR(fqy + fmy, fhy, 1e-14);
    EXPECT_NEAR(fqz + fmz, fhz, 1e-14);
}

#ifdef ATV_HAS_HF
TEST(LinkAtom, RejectsIllegalLinksAtBuild) {
    cpu::Frame f(3);
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 6;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.span(kAtom_tag)[2] = 3;
    f.box = cpu::Box::cubic(20);
    f.blocks.at("atoms").insert(
        "qm_z", cpu::ColumnSpec{ColDType::Int, 1, ColGroup::Core, 0, kColNone,
                           kDynamicFieldId});
    auto qm = f.span(cpu::kAtomQmZ);
    qm[0] = 1;
    qm[1] = 0;
    qm[2] = 0;

    cpu::Context ctx;
    {
        cpu::QmmmElectrostatic q{{{1, 1, 0.01, 3.4}}, 10.0, 12.0};
        q.links() = {{1, 1}}; // same tags
        q.set_engine(cpu::scf::RHF{});
        EXPECT_THROW(q.build({ctx, f}), std::invalid_argument);
    }
    {
        cpu::QmmmElectrostatic q{{{1, 1, 0.01, 3.4}}, 10.0, 12.0};
        q.links() = {{1, 2}};
        // mark mm as QM
        qm[1] = 6;
        q.set_engine(cpu::scf::RHF{});
        EXPECT_THROW(q.build({ctx, f}), std::invalid_argument);
        qm[1] = 0;
    }
}

TEST(LinkAtom, QmmmPlacesLinkHydrogen) {
    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_x)[1] = 1.5;
    f.span(kAtom_x)[2] = 4.0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 6;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.span(kAtom_tag)[2] = 3;
    f.span(kAtom_charge)[2] = 0.0;
    f.box = cpu::Box::cubic(20);
    f.meta.get(kEnergy) = 0;
    f.blocks.at("atoms").insert(
        "qm_z", cpu::ColumnSpec{ColDType::Int, 1, ColGroup::Core, 0, kColNone,
                           kDynamicFieldId});
    auto qmcol = f.span(cpu::kAtomQmZ);
    qmcol[0] = 1; // H QM
    qmcol[1] = 0;
    qmcol[2] = 0;

    cpu::Context ctx;
    cpu::QmmmElectrostatic qm{{{1, 1, 0.01, 3.4}}, 10.0, 12.0};
    qm.links() = {{1, 2}};
    qm.link_scale() = 0.709;
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-8;
    rhf.max_iter = 40;
    qm.set_engine(std::move(rhf));
    qm.build({ctx, f});
    cpu::LaunchContext lc{ctx, f, 0, 0};
    lc.forces_requested = false;
    qm.launch(lc);
    EXPECT_EQ(qm.last_n_links(), 1);
    EXPECT_TRUE(std::isfinite(f.meta.get(kEnergy)));
    EXPECT_TRUE(std::isfinite(qm.last_qm_energy()));
}
#endif
