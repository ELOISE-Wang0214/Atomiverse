/// @file tests/cpu/multiscale/test_qmmm_electrostatic.cpp
/// @brief Unit: rewritten QmmmElectrostatic (esmm-qmmm-03).

#ifdef ATV_HAS_HF

#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/column.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/qm_mask.h"
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"
#include "atomiverse/cpu/scf/methods.h"
#include "atomiverse/cpu/semiempirical/am1_solver.h"

using namespace ATV;

namespace {

void insert_qm_z(cpu::Frame& f) {
    f.blocks.at("atoms").insert(
        "qm_z", cpu::ColumnSpec{ColDType::Int, 1, ColGroup::Core, 0, kColNone,
                           kDynamicFieldId});
}

cpu::Frame make_embed_h2o() {
    // Same layout as embed_h2o_mm_forces (3 QM + 2 MM charges as atoms)
    cpu::Frame f(5);
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_z)[0] = 0.1173;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_y)[1] = 0.7572;
    f.span(kAtom_z)[1] = -0.4692;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_y)[2] = -0.7572;
    f.span(kAtom_z)[2] = -0.4692;
    f.span(kAtom_type)[3] = 99;
    f.span(kAtom_z)[3] = 3.0;
    f.span(kAtom_charge)[3] = 0.5;
    f.span(kAtom_type)[4] = 99;
    f.span(kAtom_z)[4] = -3.0;
    f.span(kAtom_charge)[4] = -0.5;
    for (int i = 0; i < 5; ++i)
        f.span(kAtom_tag)[i] = i + 1;
    f.box = cpu::Box::cubic(20.0);
    insert_qm_z(f);
    auto qm = f.span(cpu::kAtomQmZ);
    qm[0] = 8;
    qm[1] = 1;
    qm[2] = 1;
    qm[3] = 0;
    qm[4] = 0;
    return f;
}

} // namespace

TEST(QmmmElectrostatic, Am1EngineThrowsNamed) {
    auto f = make_embed_h2o();
    // Empty LJ table (no cross terms)
    cpu::QmmmElectrostatic qmmm{{}, 10.0, 20.0};
    qmmm.set_engine(cpu::semiempirical::AM1{});
    cpu::Context ctx;
    cpu::BuildContext b{ctx, f};
    EXPECT_THROW(
        {
            try {
                qmmm.build(b);
            } catch (const std::runtime_error& e) {
                EXPECT_NE(std::string(e.what()).find("am1"), std::string::npos);
                throw;
            }
        },
        std::runtime_error);
}

TEST(QmmmElectrostatic, RhfRunsAndAccumulatesEnergy) {
    auto f = make_embed_h2o();
    f.meta.get(kEnergy) = 1.0; // pre-seed
    f.meta.get(kEnergyElectronic) = 0.5;

    cpu::QmmmElectrostatic qmmm{{}, 10.0, 20.0};
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-10;
    rhf.max_iter = 80;
    qmmm.set_engine(std::move(rhf));

    cpu::Context ctx;
    cpu::BuildContext b{ctx, f};
    qmmm.build(b);
    cpu::LaunchContext l{ctx, f, 0, 0};
    l.forces_requested = false;
    qmmm.launch(l);

    // Pre-seed preserved: E = 1.0 + E_qm (E_qm ~ −2000 eV)
    EXPECT_NEAR(f.meta.get(kEnergy), 1.0 + qmmm.last_qm_energy(), 1e-9);
    EXPECT_NEAR(f.meta.get(kEnergyElectronic), 0.5 + qmmm.last_qm_energy(),
                1e-9);
    EXPECT_TRUE(std::isfinite(qmmm.last_qm_energy()));
    EXPECT_LT(qmmm.last_qm_energy(), 0.0);
}

TEST(QmmmElectrostatic, ForcesRequestedFalseLeavesForces) {
    auto f = make_embed_h2o();
    f.span(kAtom_fx)[0] = 42.0;
    cpu::QmmmElectrostatic qmmm{{}, 10.0, 20.0};
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-8;
    rhf.max_iter = 40;
    qmmm.set_engine(std::move(rhf));
    cpu::Context ctx;
    cpu::BuildContext b{ctx, f};
    qmmm.build(b);
    cpu::LaunchContext l{ctx, f, 0, 0};
    l.forces_requested = false;
    qmmm.launch(l);
    EXPECT_DOUBLE_EQ(f.span(kAtom_fx)[0], 42.0);
}

TEST(QmmmElectrostatic, EngineSwapRhfRks) {
    auto f = make_embed_h2o();
    cpu::Context ctx;
    {
        cpu::QmmmElectrostatic qmmm{{}, 10.0, 20.0};
        cpu::scf::RHF rhf;
        rhf.basis = cpu::scf::BasisSet::STO3G;
        rhf.tol = 1e-8;
        rhf.max_iter = 40;
        qmmm.set_engine(std::move(rhf));
        cpu::BuildContext b{ctx, f};
        qmmm.build(b);
        cpu::LaunchContext l{ctx, f, 0, 0};
        l.forces_requested = false;
        qmmm.launch(l);
        EXPECT_TRUE(std::isfinite(qmmm.last_qm_energy()));
    }
    {
        auto f2 = make_embed_h2o();
        cpu::QmmmElectrostatic qmmm{{}, 10.0, 20.0};
        cpu::scf::RKS rks;
        rks.basis = cpu::scf::BasisSet::STO3G;
        rks.tol = 1e-6;
        rks.max_iter = 40;
        qmmm.set_engine(std::move(rks));
        cpu::BuildContext b{ctx, f2};
        qmmm.build(b);
        cpu::LaunchContext l{ctx, f2, 0, 0};
        l.forces_requested = false;
        qmmm.launch(l);
        EXPECT_TRUE(std::isfinite(qmmm.last_qm_energy()));
    }
}

#else
TEST(QmmmElectrostatic, SkippedWithoutHF) { GTEST_SKIP() << "ATV_HAS_HF off"; }
#endif
