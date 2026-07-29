/// @file tests/cpu/component/test_qm_engine.cpp
/// @brief Unit: QmEngine / EmbeddingEngine protocol (esmm-qmmm-02).

#ifdef ATV_HAS_HF

#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/component/qm_engine.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"
#include "atomiverse/cpu/semiempirical/am1_solver.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;

namespace {

cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_x)[1] = 0.74;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_z)[0] = 0.1173;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_y)[1] = 0.7572;
    f.span(kAtom_z)[1] = -0.4692;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_y)[2] = -0.7572;
    f.span(kAtom_z)[2] = -0.4692;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

} // namespace

TEST(QmEngine, GaussianMethodsExposeEmbedding) {
    cpu::scf::RHF rhf;
    cpu::scf::UHF uhf;
    cpu::scf::RKS rks;
    cpu::QmEngine* engines[] = {&rhf, &uhf, &rks};
    for (cpu::QmEngine* e : engines) {
        ASSERT_NE(e->as_embedding(), nullptr);
        EXPECT_FALSE(e->qm_engine_name().empty());
    }
}

TEST(QmEngine, Am1AndDftb2NoEmbedding) {
    cpu::semiempirical::AM1 am1;
    cpu::semiempirical::DFTB2 dftb;
    cpu::QmEngine* engines[] = {&am1, &dftb};
    for (cpu::QmEngine* e : engines) {
        EXPECT_EQ(e->as_embedding(), nullptr);
        EXPECT_FALSE(e->qm_engine_name().empty());
    }
}

TEST(QmEngine, ErasedRhfEnergyMatchesDirectBitForBit) {
    auto frame = make_h2();
    cpu::scf::RHF direct;
    direct.basis = cpu::scf::BasisSet::STO3G;
    direct.tol = 1e-12;
    direct.max_iter = 80;
    direct.run(frame);
    const double e_direct = direct.energy();

    cpu::scf::RHF eng;
    eng.basis = cpu::scf::BasisSet::STO3G;
    eng.tol = 1e-12;
    eng.max_iter = 80;
    cpu::QmEngine* pe = &eng;
    auto frame2 = make_h2();
    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame2};
    pe->prepare(bctx);
    // Zero energy before evaluate (private-frame contract).
    frame2.meta.get(kEnergy) = 0.0;
    frame2.meta.get(kEnergyElectronic) = 0.0;
    cpu::LaunchContext lctx{ctx, frame2, 0, 0};
    lctx.forces_requested = false;
    pe->evaluate(lctx);
    const double e_erased = frame2.meta.get(kEnergy) / units::kHartreeInEv;
    EXPECT_EQ(e_erased, e_direct);
}

TEST(QmEngine, ForcesRequestedFalseLeavesForcesUntouched) {
    auto frame = make_h2();
    // Poison force channels
    frame.span(kAtom_fx)[0] = 1.23;
    frame.span(kAtom_fy)[0] = 4.56;
    frame.span(kAtom_fz)[0] = 7.89;
    frame.span(kAtom_fx)[1] = -1.0;
    frame.span(kAtom_fy)[1] = -2.0;
    frame.span(kAtom_fz)[1] = -3.0;

    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-10;
    rhf.max_iter = 80;
    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};
    rhf.prepare(bctx);
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    lctx.forces_requested = false;
    rhf.evaluate(lctx);

    EXPECT_DOUBLE_EQ(frame.span(kAtom_fx)[0], 1.23);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_fy)[0], 4.56);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_fz)[0], 7.89);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_fx)[1], -1.0);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_fy)[1], -2.0);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_fz)[1], -3.0);
}

TEST(QmEngine, ResetHistoryClearsWarmStart) {
    auto frame = make_h2o();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-10;
    rhf.max_iter = 100;
    // Enable ASPC so a second evaluate without reset would warm-start.
    rhf.predictor(3);

    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};
    rhf.prepare(bctx);
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    lctx.forces_requested = false;

    frame.meta.get(kEnergy) = 0.0;
    rhf.evaluate(lctx);
    const int cold_iters = rhf.iterations();
    ASSERT_GT(cold_iters, 0);

    // Second evaluate without reset: warm ASPC path → fewer iterations.
    frame.meta.get(kEnergy) = 0.0;
    rhf.evaluate(lctx);
    const int warm_iters = rhf.iterations();

    rhf.reset_history();
    frame.meta.get(kEnergy) = 0.0;
    rhf.evaluate(lctx);
    const int after_reset = rhf.iterations();

    // After reset, iteration count matches cold (or is within 1 of it —
    // ASPC engagement can still change the first-iter residual slightly).
    EXPECT_GE(after_reset, warm_iters);
    EXPECT_NEAR(static_cast<double>(after_reset),
                static_cast<double>(cold_iters), 2.0)
        << "cold=" << cold_iters << " warm=" << warm_iters
        << " after_reset=" << after_reset;
}

TEST(QmEngine, EmbeddingSitesRoundTrip) {
    cpu::scf::RHF rhf;
    cpu::EmbeddingEngine* emb = rhf.as_embedding();
    ASSERT_NE(emb, nullptr);
    emb->set_site_self_energy(false);
    emb->set_embedding_sites({{0.0, 0.0, 3.0, 0.5}});
    EXPECT_EQ(rhf.point_charges.size(), 1u);
    EXPECT_FALSE(rhf.embed_self_energy);
}

#else
TEST(QmEngine, SkippedWithoutHF) { GTEST_SKIP() << "ATV_HAS_HF off"; }
#endif
