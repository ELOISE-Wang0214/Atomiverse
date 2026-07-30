/// @file tests/cpu/scf/test_mp2.cpp
/// @brief Unit checks for RHF::mp2 / UHF::mp2 (no public internal::IntegralWorkspace).

#ifdef ATV_HAS_HF

#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

namespace {

cpu::Frame make_he() {
    cpu::Frame f(1);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[0] = 2;
    f.span(kAtom_mass)[0] = 4.0026;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_mass)[1] = 1.008;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

cpu::Frame make_h() {
    cpu::Frame f(1);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

} // namespace

TEST(MP2, HeSto3gZeroCorrelation) {
    // He/STO-3G: no virtuals → E_corr == 0
    auto frame = make_he();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-12;
    auto r = rhf.mp2(frame);
    EXPECT_TRUE(rhf.converged());
    EXPECT_NEAR(r.e_corr, 0.0, 1e-14);
    EXPECT_NEAR(r.e_tot, rhf.energy(), 1e-14);
}

TEST(MP2, H2Sto3gNegativeCorrelation) {
    auto frame = make_h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-12;
    auto r = rhf.mp2(frame);
    EXPECT_TRUE(rhf.converged());
    EXPECT_LT(r.e_corr, 0.0);
    EXPECT_NEAR(r.e_tot, rhf.energy() + r.e_corr, 1e-12);
    // Unscaled channels sum to plain MP2 (default scales 1/1).
    EXPECT_NEAR(r.e_corr, r.e_corr_os + r.e_corr_ss, 1e-12);
}

TEST(MP2, H2ScsScalesOsSs) {
    // Unit-only: SCS scales reweight OS/SS; e_corr = 1.2*OS + (1/3)*SS.
    auto frame = make_h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-12;
    rhf.mp2_os_scale = 1.2;
    rhf.mp2_ss_scale = 1.0 / 3.0;
    auto r = rhf.mp2(frame);
    EXPECT_TRUE(rhf.converged());
    EXPECT_LT(r.e_corr, 0.0);
    EXPECT_NEAR(r.e_corr, 1.2 * r.e_corr_os + (1.0 / 3.0) * r.e_corr_ss, 1e-12);
    // Two-electron H2: same-spin channel vanishes (only one occupied MO).
    EXPECT_NEAR(r.e_corr_ss, 0.0, 1e-14);
}

TEST(UMP2, HAtomSto3gZeroCorrelation) {
    // H/STO-3G UHF(1,0): single α electron, no doubles → E_corr == 0
    auto frame = make_h();
    cpu::scf::UHF uhf(1, 0);
    uhf.basis = cpu::scf::BasisSet::STO3G;
    uhf.tol = 1e-12;
    auto r = uhf.mp2(frame);
    EXPECT_TRUE(uhf.converged());
    EXPECT_NEAR(r.e_corr, 0.0, 1e-14);
    EXPECT_NEAR(r.e_tot, uhf.energy(), 1e-14);
    EXPECT_FALSE(r.ri);
}

TEST(UMP2, H2ClosedShellMatchesRMP2) {
    // Equilibrium H2: UHF(1,1) recovers the RHF solution → UMP2 == RMP2.
    auto frame_r = make_h2();
    auto frame_u = make_h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-12;
    const auto rmp2 = rhf.mp2(frame_r);

    cpu::scf::UHF uhf(1, 1);
    uhf.basis = cpu::scf::BasisSet::STO3G;
    uhf.tol = 1e-12;
    const auto ump2 = uhf.mp2(frame_u);

    EXPECT_TRUE(uhf.converged());
    EXPECT_NEAR(ump2.e_corr, rmp2.e_corr, 1e-9);
    EXPECT_NEAR(ump2.e_tot, rmp2.e_tot, 1e-9);
}

TEST(MP2, FrozenCoreRaisesEcorrVsFull) {
    // H2/STO-3G has 1 occupied spatial MO. Freezing it must zero E_corr;
    // unfrozen E_corr is strictly negative. Unit-only (no e2e pin).
    auto frame_full = make_h2();
    auto frame_fc = make_h2();
    cpu::scf::RHF full;
    full.basis = cpu::scf::BasisSet::STO3G;
    full.tol = 1e-12;
    const auto r_full = full.mp2(frame_full);
    EXPECT_TRUE(full.converged());
    EXPECT_LT(r_full.e_corr, 0.0);
    EXPECT_EQ(r_full.frozen, 0);

    cpu::scf::RHF fc;
    fc.basis = cpu::scf::BasisSet::STO3G;
    fc.tol = 1e-12;
    fc.frozen = 1;
    const auto r_fc = fc.mp2(frame_fc);
    EXPECT_TRUE(fc.converged());
    EXPECT_NEAR(r_fc.e_corr, 0.0, 1e-14);
    EXPECT_EQ(r_fc.frozen, 1);
    EXPECT_NEAR(r_fc.e_tot, fc.energy(), 1e-14);
}

TEST(MP2, FrozenCoreRejectsNegative) {
    auto frame = make_h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-12;
    rhf.frozen = -1;
    EXPECT_THROW(rhf.mp2(frame), std::invalid_argument);
}

TEST(ROMP2, LiSto3gNegativeCorrelation) {
    // Li atom doublet: ROHF(2,1)/STO-3G has virtuals → ROMP2 e_corr < 0.
    // Unit-only qualitative check (no literature pin; e2e is romp2-o2-sto3g).
    cpu::Frame frame(1);
    frame.span(kAtom_x)[0] = 0.0;
    frame.span(kAtom_y)[0] = 0.0;
    frame.span(kAtom_z)[0] = 0.0;
    frame.span(kAtom_type)[0] = 3; // Li
    frame.span(kAtom_mass)[0] = 6.94;
    frame.box = cpu::Box::cubic(20.0);

    cpu::scf::ROHF rohf(2, 1);
    rohf.basis = cpu::scf::BasisSet::STO3G;
    rohf.tol = 1e-12;
    const auto r = rohf.mp2(frame);
    EXPECT_TRUE(rohf.converged());
    EXPECT_LT(r.e_corr, 0.0);
    EXPECT_NEAR(r.e_tot, rohf.energy() + r.e_corr, 1e-12);
    EXPECT_FALSE(r.ri);
}

#endif
