/// @file tests/cpu/scf/test_ri_mp2.cpp
/// @brief RI-MP2 vs canonical MP2 on H2/STO-3G.
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

TEST(RiMp2, H2Sto3gCloseToCanonical) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[0] = f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[0] = f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.simbox = cpu::Box::cubic(20);

    cpu::scf::RHF can;
    can.basis = cpu::scf::BasisSet::STO3G;
    can.tol = 1e-10;
    auto rc = can.mp2(f);
    EXPECT_FALSE(rc.ri);
    EXPECT_LT(rc.e_corr, 0.0);

    cpu::scf::RHF ri;
    ri.basis = cpu::scf::BasisSet::STO3G;
    ri.tol = 1e-10;
    ri.ri = cpu::scf::RI::JK;
    ri.aux_basis = cpu::scf::AuxBasis::Def2JKFIT;
    auto rr = ri.mp2(f);
    EXPECT_TRUE(rr.ri);
    EXPECT_LT(rr.e_corr, 0.0);
    EXPECT_NEAR(rr.e_corr, rc.e_corr, 1e-4); // Ha; DF error small on H2
}

// Equivalence guard for the RI-MP2 GEMM rewrite (spec scf-review-07-posthf-perf,
// item 18): ri_mp2_energy's per-P dense-M reconstruction + scalar inner-P
// accumulation are replaced by a single BLAS GEMM (ia|jb)=B_ov·B_ovᵀ. The
// fitting metric (Coulomb V=(P|Q), Cholesky-whitened W=L⁻¹ ⇒ V⁻¹) and the
// df.b_tensor() output are FROZEN — the rewrite only reassociates the
// contraction, so e_corr must be unchanged to within the FP reassociation bound
// N·ε ≈ 1e-10 Ha (NOT machine precision). Golden captured from the pre-GEMM
// scalar implementation 2026-07-24 (GCC 14.3, def2-JKFIT aux, double).
TEST(RiMp2, RiMp2GoldenUnchanged) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[0] = f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[0] = f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.simbox = cpu::Box::cubic(20);

    cpu::scf::RHF ri;
    ri.basis = cpu::scf::BasisSet::STO3G;
    ri.tol = 1e-10;
    ri.ri = cpu::scf::RI::JK;
    ri.aux_basis = cpu::scf::AuxBasis::Def2JKFIT;
    auto rr = ri.mp2(f);

    // Pre-GEMM scalar golden (Ha). H2 is closed-shell 1-occupied, so the
    // same-spin channel is exactly zero and e_corr == e_corr_os.
    constexpr double kGoldenECorr = -1.3126142949017863e-02;
    EXPECT_TRUE(rr.ri);
    EXPECT_NEAR(rr.e_corr, kGoldenECorr, 1e-10)
        << "RI-MP2 e_corr drifted from the pre-GEMM scalar golden";
    EXPECT_NEAR(rr.e_corr_os, kGoldenECorr, 1e-10);
    EXPECT_NEAR(rr.e_corr_ss, 0.0, 1e-12);
}
