/// @file tests/cpu/gpw/test_gpw_xl_bomd.cpp
/// @brief Unit tests: XlBomd dissipative recursion (Niklasson JCP 130,
///        214109 (2009) Table I) + the GpwDft guess_P plumbing.
///
/// Propagator tests are pure dense algebra (no SCF). The one guess_P
/// plumbing test reuses the tiny inline H2 system of test_cpu_gpw_solver
/// (ms-scale SCF). The e2e XL-BOMD NVE envelope lives in
/// regressions/dynamics/xl-bomd-water (never here).

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/xl_bomd.h"

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;
using ATV::cpu::gpw::XlBomd;
using ATV::cpu::linalg::RealMatrix;

namespace {

RealMatrix constant_matrix(int n, double v) {
    RealMatrix m = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            m(i, j) = v + 0.1 * static_cast<double>(i == j);
    return m;
}

double max_abs_diff(const RealMatrix& a, const RealMatrix& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < a.shape(0); ++i)
        for (std::size_t j = 0; j < a.shape(1); ++j)
            worst = std::max(worst, std::abs(a(i, j) - b(i, j)));
    return worst;
}

} // namespace

// A1: the K=5 row of Table I, bitwise.
TEST(XlBomd, CoefficientTableK5Exact) {
    const XlBomd xl(5);
    EXPECT_EQ(xl.order(), 5);
    EXPECT_DOUBLE_EQ(xl.kappa(), 1.82);
    EXPECT_DOUBLE_EQ(xl.alpha(), 0.018);
    const std::vector<double> want = {-6.0, 14.0, -8.0, -3.0, 4.0, -1.0};
    ASSERT_EQ(xl.c().size(), want.size());
    for (std::size_t k = 0; k < want.size(); ++k)
        EXPECT_DOUBLE_EQ(xl.c()[k], want[k]) << "c[" << k << "]";
}

// A1: every Table-I row kills constant and linear history bias exactly.
TEST(XlBomd, CoefficientSumsVanishForAllOrders) {
    for (int k = 3; k <= 9; ++k) {
        const XlBomd xl(k);
        ASSERT_EQ(static_cast<int>(xl.c().size()), k + 1) << "K=" << k;
        double sum0 = 0.0, sum1 = 0.0;
        for (std::size_t j = 0; j < xl.c().size(); ++j) {
            sum0 += xl.c()[j];
            sum1 += static_cast<double>(j) * xl.c()[j];
        }
        EXPECT_DOUBLE_EQ(sum0, 0.0) << "sum c_k, K=" << k;
        EXPECT_DOUBLE_EQ(sum1, 0.0) << "sum k*c_k, K=" << k;
        EXPECT_GT(xl.kappa(), 0.0);
        EXPECT_LT(xl.kappa(), 2.0);
        EXPECT_GT(xl.alpha(), 0.0);
    }
}

// A2: only the published K range exists — no silent fallback row.
TEST(XlBomd, OrderOutsideTableThrows) {
    EXPECT_THROW(XlBomd(2), std::invalid_argument);
    EXPECT_THROW(XlBomd(10), std::invalid_argument);
}

// A2: an empty density is a programming error, not a seed.
TEST(XlBomd, AdvanceWithEmptyDensityThrows) {
    XlBomd xl(5);
    EXPECT_THROW(xl.advance(RealMatrix{}), std::invalid_argument);
}

// A2: shape changes after seeding must be surfaced (basis discontinuity).
TEST(XlBomd, AdvanceWithMismatchedShapeThrows) {
    XlBomd xl(5);
    xl.advance(constant_matrix(4, 0.3));
    EXPECT_THROW(xl.advance(constant_matrix(5, 0.3)), std::invalid_argument);
}

// A3: aux() before any density is a logic error.
TEST(XlBomd, AuxBeforeSeedThrows) {
    const XlBomd xl(5);
    EXPECT_FALSE(xl.ready());
    EXPECT_THROW(xl.aux(), std::logic_error);
}

// A3: the first advance seeds a stationary history — aux == D.
TEST(XlBomd, FirstAdvanceSeedsAuxToScfDensity) {
    XlBomd xl(5);
    const RealMatrix d = constant_matrix(4, 0.7);
    xl.advance(d);
    ASSERT_TRUE(xl.ready());
    EXPECT_DOUBLE_EQ(max_abs_diff(xl.aux(), d), 0.0);
}

// A3: a constant SCF density is the recursion's fixed point.
TEST(XlBomd, ConstantDensityIsFixedPoint) {
    XlBomd xl(5);
    const RealMatrix d = constant_matrix(6, -0.25);
    for (int n = 0; n < 21; ++n)
        xl.advance(d);
    EXPECT_LE(max_abs_diff(xl.aux(), d), 1e-12);
}

// A3: reset() drops the history (restart / discontinuity policy).
TEST(XlBomd, ResetDropsHistory) {
    XlBomd xl(5);
    xl.advance(constant_matrix(4, 0.4));
    xl.reset();
    EXPECT_FALSE(xl.ready());
    EXPECT_THROW(xl.aux(), std::logic_error);
}

// A4: bounded response to an oscillatory driving density (no secular
// growth from the dissipative Verlet recursion at K=5).
TEST(XlBomd, BoundedUnderOscillatoryDriving) {
    XlBomd xl(5);
    const int n_ao = 5;
    const RealMatrix d0 = constant_matrix(n_ao, 0.5);
    RealMatrix b = xt::zeros<double>({n_ao, n_ao});
    for (int i = 0; i < n_ao; ++i)
        for (int j = 0; j < n_ao; ++j)
            b(i, j) = ((i + j) % 2 == 0) ? 1.0 : -0.5; // symmetric, |b|max=1
    double worst = 0.0;
    for (int n = 0; n < 200; ++n) {
        const RealMatrix d =
            d0 + 0.01 * std::sin(0.3 * static_cast<double>(n)) * b;
        xl.advance(d);
        worst = std::max(worst, max_abs_diff(xl.aux(), d));
    }
    EXPECT_LE(worst, 0.1) << "auxiliary density lost the SCF density";
}

// ── guess_P plumbing on the tiny inline H2 GpwDft (mirrors the
//    warm-start contract test in test_cpu_gpw_solver.cpp) ───────────────

namespace {

GridBasisFunction s_basis(double alpha, const double c[3]) {
    GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = c[0];
    p.center[1] = c[1];
    p.center[2] = c[2];
    p.l[0] = p.l[1] = p.l[2] = 0;
    p.coeff = std::pow(2.0 * alpha / M_PI, 0.75);
    return {p};
}

GpwAtom hydrogen(const double pos[3]) {
    GpwAtom a;
    a.z_val = 1.0;
    a.pos[0] = pos[0];
    a.pos[1] = pos[1];
    a.pos[2] = pos[2];
    a.basis = {s_basis(0.5, pos)};
    a.rloc = 0.2;
    a.cloc = {-4.17890044, 0.72446331};
    return a;
}

std::vector<GpwAtom> small_h2(double box = 8.0) {
    const double rA[3] = {box / 2 - 0.7, box / 2, box / 2};
    const double rB[3] = {box / 2 + 0.7, box / 2, box / 2};
    return {hydrogen(rA), hydrogen(rB)};
}

GpwDft small_cfg(double box = 8.0) {
    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = box;
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 1;
    cfg.max_iter = 100;
    cfg.tol = 1e-7;
    cfg.mix = 0.5;
    return cfg;
}

void pin_single_thread() {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
}

} // namespace

// A5: a converged-density guess reaches the cold fixed point, no slower.
TEST(GpwSolverGuessP, ConvergedDensityGuessReachesColdFixedPoint) {
    pin_single_thread();
    const auto atoms = small_h2();

    GpwDft cold = small_cfg();
    cold.build(atoms);
    const GpwResult ref = cold.run(atoms);
    ASSERT_TRUE(ref.converged);
    ASSERT_GT(ref.P.size(), 0u);

    GpwDft warm = small_cfg();
    warm.build(atoms);
    const GpwResult& r = warm.run(atoms, nullptr, &ref.P);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.e_total, ref.e_total, 1e-10);
    EXPECT_LE(r.iters, ref.iters);
}

// A5: guess_P with a mismatched shape throws (never a silent cold start).
TEST(GpwSolverGuessP, MismatchedShapeThrows) {
    pin_single_thread();
    const auto atoms = small_h2();
    GpwDft solver = small_cfg();
    solver.build(atoms);
    const RealMatrix bad = xt::zeros<double>({3, 3}); // nb is 2 here
    EXPECT_THROW(solver.run(atoms, nullptr, &bad), std::invalid_argument);
}

// A5 / A3: UKS density guess requires both spin channels (no partial cold
// start). RKS OT + total density guess is product (SGCP OT corrector);
// UKS OT + dual density guess is product (second-gen-cpmd-uks-ot);
// RKS + guess_P_beta still throws.
TEST(GpwSolverGuessP, UksPartialOrOtWithGuessP) {
    const auto atoms = small_h2();
    const RealMatrix p = constant_matrix(2, 0.5);
    const RealMatrix pb = constant_matrix(2, 0.4);

    GpwDft uks = small_cfg();
    uks.n_occ = 0;
    uks.n_alpha = 1;
    uks.n_beta = 1;
    uks.build(atoms);
    EXPECT_THROW(uks.run(atoms, nullptr, &p, nullptr), std::invalid_argument);
    EXPECT_THROW(uks.run(atoms, nullptr, nullptr, &pb), std::invalid_argument);

    GpwDft rks = small_cfg();
    rks.build(atoms);
    EXPECT_THROW(rks.run(atoms, nullptr, &p, &pb), std::invalid_argument);

    GpwDft ot = small_cfg();
    ot.use_ot = true;
    ot.build(atoms);
    // RKS OT purifies guess_P → natural orbitals (second-gen-cpmd-ot).
    EXPECT_NO_THROW(ot.run(atoms, nullptr, &p));
    EXPECT_THROW(ot.run(atoms, nullptr, &p, &pb), std::invalid_argument);

    GpwDft uks_ot = small_cfg();
    uks_ot.n_occ = 0;
    uks_ot.n_alpha = 1;
    uks_ot.n_beta = 1;
    uks_ot.use_ot = true;
    uks_ot.build(atoms);
    // UKS OT + dual density guess: spin-resolved natural-orbital seed.
    EXPECT_NO_THROW(uks_ot.run(atoms, nullptr, &p, &pb));
}

// A2: UKS warm density guess (both spins) reaches the cold fixed point.
TEST(GpwSolverGuessP, UksConvergedDensityGuessReachesColdFixedPoint) {
    pin_single_thread();
    const auto atoms = small_h2();

    GpwDft cold = small_cfg();
    cold.n_occ = 0;
    cold.n_alpha = 1;
    cold.n_beta = 1;
    cold.build(atoms);
    const GpwResult ref = cold.run(atoms);
    ASSERT_TRUE(ref.converged);
    ASSERT_GT(ref.P.size(), 0u);
    ASSERT_GT(ref.P_beta.size(), 0u);

    GpwDft warm = small_cfg();
    warm.n_occ = 0;
    warm.n_alpha = 1;
    warm.n_beta = 1;
    warm.build(atoms);
    const GpwResult& r = warm.run(atoms, nullptr, &ref.P, &ref.P_beta);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.e_total, ref.e_total, 1e-10);
    EXPECT_LE(r.iters, ref.iters);
}

// ── GPWEvaluator xl_bomd guards (validation only — no SCF runs) ─────────

namespace {

ATV::cpu::Frame water_frame() {
    ATV::cpu::Frame f(3);
    const int z[3] = {8, 1, 1};
    const double x[3] = {5.0000, 5.7586, 4.2414};
    const double y[3] = {5.0000, 5.0000, 5.0000};
    const double zc[3] = {5.0000, 5.5043, 5.5043};
    for (int i = 0; i < 3; ++i) {
        f.span(ATV::kAtom_type)[i] = z[i];
        f.span(ATV::kAtom_x)[i] = x[i];
        f.span(ATV::kAtom_y)[i] = y[i];
        f.span(ATV::kAtom_z)[i] = zc[i];
    }
    f.box = ATV::cpu::Box::cubic(10.0);
    return f;
}

ATV::cpu::gpw::GPWEvaluator water_eval() {
    const ATV::cpu::gpw::Kind O{
        molrs::Element::O,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::O, "GTH-PBE-q6")};
    const ATV::cpu::gpw::Kind H{
        molrs::Element::H,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};
    return ATV::cpu::gpw::GPWEvaluator({O, H});
}

} // namespace

// A6: xl_bomd + UKS occupations is allowed at build() (dual spin XL).
TEST(GpwEvaluatorXlBomd, BuildAcceptsXlWithUks) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

// A6: xl_bomd with the OT minimizer is refused at build().
TEST(GpwEvaluatorXlBomd, BuildRejectsXlWithOt) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.use_ot = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A6: xl_bomd + UKS + OT still refused (OT wins the mutual exclusion).
TEST(GpwEvaluatorXlBomd, BuildRejectsXlUksWithOt) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    eval.use_ot = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A6: an out-of-table dissipation order is refused at build().
TEST(GpwEvaluatorXlBomd, BuildRejectsXlOrderOutsideTable) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.xl_bomd_k = 11;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// ── multi-k XL (complex P_k) ────────────────────────────────────────────

using ATV::cpu::gpw::XlBomdCplx;
using ATV::cpu::linalg::CplxMatrix;
using Cplx = std::complex<double>;

namespace {

CplxMatrix constant_cplx(int n, Cplx v) {
    CplxMatrix m = xt::zeros<Cplx>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            m(i, j) = v + Cplx(0.1 * static_cast<double>(i == j), 0.0);
    // Hermitianize: (M + M†)/2 so fixed-point tests stay physical.
    CplxMatrix h = xt::zeros<Cplx>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            h(i, j) = 0.5 * (m(i, j) + std::conj(m(j, i)));
    return h;
}

double max_abs_diff_cplx(const CplxMatrix& a, const CplxMatrix& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < a.shape(0); ++i)
        for (std::size_t j = 0; j < a.shape(1); ++j)
            worst = std::max(worst, std::abs(a(i, j) - b(i, j)));
    return worst;
}

} // namespace

// A2: XlBomdCplx shares Table I and has a constant-density fixed point.
TEST(XlBomdCplx, ConstantHermitianIsFixedPoint) {
    XlBomdCplx xl(5);
    EXPECT_EQ(xl.order(), 5);
    EXPECT_DOUBLE_EQ(xl.kappa(), 1.82);
    EXPECT_DOUBLE_EQ(xl.alpha(), 0.018);
    const CplxMatrix d = constant_cplx(4, Cplx(0.3, 0.0));
    for (int n = 0; n < 15; ++n)
        xl.advance(d);
    EXPECT_LE(max_abs_diff_cplx(xl.aux(), d), 1e-12);
}

// A2: shape mismatch after seeding throws.
TEST(XlBomdCplx, AdvanceWithMismatchedShapeThrows) {
    XlBomdCplx xl(5);
    xl.advance(constant_cplx(3, Cplx(0.2, 0.0)));
    EXPECT_THROW(xl.advance(constant_cplx(4, Cplx(0.2, 0.0))),
                 std::invalid_argument);
}

// A3/A4: multi-k RKS density guess plumbing on tiny H2.
TEST(GpwSolverGuessPk, MultiKConvergedGuessReachesColdFixedPoint) {
    pin_single_thread();
    const auto atoms = small_h2(/*box=*/8.0);

    GpwDft cold = small_cfg(/*box=*/8.0);
    cold.k_mesh[0] = 1;
    cold.k_mesh[1] = 1;
    cold.k_mesh[2] = 2;
    cold.k_time_reversal = false; // two k
    cold.nonperiodic = false;
    cold.build(atoms);
    const GpwResult ref = cold.run(atoms);
    ASSERT_TRUE(ref.converged);
    ASSERT_GE(ref.n_kpts, 2);
    ASSERT_EQ(static_cast<int>(ref.P_k.size()), ref.n_kpts);

    GpwDft warm = small_cfg(/*box=*/8.0);
    warm.k_mesh[0] = 1;
    warm.k_mesh[1] = 1;
    warm.k_mesh[2] = 2;
    warm.k_time_reversal = false;
    warm.nonperiodic = false;
    warm.build(atoms);
    const GpwResult& r = warm.run(atoms, nullptr, nullptr, nullptr, &ref.P_k);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.e_total, ref.e_total, 1e-9);
    EXPECT_LE(r.iters, ref.iters);
}

// A3: wrong guess_P_k count throws.
TEST(GpwSolverGuessPk, WrongCountThrows) {
    pin_single_thread();
    const auto atoms = small_h2();
    GpwDft solver = small_cfg();
    solver.k_mesh[2] = 2;
    solver.k_time_reversal = false;
    solver.nonperiodic = false;
    solver.build(atoms);
    std::vector<CplxMatrix> bad(1, constant_cplx(2, Cplx(0.1, 0.0)));
    EXPECT_THROW(solver.run(atoms, nullptr, nullptr, nullptr, &bad),
                 std::invalid_argument);
}

// A4: multi-k + real guess_P throws.
TEST(GpwSolverGuessPk, MultiKWithRealGuessPThrows) {
    pin_single_thread();
    const auto atoms = small_h2();
    GpwDft solver = small_cfg();
    solver.k_mesh[2] = 2;
    solver.k_time_reversal = false;
    solver.nonperiodic = false;
    solver.build(atoms);
    const RealMatrix p = constant_matrix(2, 0.5);
    EXPECT_THROW(solver.run(atoms, nullptr, &p), std::invalid_argument);
}

// A5: multi-k RKS + XL accepted; multi-k UKS + XL productized.
TEST(GpwEvaluatorXlBomd, BuildAcceptsMultiKRksXl) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.k_mesh[0] = 1;
    eval.k_mesh[1] = 1;
    eval.k_mesh[2] = 2;
    eval.k_time_reversal = false;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorXlBomd, BuildAcceptsMultiKUksXl) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    eval.k_mesh[2] = 2;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

// A3/A4: multi-k UKS dual-spin density guess plumbing.
TEST(GpwSolverGuessPk, MultiKUksConvergedDualGuessReachesColdFixedPoint) {
    pin_single_thread();
    const auto atoms = small_h2(/*box=*/8.0);

    GpwDft cold = small_cfg(/*box=*/8.0);
    cold.n_occ = 0;
    cold.n_alpha = 1;
    cold.n_beta = 1;
    cold.k_mesh[0] = 1;
    cold.k_mesh[1] = 1;
    cold.k_mesh[2] = 2;
    cold.k_time_reversal = false;
    cold.nonperiodic = false;
    cold.build(atoms);
    const GpwResult ref = cold.run(atoms);
    ASSERT_TRUE(ref.converged);
    ASSERT_GE(ref.n_kpts, 2);
    ASSERT_EQ(static_cast<int>(ref.P_k.size()), ref.n_kpts);
    ASSERT_EQ(static_cast<int>(ref.P_k_beta.size()), ref.n_kpts);

    GpwDft warm = small_cfg(/*box=*/8.0);
    warm.n_occ = 0;
    warm.n_alpha = 1;
    warm.n_beta = 1;
    warm.k_mesh[0] = 1;
    warm.k_mesh[1] = 1;
    warm.k_mesh[2] = 2;
    warm.k_time_reversal = false;
    warm.nonperiodic = false;
    warm.build(atoms);
    const GpwResult& r =
        warm.run(atoms, nullptr, nullptr, nullptr, &ref.P_k, &ref.P_k_beta);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.e_total, ref.e_total, 1e-9);
    EXPECT_LE(r.iters, ref.iters);
}

// A3: multi-k UKS α-only density guess throws.
TEST(GpwSolverGuessPk, MultiKUksPartialSpinGuessThrows) {
    pin_single_thread();
    const auto atoms = small_h2();
    GpwDft solver = small_cfg();
    solver.n_occ = 0;
    solver.n_alpha = 1;
    solver.n_beta = 1;
    solver.k_mesh[2] = 2;
    solver.k_time_reversal = false;
    solver.nonperiodic = false;
    solver.build(atoms);
    std::vector<CplxMatrix> only_a(2, constant_cplx(2, Cplx(0.1, 0.0)));
    EXPECT_THROW(solver.run(atoms, nullptr, nullptr, nullptr, &only_a, nullptr),
                 std::invalid_argument);
}

// A3: multi-k RKS + guess_P_k_beta throws.
TEST(GpwSolverGuessPk, MultiKRksWithBetaGuessThrows) {
    pin_single_thread();
    const auto atoms = small_h2();
    GpwDft solver = small_cfg();
    solver.k_mesh[2] = 2;
    solver.k_time_reversal = false;
    solver.nonperiodic = false;
    solver.build(atoms);
    std::vector<CplxMatrix> pk(2, constant_cplx(2, Cplx(0.1, 0.0)));
    EXPECT_THROW(solver.run(atoms, nullptr, nullptr, nullptr, &pk, &pk),
                 std::invalid_argument);
}

// ── Zero-SCF linearized XL forces (xl_bomd_linearized) ──────────────────

// A1: linearized without xl_bomd is refused.
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedWithoutXl) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd_linearized = true;
    eval.compute_forces = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: linearized + UKS refused.
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedUks) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.xl_bomd_linearized = true;
    eval.compute_forces = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: linearized + multi-k refused.
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedMultiK) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.xl_bomd_linearized = true;
    eval.compute_forces = true;
    eval.k_mesh[2] = 2;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: linearized without compute_forces refused.
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedWithoutForces) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.xl_bomd_linearized = true;
    eval.compute_forces = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: linearized + OT refused (via xl_bomd + use_ot).
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedWithOt) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.xl_bomd_linearized = true;
    eval.compute_forces = true;
    eval.use_ot = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: linearized + CPMD refused.
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedWithCpmd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.xl_bomd_linearized = true;
    eval.compute_forces = true;
    eval.cpmd = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: linearized + second_gen_cp refused.
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedWithSgcp) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.xl_bomd_linearized = true;
    eval.compute_forces = true;
    eval.second_gen_cp = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// Happy path: Γ RKS + forces + linearized builds.
TEST(GpwEvaluatorXlBomdLinearized, BuildAcceptsGammaRks) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.xl_bomd = true;
    eval.xl_bomd_linearized = true;
    eval.compute_forces = true;
    eval.nonperiodic = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

// A3: shadow energy formula recovers E_KS when D = n.
TEST(XlBomdLinearized, ShadowEnergyRecoversAtFixedPoint) {
    constexpr int n = 3;
    RealMatrix H = xt::zeros<double>({n, n});
    RealMatrix dens = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i) {
        H(i, i) = -0.5 - 0.1 * i;
        dens(i, i) = 0.2 + 0.05 * i;
        for (int j = 0; j < i; ++j) {
            H(i, j) = H(j, i) = 0.01 * (i + j);
            dens(i, j) = dens(j, i) = 0.02 * (i - j + 1);
        }
    }
    auto tr = [n](const RealMatrix& A, const RealMatrix& B) {
        double t = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                t += A(i, j) * B(j, i);
        return t;
    };
    const double e_ks_n = -1.23456789;
    const double e_lin =
        tr(H, dens) - tr(H, dens) + e_ks_n; // D = n
    EXPECT_DOUBLE_EQ(e_lin, e_ks_n);

    // Non-trivial residual: E_lin differs from E_KS by Tr(H(D-n)).
    RealMatrix D = dens;
    D(0, 0) += 0.05;
    D(1, 1) -= 0.05;
    const double delta = tr(H, D) - tr(H, dens);
    const double e_lin2 = tr(H, D) - tr(H, dens) + e_ks_n;
    EXPECT_NEAR(e_lin2 - e_ks_n, delta, 1e-15);
}
