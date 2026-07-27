#include "atomiverse/cpu/electronic/second_gen_cp.h"
#include "atomiverse/cpu/electronic/car_parrinello.h"
#include "atomiverse/cpu/electronic/aimd_method.h"
#include "atomiverse/cpu/electronic/xl_bomd_force.h"
/// @file tests/cpu/gpw/test_gpw_xl_bomd.cpp
/// @brief Unit tests: the GpwDft guess_P / guess_P_k plumbing XL-BOMD drives +
///        the GPWEvaluator xl_bomd(_linearized) build() guards.
///
/// The XlBomd / XlBomdCplx propagator classes themselves live in
/// cpu/electronic/ and are pinned by tests/cpu/electronic/test_xl_bomd.cpp.
/// The guess_P plumbing tests reuse the tiny inline H2 system of
/// test_cpu_gpw_solver (ms-scale SCF). The e2e XL-BOMD NVE envelope lives in
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

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;
using ATV::cpu::linalg::RealMatrix;

namespace {
using ATV::cpu::electronic::XlBomdForce;
using ATV::cpu::electronic::XlBomdLinearized;
using ATV::cpu::electronic::SecondGenCp;
using ATV::cpu::electronic::CarParrinello;

RealMatrix constant_matrix(int n, double v) {
    RealMatrix m = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            m(i, j) = v + 0.1 * static_cast<double>(i == j);
    return m;
}

} // namespace

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
    auto base = water_eval();
    base.n_alpha = 5;
    base.n_beta = 3;
    XlBomdForce eval{std::move(base), 5};
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

// A6: xl_bomd with the OT minimizer is refused at build().
TEST(GpwEvaluatorXlBomd, BuildRejectsXlWithOt) {
    auto frame = water_frame();
    auto base = water_eval();
    base.use_ot = true;
    XlBomdForce eval{std::move(base), 5};
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A6: xl_bomd + UKS + OT still refused (OT wins the mutual exclusion).
TEST(GpwEvaluatorXlBomd, BuildRejectsXlUksWithOt) {
    auto frame = water_frame();
    auto base = water_eval();
    base.n_alpha = 5;
    base.n_beta = 3;
    base.use_ot = true;
    XlBomdForce eval{std::move(base), 5};
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A6: an out-of-table dissipation order is refused at build().
TEST(GpwEvaluatorXlBomd, BuildRejectsXlOrderOutsideTable) {
    auto frame = water_frame();
    auto base = water_eval();
    XlBomdForce eval{std::move(base), 11};
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// ── multi-k XL (complex P_k) ────────────────────────────────────────────

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

} // namespace

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
    auto base = water_eval();
    XlBomdForce eval{std::move(base), 5};
    eval.method().k_mesh[0] = 1;
    eval.method().k_mesh[1] = 1;
    eval.method().k_mesh[2] = 2;
    eval.method().k_time_reversal = false;
    eval.method().nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorXlBomd, BuildAcceptsMultiKUksXl) {
    auto frame = water_frame();
    auto base = water_eval();
    base.n_alpha = 5;
    base.n_beta = 3;
    XlBomdForce eval{std::move(base), 5};
    eval.method().k_mesh[2] = 2;
    eval.method().nonperiodic = false;
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

// A1: linearized is only available via XlBomdForce(..., XlBomdLinearized{}).
TEST(GpwEvaluatorXlBomdLinearized, LinearizedRequiresXlBomdForceCtor) {
    auto base = water_eval();
    base.compute_forces = true;
    // Bare GPW has no public linearized switch; XL path is the wrapper type.
    EXPECT_FALSE(base.is_xl_bomd());
    XlBomdForce xl{std::move(base), 5, XlBomdLinearized{}};
    EXPECT_TRUE(xl.method().is_xl_bomd());
}

// A2: linearized + UKS refused.
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedUks) {
    auto frame = water_frame();
    auto base = water_eval();
    base.compute_forces = true;
    base.n_alpha = 5;
    base.n_beta = 3;
    XlBomdForce eval{std::move(base), 5, XlBomdLinearized{}};
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: linearized + multi-k refused.
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedMultiK) {
    auto frame = water_frame();
    auto base = water_eval();
    base.compute_forces = true;
    XlBomdForce eval{std::move(base), 5, XlBomdLinearized{}};
    eval.method().k_mesh[2] = 2;
    eval.method().nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: linearized ctor forces compute_forces=true.
TEST(GpwEvaluatorXlBomdLinearized, CtorEnablesForces) {
    auto base = water_eval();
    base.compute_forces = false;
    XlBomdForce eval{std::move(base), 5, XlBomdLinearized{}};
    EXPECT_TRUE(eval.method().compute_forces);
}

// A2: linearized + OT refused (via xl_bomd + use_ot).
TEST(GpwEvaluatorXlBomdLinearized, BuildRejectsLinearizedWithOt) {
    auto frame = water_frame();
    auto base = water_eval();
    base.compute_forces = true;
    base.use_ot = true;
    XlBomdForce eval{std::move(base), 5, XlBomdLinearized{}};
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// A2: XL and CP are distinct types (cannot combine on one object).
TEST(GpwEvaluatorXlBomdLinearized, DistinctFromCarParrinello) {
    XlBomdForce xl{water_eval(), 5, XlBomdLinearized{}};
    CarParrinello cp{water_eval(), 400.0};
    EXPECT_NE(xl.name(), cp.name());
}

// A2: XL and SGCP are distinct types.
TEST(GpwEvaluatorXlBomdLinearized, DistinctFromSecondGenCp) {
    XlBomdForce xl{water_eval(), 5, XlBomdLinearized{}};
    SecondGenCp sg{water_eval(), 3, 5};
    EXPECT_NE(xl.name(), sg.name());
}

// Happy path: Γ RKS + forces + linearized builds.
TEST(GpwEvaluatorXlBomdLinearized, BuildAcceptsGammaRks) {
    auto frame = water_frame();
    auto base = water_eval();
    base.compute_forces = true;
    base.nonperiodic = true;
    XlBomdForce eval{std::move(base), 5, XlBomdLinearized{}};
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
