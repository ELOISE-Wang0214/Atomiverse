/// @file tests/cpu/scf/test_method_driver_parity.cpp
/// @brief scf-driver-refactor-03a — RHF/UHF/RKS/UKS become
///        `cpu::CPUDriver<ATV::scf::Phase>` and stay bit-parity with the
///        RETAINED in-tree `SCFSolver` reference.
///
/// ## RED-first (this file fails to COMPILE until 03a lands)
///
/// The four `static_assert`s below require `MeanFieldMethod<D>` to derive
/// `cpu::CPUDriver<ATV::scf::Phase>`. Today it derives only `HF`, so each
/// assertion is `false` and the translation unit does not compile — that is the
/// intended RED (spec Testing strategy: "断言 `MeanFieldMethod` 派生
/// `CPUDriver<scf::Phase>`"). The wiring test additionally calls the loop-node
/// API `converged()` / `stop_reason()`, which the concrete methods only gain by
/// deriving `CPUDriver`; those missing members are the same not-yet-implemented
/// contract, not harness breakage.
///
/// ## What the GREEN suite pins
///
///   - **Parity (the domain validation).** For each of RHF/UHF/RKS/UKS on the
///     same-scale small molecules the existing `tests/cpu/scf/*` use
///     (H₂ / H₂O / He / H, STO-3G), the NEW driver path (`method.run(frame)`)
///     matches the in-tree `SCFSolver` reference — built on an IDENTICAL
///     `IntegralPack` + `SCFConfig` mirroring the method's private
///     `make_config` — to `1e-12` Ha, plus eig / P / spin-contamination
///     spot-checks (acceptance ac-001..004). This is a unit-tier check: same
///     scale as `test_scf_solver`, seconds, no out-of-library consumer, no
///     literature value (that e2e is 03b's `regressions/scf/<case>`).
///   - **One funnel, two entries.** The same molecule via `run(frame)` and via
///     binding the method as an `HF` component into `cpu::Evaluate` (the
///     `launch` path) yields the identical SCF energy — proving both route
///     through the private `drive_scf_` (ac-006).
///   - **Wiring.** After a converging run the loop node reports `converged()`
///     and `stop_reason() == "scf_converged"` — evidence that the ctor
///     registered `ScfConverged` via `add_criterion` and the driver polled it
///     (ac-005).
///   - **Edges.** `level_shift ≠ 0` converges and reports de-shifted orbital
///     energies matching the reference; `mp2()` / `fci()` still route through
///     the driven method unchanged (ac-007).
///
/// The reference side deliberately uses `SCFConfig` / `SCFSolver` /
/// `IntegralPack` — the RETAINED living-parity machinery (spec 03a). The NEW
/// path (`method.run`, `Evaluate` + the method component) never names them.
/// Threads are pinned to 1 in CMake so the `1e-12` bit-parity is deterministic
/// (see `test_cpu_scf_stages` for the same reason).

#ifdef ATV_HAS_RHF

#include <cmath>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>

#include <gtest/gtest.h>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h" // RHF/UHF/RKS/UKS umbrella (NEW path)
// Retained living-parity reference only — SCFSolver + SCFConfig + IntegralPack.
#include "atomiverse/cpu/scf/internal/scf_solver.h"
#include "atomiverse/facade/stage.h" // ATV::scf::Phase, ATV::evaluate::Phase

using namespace ATV;

// ── Compile-time RED: each concrete method IS a CPUDriver<scf::Phase> ──────────
//
// TODAY these are `false` (RHF/UHF/RKS/UKS derive `HF`, not the loop node), so
// the file fails to compile. After 03a `MeanFieldMethod<D>` derives BOTH
// `CPUDriver<ATV::scf::Phase>` and `HF`, and all four flip to `true`.
static_assert(std::is_base_of_v<cpu::CPUDriver<ATV::scf::Phase>, cpu::scf::RHF>,
              "03a: RHF must derive cpu::CPUDriver<ATV::scf::Phase>");
static_assert(std::is_base_of_v<cpu::CPUDriver<ATV::scf::Phase>, cpu::scf::UHF>,
              "03a: UHF must derive cpu::CPUDriver<ATV::scf::Phase>");
static_assert(std::is_base_of_v<cpu::CPUDriver<ATV::scf::Phase>, cpu::scf::RKS>,
              "03a: RKS must derive cpu::CPUDriver<ATV::scf::Phase>");
static_assert(std::is_base_of_v<cpu::CPUDriver<ATV::scf::Phase>, cpu::scf::UKS>,
              "03a: UKS must derive cpu::CPUDriver<ATV::scf::Phase>");

namespace {

// ── Geometries reused verbatim from the existing SCF solver tests ─────────────

/// @brief H₂O / STO-3G (10 electrons) — test_scf_solver.cpp:17-36 and
///        test_uks_solver.cpp:27-46 (identical).
cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_x)[1] = 0.7572;
    f.span(kAtom_x)[2] = -0.7572;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_y)[2] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_z)[1] = -0.5860;
    f.span(kAtom_z)[2] = -0.5860;
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_mass)[0] = 15.999;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_mass)[2] = 1.008;
    f.simbox = cpu::Box::cubic(50.0);
    return f;
}

/// @brief Single H atom (1 electron) — test_uhf_solver.cpp:16-25.
cpu::Frame make_h() {
    cpu::Frame f(1);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.simbox = cpu::Box::cubic(20.0);
    return f;
}

/// @brief He atom (2 electrons) — test_uks_solver.cpp:16-25 / test_mp2.cpp:16-25.
cpu::Frame make_he() {
    cpu::Frame f(1);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[0] = 2;
    f.span(kAtom_mass)[0] = 4.0026;
    f.simbox = cpu::Box::cubic(20.0);
    return f;
}

/// @brief H₂ at r = 0.74 Å (2 electrons) — test_mp2.cpp:27-41 /
///        test_level_shift.cpp:6-20.
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
    f.simbox = cpu::Box::cubic(20.0);
    return f;
}

// ── FP tolerances (project scientific-test table) ─────────────────────────────
constexpr double kEnergyExact = 1e-12;  ///< Ha — spec 03a parity mandate.
constexpr double kEigExact = 1e-10;     ///< Ha — orbital-energy spot-check.
constexpr double kDensityExact = 1e-10; ///< density-matrix / spin spot-check.

void expect_vector_near(const xt::xtensor<double, 1>& got,
                        const xt::xtensor<double, 1>& ref, double tol) {
    ASSERT_EQ(got.shape(0), ref.shape(0));
    for (std::size_t i = 0; i < got.shape(0); ++i)
        EXPECT_NEAR(got(i), ref(i), tol) << "index " << i;
}

void expect_matrix_near(const xt::xtensor<double, 2>& got,
                        const xt::xtensor<double, 2>& ref, double tol) {
    ASSERT_EQ(got.shape(0), ref.shape(0));
    ASSERT_EQ(got.shape(1), ref.shape(1));
    for (std::size_t i = 0; i < got.shape(0); ++i)
        for (std::size_t j = 0; j < got.shape(1); ++j)
            EXPECT_NEAR(got(i, j), ref(i, j), tol)
                << "(" << i << "," << j << ")";
}

// ── Living-parity reference (SCFSolver on an identical pack + config) ──────────

/// @brief The `SCFConfig` the method's private `make_config` produces: every
///        shared knob copied from the method's PUBLIC members, and only the
///        concrete (spin, mean-field type, xc, occupation) tuple supplied per
///        method. Copying from members keeps the reference tracking any knob the
///        test overrides (e.g. `level_shift`, the RKS grid) so parity is exact.
template <typename Method>
cpu::scf::SCFConfig reference_config(const Method& m, cpu::scf::SpinType spin,
                                     cpu::scf::MeanFieldType mf,
                                     cpu::scf::XCType xc, int n_alpha,
                                     int n_beta) {
    cpu::scf::SCFConfig c;
    c.spin = spin;
    c.method = mf;
    c.xc_type = xc;
    c.n_alpha = n_alpha;
    c.n_beta = n_beta;
    c.tol = m.tol;
    c.conv_tol_grad = m.conv_tol_grad;
    c.max_iter = m.max_iter;
    c.use_diis = m.use_diis;
    c.diis_max_vectors = m.diis_max_vectors;
    c.diis_start_iter = m.diis_start_iter;
    c.damping_factor = m.damping_factor;
    c.init_guess = m.init_guess;
    c.level_shift = m.level_shift;
    c.lindep_threshold = m.lindep_threshold;
    c.grid_level = m.grid_level;
    c.n_radial = m.n_radial;
    c.n_angular = m.n_angular;
    c.force_native_two_electron_deriv = m.force_native_two_electron_deriv;
    c.xc_grad_project_translation = m.xc_grad_project_translation;
    c.hybrid_alpha_override = m.hybrid_alpha_override;
    c.fixed_grid = m.fixed_grid;
    return c;
}

/// @brief Run the retained `SCFSolver` reference on the same geometry + basis as
///        `m`, returning a copy of its converged `SCFResult` [Hartree].
template <typename Method>
cpu::scf::SCFResult
reference_result(const Method& m, const cpu::Frame& frame,
                 cpu::scf::SpinType spin, cpu::scf::MeanFieldType mf,
                 cpu::scf::XCType xc, int n_alpha, int n_beta) {
    cpu::scf::IntegralPack ints{frame, m.basis, m.charge, m.compute_eri};
    cpu::scf::SCFConfig cfg =
        reference_config(m, spin, mf, xc, n_alpha, n_beta);
    cpu::scf::SCFSolver ref{cfg};
    ref.build(ints);
    return ref.run(ints);
}

} // namespace

// ── Per-method parity (ac-001..004) ───────────────────────────────────────────

TEST(MethodDriverParity, RhfWaterMatchesReference) {
    auto frame = make_h2o();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    const cpu::scf::SCFResult& got = rhf.run(frame);

    // H₂O closed-shell: 10 electrons → 5 doubly-occupied.
    const cpu::scf::SCFResult ref =
        reference_result(rhf, frame, cpu::scf::SpinType::Restricted,
                         cpu::scf::MeanFieldType::HF, cpu::scf::XCType::LSDA,
                         /*n_alpha=*/5, /*n_beta=*/5);

    ASSERT_TRUE(got.converged);
    ASSERT_TRUE(ref.converged);
    EXPECT_NEAR(got.energy, ref.energy, kEnergyExact);
    expect_vector_near(got.eigenvalues_alpha, ref.eigenvalues_alpha, kEigExact);
    expect_matrix_near(got.P_alpha, ref.P_alpha, kDensityExact);
}

TEST(MethodDriverParity, UhfHydrogenMatchesReference) {
    auto frame = make_h();
    cpu::scf::UHF uhf(/*n_alpha=*/1, /*n_beta=*/0);
    uhf.basis = cpu::scf::BasisSet::STO3G;
    const cpu::scf::SCFResult& got = uhf.run(frame);

    const cpu::scf::SCFResult ref =
        reference_result(uhf, frame, cpu::scf::SpinType::Unrestricted,
                         cpu::scf::MeanFieldType::HF, cpu::scf::XCType::LSDA,
                         uhf.n_alpha, uhf.n_beta);

    ASSERT_TRUE(got.converged);
    ASSERT_TRUE(ref.converged);
    EXPECT_NEAR(got.energy, ref.energy, kEnergyExact);
    EXPECT_NEAR(got.spin_contamination, ref.spin_contamination, kDensityExact);
    expect_matrix_near(got.P_alpha, ref.P_alpha, kDensityExact);
}

TEST(MethodDriverParity, RksWaterMatchesReference) {
    auto frame = make_h2o();
    cpu::scf::RKS rks(cpu::scf::XCType::PBE);
    rks.basis = cpu::scf::BasisSet::STO3G;
    // Light grid for unit-test speed (test_uks_solver.cpp:69-80); the reference
    // tracks it because reference_config copies grid_level/n_radial/n_angular.
    rks.grid_level = cpu::scf::GridLevel::SG0;
    rks.n_radial = 30;
    rks.n_angular = 50;
    const cpu::scf::SCFResult& got = rks.run(frame);

    const cpu::scf::SCFResult ref = reference_result(
        rks, frame, cpu::scf::SpinType::Restricted, cpu::scf::MeanFieldType::KS,
        rks.xc, /*n_alpha=*/5, /*n_beta=*/5);

    ASSERT_TRUE(got.converged);
    ASSERT_TRUE(ref.converged);
    EXPECT_NEAR(got.energy, ref.energy, kEnergyExact);
    expect_vector_near(got.eigenvalues_alpha, ref.eigenvalues_alpha, kEigExact);
}

TEST(MethodDriverParity, UksHeliumMatchesReference) {
    auto frame = make_he();
    cpu::scf::UKS uks(/*n_alpha=*/1, /*n_beta=*/1, cpu::scf::XCType::LSDA);
    uks.basis = cpu::scf::BasisSet::STO3G;
    const cpu::scf::SCFResult& got = uks.run(frame);

    const cpu::scf::SCFResult ref = reference_result(
        uks, frame, cpu::scf::SpinType::Unrestricted,
        cpu::scf::MeanFieldType::KS, uks.xc, uks.n_alpha, uks.n_beta);

    ASSERT_TRUE(got.converged);
    ASSERT_TRUE(ref.converged);
    EXPECT_NEAR(got.energy, ref.energy, kEnergyExact);
    expect_matrix_near(got.P_alpha, ref.P_alpha, kDensityExact);
    expect_matrix_near(got.P_beta, ref.P_beta, kDensityExact);
}

// ── Wiring: the loop node knows convergence via the ScfConverged criterion ─────
// (ac-005) — probed behaviorally through the public CPUDriver API, not by
// introspecting the Plan.

TEST(MethodDriverParity, ConvergenceWiringFlipsDriverStatus) {
    auto frame = make_h2o();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    const cpu::scf::SCFResult& got = rhf.run(frame);

    EXPECT_TRUE(got.converged);   // SCFResult flag.
    EXPECT_TRUE(rhf.converged()); // the LOOP knows (CPUDriver::converged).
    EXPECT_EQ(rhf.stop_reason(), "scf_converged"); // ScfConverged ended it.
}

// ── One funnel, two entries: run(frame) vs Evaluate + method component ─────────
// (ac-006)

TEST(MethodDriverParity, RunAndComponentLaunchAgree) {
    auto frame_run = make_h2o();
    cpu::scf::RHF direct;
    direct.basis = cpu::scf::BasisSet::STO3G;
    const double e_run = direct.run(frame_run).energy; // [Hartree]

    // Same molecule, but driven through the Component `launch` path by binding
    // the method into Evaluate. compute_force=false so launch writes the pure
    // SCF energy (no gradient recompute) and last_result() carries it.
    auto frame_eval = make_h2o();
    cpu::scf::RHF component;
    component.basis = cpu::scf::BasisSet::STO3G;
    component.with_gradient(false);

    cpu::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate).add_ref(component, /*order=*/0);
    eval.build();
    eval.run(frame_eval);

    // Both entries route through the same private drive_scf_ → identical energy.
    EXPECT_NEAR(component.last_result().energy, e_run, kEigExact);
    // Frame is the currency: the component wrote a bound electronic energy [eV].
    const double e_meta = static_cast<double>(frame_eval.meta.get(kEnergy));
    EXPECT_TRUE(std::isfinite(e_meta));
    EXPECT_LT(e_meta, 0.0);
    EXPECT_NEAR(e_meta, e_run * units::kHartreeInEv, 1e-6);
}

// ── Edge: level shift converges and reports de-shifted orbital energies ────────

TEST(MethodDriverParity, LevelShiftMatchesReferenceDeShifted) {
    auto frame = make_h2(); // test_level_shift.cpp geometry (H₂, r = 0.74).
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.level_shift = 0.5;
    const cpu::scf::SCFResult& got = rhf.run(frame);

    const cpu::scf::SCFResult ref =
        reference_result(rhf, frame, cpu::scf::SpinType::Restricted,
                         cpu::scf::MeanFieldType::HF, cpu::scf::XCType::LSDA,
                         /*n_alpha=*/1, /*n_beta=*/1);

    ASSERT_TRUE(got.converged);
    ASSERT_TRUE(ref.converged);
    EXPECT_NEAR(got.energy, ref.energy, kEnergyExact);
    // De-shift post-processing must reproduce the reference's unshifted eigs.
    expect_vector_near(got.eigenvalues_alpha, ref.eigenvalues_alpha, kEigExact);
}

// ── Edge: mp2()/fci() still route through the driven method (ac-007) ───────────

TEST(MethodDriverParity, Mp2AndFciStillRouteThroughDrivenMethod) {
    // He / STO-3G has no virtuals → both correlation corrections vanish, and
    // both must sit on top of the driven-method SCF energy (test_mp2.cpp,
    // test_fci.cpp inputs).
    auto frame_mp2 = make_he();
    cpu::scf::RHF rhf_mp2;
    rhf_mp2.basis = cpu::scf::BasisSet::STO3G;
    rhf_mp2.tol = 1e-12;
    const cpu::scf::MP2Result mp2 = rhf_mp2.mp2(frame_mp2);
    EXPECT_TRUE(rhf_mp2.last_result().converged);
    EXPECT_NEAR(mp2.e_corr, 0.0, 1e-14);
    EXPECT_NEAR(mp2.e_tot, rhf_mp2.last_result().energy, 1e-12);

    auto frame_fci = make_he();
    cpu::scf::RHF rhf_fci;
    rhf_fci.basis = cpu::scf::BasisSet::STO3G;
    rhf_fci.tol = 1e-12;
    const cpu::scf::FCIResult fci = rhf_fci.fci(frame_fci);
    EXPECT_TRUE(rhf_fci.last_result().converged);
    EXPECT_NEAR(fci.e_corr, 0.0, 1e-12);
    EXPECT_NEAR(fci.e_fci, rhf_fci.last_result().energy, 1e-12);
}

#endif // ATV_HAS_RHF
