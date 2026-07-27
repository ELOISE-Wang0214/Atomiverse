/// @file tests/cpu/scf/test_scf_mc_predictor_guard.cpp
/// @brief esmm-contract-04-forces-requested (U6) — the ASPC-predictor guard:
///        an extrapolating density guess is MD-only and must refuse to run
///        under Monte Carlo sampling.
///
/// Mirrors one production unit: `src/cpu/scf/mean_field_method.cpp`
/// (`MeanFieldMethod<D>::run_nested_erased`, the nested entry every outer
/// driver reaches a Gaussian method through). One seam per file, as the rest of
/// this directory does it (`test_scf_cpmd.cpp` for the Car-Parrinello refusal
/// lattice, `test_scf_xl_bomd.cpp` for the XL-BOMD guess ladder).
///
/// ## Why a throw is the correct minimal thing
///
/// MC's `mc::Transaction::rollback` restores the geometry of a rejected trial
/// bitwise, but it cannot notify an evaluator — so the rejected move's density
/// stays in the method's warm-start state. Two warm-start flavours, two very
/// different consequences:
///
///   - **plain last-density reuse** is correctness-neutral: it is a *guess*,
///     and the SCF converges to the geometry's own fixed point to `tol` no
///     matter how it was seeded. Nothing to forbid;
///   - the **ASPC predictor** extrapolates a *time-continuous* density history
///     and then runs only `corrector_steps_` corrector passes at the loosened
///     `corrector_tol_`. Accept/reject jumps violate the continuity the
///     extrapolation assumes, so the returned energy becomes a function of the
///     rejected history — the Metropolis test is then applied to
///     history-dependent energies and detailed balance is broken silently.
///
/// A rollback-notification hook was considered and rejected (cross-cutting
/// machinery for one consumer, and the loosened corrector would still bias the
/// energies even with notification). What is pinned here is the alternative: a
/// **typed, actionable refusal** naming `predictor()` and the sampling context.
///
/// ## Cases
///
///   - `PredictorUnderSamplingThrows` (ac-003, RED today: the sweep completes
///     silently) — `run_nested_erased` currently opens with `(void)ctx;`
///     (`src/cpu/scf/mean_field_method.cpp:346`), so it cannot see that it is
///     inside a trial at all;
///   - `LastDensityWarmStartUnderSamplingIsFine` — the negative leg: the
///     guard must key on the predictor, not on "nested under MC" (control,
///     passes today);
///   - `PredictorUnderMdStillRuns` — the guard must not cost MD its ASPC warm
///     start, which is the whole reason the predictor exists (control).
///
/// The Context channel the guard reads (`ctx.mc.proposal != nullptr`) is
/// published by `McBase::run` before the first nested drive, so the real MC
/// scaffold is both the cheapest and the only honest harness: a hand-built
/// Context would pin the guard's *condition* rather than the situation it
/// exists for.
///
/// **Scope.** Unit only. H₂/STO-3G (2 AOs) at a loose tolerance, one sweep;
/// the sampling physics (a QM-in-MC bond-length distribution pinned against
/// the analytic stiff-spring variance) is `regressions/mc/mc_bias_es`.
///
/// **Determinism.** Fixed geometry, fixed MC seed, no wall clock, no
/// filesystem, no network, no third-party oracle. No energy value is asserted
/// — only refusal / non-refusal.
///
/// Spec: `.claude/specs/esmm-contract-04-forces-requested/spec.md` (ac-003).

#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/mover.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/velocity_verlet.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/montecarlo/mc.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

using cpu::scf::BasisSet;

namespace {

/// Loose SCF convergence [Hartree] — this file asserts refusal, never a value.
constexpr double kLooseTol = 1e-6;

/// ASPC order used by every predictor leg below (the `predictor()` default).
constexpr int kAspcOrder = 3;

/// Nuclear timestep [fs] for the MD control leg.
constexpr double kDtFs = 0.5;

bool mentions(const std::exception& e, const std::string& needle) {
    return std::string(e.what()).find(needle) != std::string::npos;
}

/// @brief H₂ at r = 0.74 Å in a 20 Å vacuum box (Z = 1, mass 1.008 amu).
cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_mass)[1] = 1.008;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

/// @brief RHF/STO-3G, loose tolerance, no gradient (nothing here reads forces).
cpu::scf::RHF h2_rhf() {
    cpu::scf::RHF rhf;
    rhf.basis = BasisSet::STO3G;
    rhf.tol = kLooseTol;
    rhf.max_iter = 50;
    rhf.compute_force = false;
    return rhf;
}

/// @brief Minimal `Mover`-kind stub — MC requires one and this file does not
///        care what a trial move actually does.
struct StubMover final : cpu::Mover {
    std::string name() const override { return "stub_mover"; }
    void launch(const cpu::LaunchContext&) override {}
};

/// @brief Bind `method` as MC's nested energy source and take one sweep.
///
/// Validation may fire at the nested `build_on` or at the first nested drive,
/// and both sites are inside this call — so a refusal assertion does not depend
/// on which one the implementation picks (the same helper shape
/// `tests/cpu/scf/test_scf_xl_bomd.cpp` uses for its refusal rows).
/// @param method Borrowed concrete method (the driver does not own it).
/// @param frame  Geometry driven for one sweep.
void drive_nested_mc_sweep(cpu::scf::RHF& method, cpu::Frame& frame) {
    cpu::mc::MC mc{/*temperature_K=*/300.0, /*seed=*/12345u};
    mc.on(sampling::Phase::propose).add(StubMover{}, /*order=*/0);
    mc.on(sampling::Phase::energy).add_ref(method, /*order=*/0);
    mc.build();
    mc.run(frame, /*n_sweeps=*/1);
}

/// @brief Bind `method` at the MD force phase and take one step.
void drive_nested_md_step(cpu::scf::RHF& method, cpu::Frame& frame) {
    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{kDtFs});
    sim.on(md::Phase::force).add_ref(method);
    sim.build();
    sim.run(frame, /*n_steps=*/1);
}

} // namespace

// ── ac-003: the extrapolating guess refuses to sample ───────────────────────

/// RED today: the sweep runs to completion and the ASPC history quietly
/// accumulates whatever the rejected trials proposed.
TEST(MeanFieldMethodSamplingGuard, PredictorUnderSamplingThrows) {
    cpu::Frame frame = make_h2();
    cpu::scf::RHF rhf = h2_rhf();
    rhf.predictor(kAspcOrder);

    try {
        drive_nested_mc_sweep(rhf, frame);
        FAIL() << "a Monte Carlo sweep drove an ASPC-predictor SCF without "
                  "complaint — rejected moves stay in the density history and "
                  "the loosened corrector makes the trial energies "
                  "history-dependent, so detailed balance is broken with no "
                  "diagnostic at all";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(mentions(e, "predictor"))
            << "the refusal must name `predictor()` — the one knob the user "
               "has to turn off. Got: "
            << e.what();
        EXPECT_TRUE(mentions(e, "sampling"))
            << "the refusal must name the sampling context, so the reader "
               "knows it is MC and not the method that is unsupported. Got: "
            << e.what();
        EXPECT_EQ(rhf.iterations(), 0)
            << "the guard fired only after an SCF had already been spent on a "
               "configuration it was going to refuse";
    }
}

/// The negative leg of ac-003: plain last-density reuse is a *guess*, and a
/// converged SCF does not remember it. Nesting a mean-field method under MC
/// must stay legal — that is the whole point of the wave. (Passes today; it
/// pins the GREEN against a blanket "no Gaussian method under MC" refusal.)
TEST(MeanFieldMethodSamplingGuard, LastDensityWarmStartUnderSamplingIsFine) {
    cpu::Frame frame = make_h2();
    cpu::scf::RHF rhf = h2_rhf(); // no predictor(): last-P reuse only

    EXPECT_NO_THROW(drive_nested_mc_sweep(rhf, frame));
    EXPECT_GT(rhf.iterations(), 0)
        << "the nested method never ran, so this leg proves nothing about the "
           "guard";
    EXPECT_TRUE(rhf.converged());
}

/// The guard is scoped to sampling: MD is exactly where the ASPC predictor is
/// correct and valuable (a time-continuous nuclear trajectory), so the same
/// configuration must still run there. (Control, passes today.)
TEST(MeanFieldMethodSamplingGuard, PredictorUnderMdStillRuns) {
    cpu::Frame frame = make_h2();
    cpu::scf::RHF rhf = h2_rhf();
    rhf.predictor(kAspcOrder);

    EXPECT_NO_THROW(drive_nested_md_step(rhf, frame));
    EXPECT_GT(rhf.iterations(), 0);
    EXPECT_TRUE(rhf.converged());
}
