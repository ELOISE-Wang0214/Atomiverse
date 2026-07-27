/// @file tests/cpu/scf/test_scf_aspc_fallback_accumulate.cpp
/// @brief Unit: the ONE entry that runs two SCF epilogues — the ASPC
///        corrector-miss fallback — must leave exactly ONE contribution in the
///        driver-owned accumulator channels.
///
/// Mirrors one production unit: `include/atomiverse/cpu/scf/mean_field_method.h`
/// (`MeanFieldMethod<Derived>` — the `predictor(order, corrector_steps,
/// corrector_tol, fallback)` surface and the fallback branch of `drive_scf_`,
/// `src/cpu/scf/mean_field_method.cpp:558-626`). Suite
/// `MeanFieldMethodAspcFallback`, sibling to `MeanFieldMethodXlBomd`
/// (test_scf_xl_bomd.cpp) and `MeanFieldMethodMcPredictorGuard`
/// (test_scf_mc_predictor_guard.cpp) over the same unit.
///
/// ## The path, and why it is the only one of its kind
///
/// `drive_scf_` normally calls `drive_scf_core_` exactly once, so its epilogue
/// (`:1234-1283`) writes the accumulators once. On **one** path it calls it
/// twice in a single entry: when the ASPC predictor produced a guess, the
/// short corrector solve (`max_iter = corrector_steps`,
/// `tol = corrector_tol`) failed to meet that tolerance, and
/// `fallback_full_scf_` is set — then the full solve re-runs
/// (`:600` then `:620`).
///
/// Since `esmm-contract-05` the **in-plan** entry (`run_nested_erased`) writes
/// `kEnergy` / `kEnergyElectronic` / `kAtom_f*` with `+=` — the hosting driver
/// owns zeroing. Two epilogues then mean two contributions: ~2× the electronic
/// energy and ~2× the QM force on every ASPC step whose corrector misses.
/// Overwrite (`=`, the standalone entries) was self-healing here; accumulation
/// is not, so the implementation snapshots those channels before the corrector
/// and restores them before the re-solve (`:574-599`, `:602-617`). This file is
/// the pin on that restore. It is a **guard pin, not a RED**: the fix landed
/// with the flip, and no test covered the path.
///
/// ## What a regression would look like
///
/// Measured on the trajectory below (see "Measured landscape"), the corrector's
/// own epilogue contribution is `−30.3936350457 eV` — numerically the same
/// energy as the full solve, since 2 iterations from an ASPC guess are already
/// converged in energy and only fail the deliberately unreachable *threshold*.
/// Dropping the restore therefore puts `−60.7872700915 eV` in `kEnergy`
/// (and ~2.48 eV/Å of force on each H): `kEnergy` would sit 30.4 eV — more than
/// 10¹⁶ tolerances — from `energy() × kHartreeInEv`.
///
/// Which assertion fires:
///   - `EnergyLandsExactlyOnceOnTheFallbackStep` — `kEnergy` vs the method's own
///     `energy()` readback: off by 30.4 eV. **This is the 2× catcher.**
///   - `FallbackStepReproducesThePlainBomdTwin` — off by 30.4 eV.
///   - `ForcesLandExactlyOnceOnTheFallbackStep` — off by ~1.24 eV/Å per atom.
///   - `ElectronicShareIsRestoredWithTheTotal` — would still pass (both channels
///     double together), which is exactly why it is a *separate* concern: it
///     catches a restore that heals `kEnergy` and forgets the per-term key.
///
/// ## Reaching the path from the public API only
///
/// Three legs of the same 6-step nested trajectory identify it without reading
/// a single private member:
///
/// | leg | `corrector_tol` | `fallback` | SCF iterations per step | last step |
/// |---|---|---|---|---|
/// | A | — (no predictor) | — | 3 3 3 3 3 3 | converged |
/// | B | 1e-4 (default) | on | 3 3 3 3 **2 2** | converged |
/// | C | 1e-30 | **off** | 3 3 3 3 **2! 2!** | **NOT converged** |
/// | D | 1e-30 | on | 3 3 3 3 3 3 | converged |
///
/// Leg B proves ASPC *prediction* is active from step 5 (a 2-iteration
/// converged step is the corrector budget, unreachable for the cold/warm full
/// solve, which needs 3). Leg C proves that at `corrector_tol = 1e-30` the
/// corrector *misses*: it stops at the corrector budget, unconverged. Leg D is
/// therefore the double-epilogue path — same tolerance as C, fallback on, and
/// the step ends converged after the full re-solve. `1e-30 Ha` is unreachable
/// by construction: consecutive SCF energies of an O(1 Ha) system differ by
/// ≥ 1 ulp ≈ 2 × 10⁻¹⁶ Ha unless bit-identical, and the criterion's second
/// gate (‖DIIS error‖ < √tol = 1e-15) is not met either — leg C measures the
/// miss rather than assuming it.
///
/// ## Ballistic nuclei (the twin comparison technique)
///
/// `compute_force = false` for the first five steps, so the nuclei move at
/// constant velocity and the ASPC leg and the plain-BOMD twin walk a
/// **bit-identical** geometry sequence (measured: Δx = 0.000e+00 Å after 6
/// steps). The force write is switched on for step 6 only — the step whose
/// accumulators are asserted — so both legs evaluate the gradient at the same
/// geometry and the residual difference is ULP noise, not trajectory drift.
/// Borrowed from tests/cpu/scf/test_scf_xl_bomd.cpp (same file header,
/// technique 1).
///
/// ## Measured landscape (frozen 2026-07-27, this exact configuration)
///
/// H₂/STO-3G, `tol = 1e-10` Ha, dt = 0.5 fs, symmetric stretch v = ±0.004 Å/fs,
/// 6 nested steps, `OMP_NUM_THREADS=1`:
///
/// | quantity | leg D (fallback) | leg A (plain BOMD twin) |
/// |---|---|---|
/// | `kEnergy` [eV] | −30.39363504574894 | −30.39363504574894 |
/// | `energy() × kHartreeInEv` [eV] | −30.39363504574894 | same |
/// | `kEnergyElectronic` [eV] | −30.39363504574894 | −30.39363504574894 |
/// | `f_x` on H₁ [eV/Å] | +1.24238103199911 | +1.24238103199911 |
///
/// `kEnergy − energy()×kHartreeInEv = +0.000e+00` and leg D − leg A =
/// +0.000e+00 on both the energy and the force: bit-identical. (Between other
/// leg pairs the residue reaches 1 × 10⁻¹⁴ eV/Å — one ulp — which is why the
/// tolerances below are ULP-scale rather than literally zero.)
///
/// **Determinism.** Fixed geometry and velocities, no RNG, no wall clock, no
/// filesystem, no network, no third-party package at build or test time.
/// `set_threads(1)` on every driver (and `OMP_NUM_THREADS=1` /
/// `OPENBLAS_NUM_THREADS=1` pinned in CMake) because the claims are ULP-level.
///
/// **Scope.** Unit: four 6-step H₂/STO-3G nested trajectories (2 AOs), ~20 ms
/// total. No energy-drift or iteration-count physics — that is
/// `regressions/dynamics/xl-bomd-rhf-water` territory.

#ifdef ATV_HAS_HF

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/velocity_verlet.h"
#include "atomiverse/cpu/frame/box.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

namespace {

/// @brief Nuclear MD timestep [fs].
constexpr double kDtFs = 0.5;

/// @brief Converged-SCF tolerance |ΔE| [Hartree] for every leg here.
constexpr double kScfTol = 1e-10;

/// @brief H–H separation at t = 0 [Å] — CCCBDB HF/STO-3G r_e at table
///        precision (NIST Standard Reference Database 101). Only the *identity*
///        of the geometry matters here, not its energy.
constexpr double kH2Bond = 0.7122;

/// @brief Nested steps per leg. ASPC(3) needs K+1 = 4 recorded densities, so
///        prediction starts at step 5 (measured); 6 steps puts the asserted
///        step safely inside the predicting regime with one step of margin.
constexpr int kSteps = 6;

/// @brief ASPC extrapolation order K (Kolafa; the shipped default).
constexpr int kAspcOrder = 3;

/// @brief Corrector iteration budget — the shipped default, and the iteration
///        count that fingerprints "the corrector ran" (a full solve needs 3).
constexpr int kCorrectorSteps = 2;

/// @brief A corrector tolerance [Hartree] no double-precision SCF can meet.
///
/// Consecutive energies of an O(1 Ha) system differ by ≥ 1 ulp ≈ 2 × 10⁻¹⁶ Ha
/// unless bit-identical, and `ScfConverged` also demands
/// ‖DIIS error‖ < √tol = 1e-15. Leg C *measures* the resulting miss rather than
/// trusting this argument.
constexpr double kUnreachableCorrectorTol = 1e-30;

/// @brief Energy agreement [eV] asserted between channels that must carry the
///        same single contribution.
///
/// Measured separation of every compared pair is 0.0 — bit-identical, both
/// `kEnergy` vs `energy()` and the fallback leg vs the plain-BOMD twin. One ulp
/// at 30 eV is 3.6 × 10⁻¹⁵ eV, and the failure this guards against is +30.4 eV,
/// so the pin sits at 10⁻¹² eV: three orders above the ULP floor and 10¹⁶ below
/// the bug.
constexpr double kEnergyTol = 1e-12;

/// @brief Force agreement [eV/Å] between the fallback leg and its twin.
///
/// Project FP table, force/exact tier (1e-8). Measured 0.0 for that pair (1 ulp
/// = 2 × 10⁻¹⁶ eV/Å at this magnitude, and other leg pairs do differ by up to
/// 1 × 10⁻¹⁴); the double-count signature is 1.24 eV/Å per atom.
constexpr double kForceTol = 1e-8;

/// @brief H₂ at `kH2Bond`, centred in a 20 Å cubic box, with a symmetric
///        stretch velocity so the geometry changes every step (a frozen
///        geometry would make the ASPC extrapolation exact and the corrector's
///        behaviour uninteresting).
/// @return Two-atom Frame (Z = 1, m = 1.008 amu, v = ∓0.004 Å/fs on x).
cpu::Frame make_h2() {
    cpu::Frame f(2);
    const double centre = 10.0; // box centre: no wrap in the printed geometry
    f.span(kAtom_x)[0] = centre - 0.5 * kH2Bond;
    f.span(kAtom_x)[1] = centre + 0.5 * kH2Bond;
    f.span(kAtom_type)[0] = 1; // Z = 1
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_vx)[0] = -0.004;
    f.span(kAtom_vx)[1] = 0.004;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

/// @brief Readback of one nested trajectory at its final step.
struct LegResult {
    std::vector<int> iterations;  ///< SCF iterations, per step.
    std::vector<bool> converged;  ///< Convergence verdict, per step.
    double energy_ev = 0;         ///< `kEnergy` after the last step [eV].
    double electronic_ev = 0;     ///< `kEnergyElectronic` [eV].
    double method_energy_ev = 0;  ///< `energy() × kHartreeInEv` [eV].
    double fx0 = 0;               ///< f_x on atom 0 after the last step [eV/Å].
    double fx1 = 0;               ///< f_x on atom 1 [eV/Å].
    double x0 = 0;                ///< x of atom 0 after the last step [Å].
};

/// @brief Drive `kSteps` nested MD steps of RHF/STO-3G H₂ and report the final
///        step's accumulators.
///
/// Ballistic for the first `kSteps - 1` steps (`compute_force = false`), then
/// the gradient is switched on for the last step alone — so every leg walks the
/// same geometry sequence and the asserted step is the only one that writes
/// forces. One call site per leg (A-D), in five tests.
/// @param use_predictor   Attach the ASPC density predictor.
/// @param corrector_tol   Corrector |ΔE| threshold [Hartree]; ignored when
///                        `use_predictor` is false.
/// @param fallback        Re-solve at full tolerance when the corrector misses.
/// @param gradient_on_last Write forces on the final step.
/// @return Per-step iteration trace + the final step's channels.
LegResult drive_nested_trajectory(bool use_predictor, double corrector_tol,
                                  bool fallback, bool gradient_on_last) {
    cpu::Frame frame = make_h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = kScfTol;
    rhf.max_iter = 200;
    rhf.compute_force = false;
    if (use_predictor)
        rhf.predictor(kAspcOrder, kCorrectorSteps, corrector_tol, fallback);

    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{kDtFs});
    sim.on(md::Phase::force).add_ref(rhf);
    sim.set_threads(1);
    sim.build();

    LegResult out;
    for (int step = 1; step <= kSteps; ++step) {
        if (gradient_on_last && step == kSteps)
            rhf.compute_force = true;
        sim.run(frame, 1);
        out.iterations.push_back(rhf.iterations());
        out.converged.push_back(rhf.converged());
    }
    out.energy_ev = frame.meta.get(kEnergy);
    out.electronic_ev = frame.meta.get(kEnergyElectronic);
    out.method_energy_ev = rhf.energy() * units::kHartreeInEv;
    out.fx0 = frame.span(kAtom_fx)[0];
    out.fx1 = frame.span(kAtom_fx)[1];
    out.x0 = frame.span(kAtom_x)[0];
    return out;
}

} // namespace

// ─── the path: prediction is active, and the corrector misses ──────────

/// Leg B. With the shipped `corrector_tol = 1e-4` the predicted step converges
/// **inside the corrector budget** (2 iterations, `kCorrectorSteps`), which the
/// full solve cannot do on this trajectory (3 iterations, leg A). That is the
/// public-API witness that ASPC prediction is engaged at the asserted step —
/// without it, a "miss" in leg C could just be a plain solve cut short.
TEST(MeanFieldMethodAspcFallback, PredictionIsActiveOnTheMeasuredStep) {
    const LegResult plain =
        drive_nested_trajectory(false, 0.0, false, /*gradient_on_last=*/false);
    const LegResult predicted = drive_nested_trajectory(
        true, /*corrector_tol=*/1e-4, true, /*gradient_on_last=*/false);

    ASSERT_EQ(static_cast<int>(predicted.iterations.size()), kSteps);
    EXPECT_TRUE(predicted.converged.back());
    EXPECT_EQ(predicted.iterations.back(), kCorrectorSteps)
        << "the last step did not stop at the corrector budget, so ASPC never "
           "predicted there and this file's fallback legs would be measuring a "
           "plain solve (plain-BOMD twin took "
        << plain.iterations.back() << " iterations)";
    EXPECT_GT(plain.iterations.back(), kCorrectorSteps)
        << "a full solve now also finishes within the corrector budget, so the "
           "iteration count no longer fingerprints the corrector — re-derive "
           "the leg table in this file's header";
}

/// Leg C. At `corrector_tol = 1e-30` the corrector cannot converge, and with
/// the fallback switched **off** the step is simply left unconverged at the
/// corrector budget. This is the direct measurement that the fallback branch of
/// `drive_scf_` is entered in leg D (same tolerance, fallback on) — the
/// precondition `!converged()` is observed here rather than assumed.
///
/// It also measures the *other* precondition of the double-count: the
/// epilogue is the tail of `drive_scf_core_` and runs whether or not the loop
/// converged, so the unconverged corrector still deposits a full electronic
/// energy in `kEnergy`. That deposit is what leg D's fallback must undo, and
/// measuring it here is what turns "a second epilogue would double the energy"
/// from a code-reading claim into a number (−30.39 eV, i.e. the doubled total
/// would be −60.79 eV).
TEST(MeanFieldMethodAspcFallback, CorrectorMissesTheUnreachableToleranceWithoutFallback) {
    const LegResult leg = drive_nested_trajectory(
        true, kUnreachableCorrectorTol, /*fallback=*/false,
        /*gradient_on_last=*/false);

    EXPECT_FALSE(leg.converged.back())
        << "the corrector met |dE| < " << kUnreachableCorrectorTol
        << " Ha — the tolerance is no longer unreachable, so the fallback legs "
           "of this file never exercise the double-epilogue path";
    EXPECT_EQ(leg.iterations.back(), kCorrectorSteps)
        << "the corrector did not stop at its budget";
    EXPECT_NEAR(leg.energy_ev, leg.method_energy_ev, kEnergyTol)
        << "the unconverged corrector did not write its energy to kEnergy — if "
           "the corrector's epilogue no longer accumulates, the guard in "
           "mean_field_method.cpp:574-617 is dead code and the legs below stop "
           "meaning what this file says they mean";
    EXPECT_LT(leg.energy_ev, -1.0)
        << "the corrector's contribution is not a full electronic energy";
}

// ─── the pin: exactly one contribution per entry ───────────────────────

/// **The 2× catcher.** On the fallback step `drive_scf_core_` runs twice, so its
/// epilogue accumulates twice unless the guard restores the channels between
/// the two calls. `energy()` is a plain member overwritten by the last solve, so
/// it is the independent witness of what ONE contribution is worth: `kEnergy`
/// must equal it (bit-identical in measurement), not twice it.
TEST(MeanFieldMethodAspcFallback, EnergyLandsExactlyOnceOnTheFallbackStep) {
    const LegResult leg =
        drive_nested_trajectory(true, kUnreachableCorrectorTol, /*fallback=*/true,
                                /*gradient_on_last=*/true);

    ASSERT_TRUE(leg.converged.back())
        << "the fallback did not re-solve to convergence, so this step is not "
           "the double-epilogue path";
    EXPECT_GT(leg.iterations.back(), kCorrectorSteps)
        << "the reported iteration count is still the corrector's — the full "
           "re-solve never ran";
    EXPECT_NEAR(leg.energy_ev, leg.method_energy_ev, kEnergyTol)
        << "kEnergy = " << leg.energy_ev << " eV but the method's own energy is "
        << leg.method_energy_ev
        << " eV. A ratio near 2 is the ASPC double-epilogue signature: the "
           "corrector's `+=` was not undone before the fallback re-solve "
           "(mean_field_method.cpp:574-617), so this method contributed twice "
           "to a driver-owned accumulator in one launch";
    EXPECT_LT(std::fabs(leg.energy_ev / leg.method_energy_ev - 1.0), 1e-12)
        << "kEnergy / energy() = " << leg.energy_ev / leg.method_energy_ev
        << " (2.0 would be the exact double-count)";
}

/// The per-term key must be restored in lockstep with the total: a guard that
/// heals `kEnergy` and forgets `kEnergyElectronic` leaves the decomposition
/// invariant (`kEnergy == Σ per-term + unclassified`) broken by one whole
/// electronic energy — and *no* total-energy assertion can see it, because both
/// channels would double together in the unguarded code. Hence a separate test.
TEST(MeanFieldMethodAspcFallback, ElectronicShareIsRestoredWithTheTotal) {
    const LegResult leg =
        drive_nested_trajectory(true, kUnreachableCorrectorTol, /*fallback=*/true,
                                /*gradient_on_last=*/true);

    ASSERT_TRUE(leg.converged.back());
    EXPECT_NEAR(leg.electronic_ev, leg.energy_ev, kEnergyTol)
        << "kEnergyElectronic = " << leg.electronic_ev << " eV but kEnergy = "
        << leg.energy_ev
        << " eV. RHF is the only source in this loop, so the electronic share "
           "IS the total; a mismatch means the fallback restored one channel "
           "and not the other";
}

/// The force channel is guarded by the same snapshot/restore, on the
/// `opts.forces` branch of the epilogue. Compared against the plain-BOMD twin
/// at a bit-identical geometry (ballistic nuclei), so the only admissible
/// difference is ULP noise — a doubled gradient is ~1.24 eV/Å per atom here.
TEST(MeanFieldMethodAspcFallback, ForcesLandExactlyOnceOnTheFallbackStep) {
    const LegResult leg =
        drive_nested_trajectory(true, kUnreachableCorrectorTol, /*fallback=*/true,
                                /*gradient_on_last=*/true);
    const LegResult twin =
        drive_nested_trajectory(false, 0.0, false, /*gradient_on_last=*/true);

    ASSERT_TRUE(leg.converged.back());
    ASSERT_TRUE(twin.converged.back());
    ASSERT_NEAR(leg.x0, twin.x0, 1e-14)
        << "the two legs did not walk the same geometry sequence, so their "
           "forces are not comparable (the ballistic-nuclei technique broke)";
    ASSERT_GT(std::fabs(twin.fx0), 0.1)
        << "the twin's force is ~0, so a doubled force would be invisible here";

    EXPECT_NEAR(leg.fx0, twin.fx0, kForceTol)
        << "f_x[0] = " << leg.fx0 << " eV/Å on the fallback step but "
        << twin.fx0
        << " eV/Å on the plain-BOMD twin at the same geometry — a factor ~2 is "
           "the corrector gradient left in the accumulator";
    EXPECT_NEAR(leg.fx1, twin.fx1, kForceTol);
}

/// Equivalence, the other way round: the fallback step must be
/// indistinguishable from the plain full solve it falls back *to*. The
/// implementation restores the snapshot rather than subtracting the corrector's
/// contribution precisely so that this holds bit-for-bit; a subtract-based
/// undo would show up here as ULP-scale-but-not-ULP-clean drift, and a missing
/// undo as 30.4 eV.
TEST(MeanFieldMethodAspcFallback, FallbackStepReproducesThePlainBomdTwin) {
    const LegResult leg =
        drive_nested_trajectory(true, kUnreachableCorrectorTol, /*fallback=*/true,
                                /*gradient_on_last=*/true);
    const LegResult twin =
        drive_nested_trajectory(false, 0.0, false, /*gradient_on_last=*/true);

    ASSERT_TRUE(leg.converged.back());
    ASSERT_TRUE(twin.converged.back());
    EXPECT_NEAR(leg.energy_ev, twin.energy_ev, kEnergyTol)
        << "the ASPC fallback step reported " << leg.energy_ev
        << " eV where a plain full solve at the same geometry reports "
        << twin.energy_ev << " eV";
    EXPECT_NEAR(leg.electronic_ev, twin.electronic_ev, kEnergyTol);
}

#endif // ATV_HAS_HF
