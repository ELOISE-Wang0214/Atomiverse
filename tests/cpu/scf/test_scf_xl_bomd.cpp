/// @file tests/cpu/scf/test_scf_xl_bomd.cpp
/// @brief esmm-cpmd-05-scf-xlbomd — extended-Lagrangian BOMD (Niklasson) as
///        the Gaussian-SCF initial guess: the `xl_bomd` member surface, the
///        exclusivity gates, auxiliary-history integrity, and the single
///        ordered guess ladder.
///
/// Mirrors one production unit: `include/atomiverse/cpu/scf/mean_field_method.h`
/// (`MeanFieldMethod<Derived>` — the `xl_bomd` / `xl_bomd_k` members, the
/// `xl_advances()` readback, the ladder inside `drive_scf_`, and the
/// advance-on-converged-nested-run hook). Two suites, one unit:
///   - `MeanFieldMethodXlBomd` — member surface, refusals, history integrity;
///   - `MeanFieldMethodGuessLadder` — the ac-001 resolution order
///     `explicit guess > xl aux() > ASPC predict > last-P reuse`.
///
/// **The propagator itself** (`electronic::XlBomd` — the dissipative Verlet
/// recursion and its Table-I coefficients, Niklasson et al., J. Chem. Phys.
/// **130**, 214109 (2009)) is method-agnostic dense algebra pinned by
/// tests/cpu/electronic/test_xl_bomd.cpp; nothing here re-tests it. What is
/// pinned here is only the SCF side of the seam.
///
/// **Scope.** Unit only. The XL-BOMD *physics* — the warm-step iteration
/// margin averaged over a long trajectory and the energy-drift envelope that
/// makes XL worth having — is not a seconds-scale claim and is pinned by
/// `regressions/dynamics/xl-bomd-rhf-water` against Niklasson (2009). The
/// whole file is H₂O/STO-3G (7 AOs, 10 electrons); all 14 tests together took
/// 0.23 s at OMP_NUM_THREADS=8 when dry-run against a no-op XL stub.
///
/// ## Two techniques this file uses, and why
///
/// **1. Ballistic nuclei (`compute_force = false`).** The ladder is about
/// *which density seeds the SCF*, so every other input must be held identical
/// between the legs of a comparison. With the force write suppressed the
/// nuclei move at constant velocity, so two twin runs see a **bit-identical**
/// geometry sequence (measured: max|Δx| = 0.0 Å after 9 steps) instead of
/// drifting apart through their slightly different converged forces. The SCF
/// still runs a full converged solve per step and still warm-starts from the
/// previous step — the only thing removed is the confound.
///
/// **2. The one-shot energy fingerprint (`tol = kOneShotTol`).** A converged
/// energy is guess-independent by construction, so it cannot witness *which*
/// guess was used. An |ΔE| threshold of 1e6 Ha is met by the first energy
/// difference the loop measures, so the step stops after 2 iterations and
/// `energy()` becomes a deterministic function of the seeding density alone.
/// Measured separation on this trajectory: two different warm-guess sources
/// give one-shot energies 5.6e-07 Ha apart, while re-running the *same* source
/// reproduces its value to ≤ 4e-14 Ha (1 and 8 OpenMP threads alike) — seven
/// orders of margin either side of the 1e-9 / 1e-10 Ha thresholds used below.
///
/// ## Measured landscape (frozen 2026-07-27, this exact configuration)
///
/// H₂O/STO-3G, ballistic nuclei, dt = 0.5 fs, `tol` = 1e-10, RHF SCF
/// iterations per nested step 1..9:
///   - last-P reuse:      10, 9, 7, 9, 9, 9, 9, 9, 9
///   - ASPC(3) predictor: 10, 9, 7, 9, 2, 2, 2, 2, 2
/// i.e. cold seed 10, plain warm start 9, and an *extrapolating* guess 2 once
/// its history fills. UHF{5,5} on the same trajectory: cold 21, last-P 13,
/// ASPC 2. XL-BOMD is an extrapolating guess of the same family, so the
/// assertions below are `< cold` and `<= last-P`: the strict warm-step margin
/// is a statistical property of a long run and belongs to the regression
/// (ac-004), not to a single frozen step of a unit test. Identical counts were
/// measured at OMP_NUM_THREADS=1 and 8.
///
/// ## New public readback required by ac-003
///
/// `electronic::XlBomd` is private workspace (Frame purity: no propagator on
/// the public surface), so "the aux history advanced exactly once per
/// converged nested run" needs one scalar witness, exactly as
/// `soscf_iterations()` surfaces the private CIAH stepper's macro count:
///
///     /// @brief XL-BOMD auxiliary-density propagation steps taken (one per
///     ///        converged nested run). @return Count; 0 under BOMD.
///     int xl_advances() const;
///
/// Behavioural proxies alone cannot close ac-003: an implementation that never
/// advances and one that advances on every entry (standalone runs included)
/// are indistinguishable from iteration counts.
///
/// **Determinism.** Fixed geometry and velocities, no RNG, no wall clock, no
/// filesystem, no network, no third-party oracle at test time. The single
/// hard-coded golden is the H₂O/STO-3G RHF energy already frozen in-tree
/// (tests/cpu/scf/test_method_driver_parity.cpp:156, captured 2026-07-20).

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/velocity_verlet.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/methods.h"
#include "atomiverse/cpu/scf/scf_guess.h"
#include "atomiverse/cpu/scf/slater_determinant.h"

using namespace ATV;

using cpu::linalg::RealMatrix;
using cpu::scf::BasisSet;
using cpu::scf::SCFInitialGuess;
using cpu::scf::SlaterDeterminant;
using Spin = cpu::scf::SlaterDeterminant::Spin;

namespace {

/// Nuclear MD timestep [fs].
constexpr double kDtFs = 0.5;

/// Converged-SCF tolerance |ΔE| [Hartree] for every trajectory here.
constexpr double kScfTol = 1e-10;

/// Warm-up nested steps before a measurement. K = 5 keeps a depth-6 history,
/// so 8 propagations have fully replaced the stationary seed.
constexpr int kWarmSteps = 8;

/// |ΔE| threshold [Hartree] that stops the SCF at once, turning a step into a
/// one-shot evaluation of its initial guess (see the file header).
constexpr double kOneShotTol = 1e6;

/// Two one-shot energies [Hartree] this close came from the same guess
/// (measured reproducibility of one source: ≤ 4e-14 Ha across thread counts).
constexpr double kSameGuessTol = 1e-10;

/// Two one-shot energies [Hartree] at least this far apart came from
/// different guesses (measured separation of two warm sources: 5.6e-07 Ha).
constexpr double kGuessSeparation = 1e-9;

/// Converged-energy agreement [Hartree] between two guess strategies at the
/// same geometry — the SCF fixed point does not know how it was seeded
/// (project FP table: energy, exact tier).
constexpr double kFixedPointTol = 1e-10;

/// H₂O/STO-3G RHF total energy [Hartree] at the reference geometry below.
/// Hard-coded golden, in-tree provenance: tests/cpu/scf/test_method_driver_parity.cpp:156
/// (`kRhfWaterEnergy`, frozen 2026-07-20 from the driver-ized RHF that spec
/// scf-driver-refactor-03a proved bit-identical to the retired solver).
constexpr double kPinnedWaterRhf = -74.9629648020;

/// Pinned-energy tolerance [Hartree] (that pin was captured single-threaded;
/// this binary runs at OMP_NUM_THREADS=8, whose reduction order costs ~1e-13).
constexpr double kPinnedEnergyTol = 1e-9;

/// @brief H₂O / STO-3G at the canonical in-tree geometry (r = 0.9572 Å,
///        ∠HOH = 104.52°) — the geometry of the pinned RHF energy above.
/// @param moving When true, give the hydrogens a sub-thermal symmetric-stretch
///        velocity (±0.008, 0, +0.006) Å/fs with the oxygen recoiling at
///        −0.000756 Å/fs so the net momentum is zero.
/// @return Three-atom Frame (Z = 8, 1, 1; masses in amu; 50 Å vacuum box).
cpu::Frame make_water(bool moving) {
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
    if (moving) {
        f.span(kAtom_vx)[1] = 0.008;
        f.span(kAtom_vz)[1] = 0.006;
        f.span(kAtom_vx)[2] = -0.008;
        f.span(kAtom_vz)[2] = 0.006;
        f.span(kAtom_vz)[0] = -0.000756;
    }
    f.box = cpu::Box::cubic(50.0);
    return f;
}

/// @brief RHF/STO-3G with XL-BOMD switched on and nothing else unusual.
/// @return The method by value (the caller owns it; drivers borrow it).
cpu::scf::RHF xl_rhf() {
    cpu::scf::RHF rhf;
    rhf.basis = BasisSet::STO3G;
    rhf.tol = kScfTol;
    rhf.max_iter = 200;
    rhf.xl_bomd = true;
    return rhf;
}

/// @brief The plain Born–Oppenheimer twin: same method, last-P reuse only.
cpu::scf::RHF bomd_rhf() {
    cpu::scf::RHF rhf = xl_rhf();
    rhf.xl_bomd = false;
    return rhf;
}

/// @brief Bind `method` at the MD force phase and take exactly one step.
///
/// The uniform driver for every refusal row (the same helper the CPMD suite
/// uses, tests/cpu/scf/test_scf_cpmd.cpp): validation may fire at the nested
/// `build_on` or at the first nested drive, and both sites are inside this
/// call, so a row's `EXPECT_THROW` does not depend on which one the
/// implementation picks.
/// @param method The concrete mean-field method under test (moved in).
/// @param frame  Geometry driven for one step.
template <typename Method>
void drive_nested_md_step(Method method, cpu::Frame& frame) {
    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{kDtFs});
    sim.on(md::Phase::force).add_ref(method);
    sim.build();
    sim.run(frame, 1);
}

/// Readback of one ballistic trajectory (see the file header for why the
/// nuclei are ballistic and what the one-shot energy measures).
struct TrajectoryProbe {
    int warm_iterations = 0; ///< SCF iterations of the last converged step.
    double warm_energy = 0;  ///< Converged energy of that step [Hartree].
    double probe_energy = 0; ///< One-shot energy of the extra step [Hartree];
                             ///< 0 when no probe step was taken.
};

/// @brief Drive a borrowed `method` through `warm_steps` nested MD steps on
///        the ballistic water trajectory, optionally followed by one extra
///        step at `kOneShotTol` whose energy fingerprints the guess.
///
/// One driver spans both phases, so the probe step continues the *same*
/// trajectory (and the same nested wiring) the warm steps built. Mutates two
/// members of the borrowed method: `compute_force` (off, for the ballistic
/// trajectory) and, when `probe`, `tol` (left at `kOneShotTol`).
/// @param method     Borrowed method; its warm state survives the call.
/// @param warm_steps Converged nested steps to take (≥ 1).
/// @param probe      Take the extra one-shot step afterwards.
/// @return Iteration count + energies (see `TrajectoryProbe`).
template <typename Method>
TrajectoryProbe drive_ballistic_md(Method& method, int warm_steps, bool probe) {
    method.compute_force = false; // zero forces ⇒ identical geometry sequence
    cpu::Frame frame = make_water(/*moving=*/true);
    TrajectoryProbe out;
    {
        cpu::md::MD sim;
        sim.set_integrator(cpu::md::VelocityVerlet{kDtFs});
        sim.on(md::Phase::force).add_ref(method);
        sim.build();
        sim.run(frame, warm_steps);
        out.warm_iterations = method.iterations();
        out.warm_energy = method.energy();
        if (probe) {
            method.tol = kOneShotTol;
            sim.run(frame, 1);
            out.probe_energy = method.energy();
        }
    }
    return out;
}

/// @brief Restricted total density 2·C_occ C_occᵀ of a published determinant.
/// @param det Converged wavefunction (restricted).
/// @return P [dimensionless], the convention the SCF guess expects.
RealMatrix total_density(const SlaterDeterminant& det) {
    const RealMatrix& C = det.coefficients(Spin::alpha);
    const int n_occ = det.n_electrons(Spin::alpha);
    const std::size_t n = C.shape()[0];
    RealMatrix P = xt::zeros<double>({n, n});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (int o = 0; o < n_occ; ++o)
                s += C(i, static_cast<std::size_t>(o))
                     * C(j, static_cast<std::size_t>(o));
            P(i, j) = 2.0 * s;
        }
    return P;
}

/// @brief A converged density at a *displaced* geometry — a guess that is
///        good enough to run but distinguishable from any warm state built on
///        the reference trajectory.
/// @return Restricted α/β density guess.
SCFInitialGuess displaced_geometry_guess() {
    cpu::Frame f = make_water(/*moving=*/false);
    f.span(kAtom_x)[1] += 0.05;
    cpu::scf::RHF source;
    source.basis = BasisSet::STO3G;
    source.tol = kScfTol;
    source.max_iter = 200;
    source.run(f);
    const RealMatrix P = total_density(source.determinant());
    SCFInitialGuess guess;
    guess.P_alpha = P;
    guess.P_beta = P;
    guess.has_alpha = true;
    guess.has_beta = true;
    return guess;
}

/// @brief One standalone one-shot run of `method` at the reference geometry,
///        seeded by `guess`.
/// @param method Borrowed method (its internal warm state is what we probe).
/// @param guess  Explicit density guess handed to `run`.
/// @return Energy after the one-shot step [Hartree].
template <typename Method>
double one_shot_with_guess(Method& method, const SCFInitialGuess& guess) {
    method.tol = kOneShotTol;
    cpu::Frame f = make_water(/*moving=*/false);
    method.run(f, &guess);
    return method.energy();
}

} // namespace

// ─── the XL-BOMD member surface ────────────────────────────────────────

/// Parameters are **members with defaults**, not an `XlConfig` bag and not
/// constructor arguments: a plain BOMD user must be able to ignore all of
/// them, and an XL user must get the same dissipation order the GPW evaluator
/// ships (`cpu/gpw/gpw_evaluator.h`: `xl_bomd`, `xl_bomd_k = 5`). The
/// readback starts at zero, before any auxiliary density exists.
TEST(MeanFieldMethodXlBomd, MemberSurfaceDefaults) {
    cpu::scf::RHF rhf;
    EXPECT_FALSE(rhf.xl_bomd);
    EXPECT_EQ(rhf.xl_bomd_k, 5); // Niklasson Table I row K = 5
    EXPECT_EQ(rhf.xl_advances(), 0);
}

// ─── ac-002: one guess owner at a time ─────────────────────────────────

/// ASPC extrapolates a polynomial through past densities; XL-BOMD propagates
/// an auxiliary density with its own equation of motion. Both answer "what is
/// the guess now", so running both means two answers — the same pair GPW
/// refuses (`GPWEvaluator::build`) and the same shape as the CPMD × predictor
/// row in tests/cpu/scf/test_scf_cpmd.cpp.
TEST(MeanFieldMethodXlBomd, RefusesAspcPredictor) {
    auto f = make_water(/*moving=*/false);
    auto m = xl_rhf();
    m.predictor(3);
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

/// The other half of the lattice esmm-cpmd-04 left open: its
/// `validate_cpmd_config_` refuses a competing *density* propagator, but
/// `xl_bomd` did not exist yet, so tests/cpu/scf/test_scf_cpmd.cpp records the
/// `xl_bomd` row as deliberately NOT pinned (file header, "Refusal lattice
/// coverage"). This test closes that deferred row. Both flags live on one
/// method object, so whichever validator runs first sees the pair — one throw
/// covers both directions of the exclusion.
TEST(MeanFieldMethodXlBomd, RefusesCpmd) {
    auto f = make_water(/*moving=*/false);
    auto m = xl_rhf();
    m.cpmd = true; // VelocityVerlet publishes the timestep CPMD needs, so this
                   // row fails on the exclusion and not on the missing clock
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

/// The dissipative recursion's (κ, α, c_k) are tabulated rows, not a formula:
/// Table I of Niklasson et al., JCP **130**, 214109 (2009) covers K ∈ [3, 9]
/// only. `electronic::XlBomd` throws outside it and the method must surface
/// that rather than clamp — house rule: lookup tables throw, no fallbacks
/// (GPW does the same at `build`).
TEST(MeanFieldMethodXlBomd, RefusesDissipationOrderOutsideTable) {
    auto f = make_water(/*moving=*/false);
    auto m = xl_rhf();
    m.xl_bomd_k = 11;
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

// ─── ac-003: the auxiliary history is time-continuous ──────────────────

/// The aux density is a dynamical variable on the nuclear clock: exactly one
/// propagation per nuclear step, or the recursion is integrating a different
/// trajectory than the nuclei are following.
TEST(MeanFieldMethodXlBomd, AuxAdvancesOncePerConvergedNestedStep) {
    auto m = xl_rhf();
    const TrajectoryProbe t = drive_ballistic_md(m, /*warm_steps=*/3, false);
    ASSERT_TRUE(m.converged()) << "the trajectory must converge every step for "
                                  "this test to say anything";
    EXPECT_GT(t.warm_iterations, 0);
    EXPECT_EQ(m.xl_advances(), 3);
}

/// A standalone `run()` is a single point somewhere else in configuration
/// space — feeding it into the recursion would splice a foreign geometry into
/// a time-continuous history. GPW takes the same line (XL state advances on
/// the MD force path only), and unlike CPMD a standalone entry is *allowed*
/// here: XL is a guess strategy, so the run itself must still work.
TEST(MeanFieldMethodXlBomd, StandaloneRunLeavesAuxUntouched) {
    auto m = xl_rhf();
    drive_ballistic_md(m, /*warm_steps=*/3, false);
    ASSERT_EQ(m.xl_advances(), 3);

    auto f = make_water(/*moving=*/false);
    m.compute_force = true;
    m.tol = kScfTol;
    EXPECT_NO_THROW(m.run(f));
    EXPECT_TRUE(m.converged());
    EXPECT_EQ(m.xl_advances(), 3) << "a standalone single point is not a step "
                                     "of the nuclear trajectory";
}

/// An unconverged SCF has no density worth propagating: recording it would
/// inject the iteration's residual error into every later guess. `max_iter =
/// 1` is the cheapest deterministic way to never converge (measured: 1
/// iteration, `converged() == false`, no throw, on every step).
TEST(MeanFieldMethodXlBomd, UnconvergedNestedRunLeavesAuxUntouched) {
    auto m = xl_rhf();
    m.max_iter = 1;
    drive_ballistic_md(m, /*warm_steps=*/3, false);
    ASSERT_FALSE(m.converged()) << "this test needs a genuinely unconverged "
                                   "nested run";
    EXPECT_EQ(m.xl_advances(), 0);
}

// ─── ac-001: the XL signature, in the form a unit test can own ─────────

/// Warm-starting from the propagated aux must beat the cold seed and must not
/// be *worse* than the plain last-P reuse it replaces. Measured on this exact
/// trajectory: cold seed 10 iterations, last-P warm step 9, an extrapolating
/// guess (ASPC, the same family as XL) 2. The strict XL-vs-last-P margin is a
/// mean over a long run — that claim lives in
/// regressions/dynamics/xl-bomd-rhf-water (ac-004), not in one frozen step.
TEST(MeanFieldMethodXlBomd, WarmStepIteratesFewerTimesThanTheColdSeed) {
    auto cold = xl_rhf();
    const TrajectoryProbe c = drive_ballistic_md(cold, /*warm_steps=*/1, false);
    auto warm = xl_rhf();
    const TrajectoryProbe w = drive_ballistic_md(warm, kWarmSteps + 1, false);
    auto twin = bomd_rhf();
    const TrajectoryProbe t = drive_ballistic_md(twin, kWarmSteps + 1, false);

    ASSERT_TRUE(cold.converged());
    ASSERT_TRUE(warm.converged());
    ASSERT_TRUE(twin.converged());
    EXPECT_LT(w.warm_iterations, c.warm_iterations)
        << "XL warm start is not warmer than a cold seed";
    EXPECT_LE(w.warm_iterations, t.warm_iterations)
        << "XL is slower than the plain last-P reuse it replaces — a "
           "wrong-scale aux (restricted must carry the factor-2 total "
           "density) looks exactly like this";
}

/// The guess is a guess: whichever density seeds it, the SCF must land on the
/// same stationary point. Both legs walk a bit-identical geometry sequence
/// (zero forces), so this is a same-geometry comparison and the exact-energy
/// tier applies. Measured agreement of the two guess paths: 4.3e-14 Ha.
TEST(MeanFieldMethodXlBomd, WarmTrajectoryReachesTheSameFixedPointAsLastP) {
    auto xl = xl_rhf();
    const TrajectoryProbe x = drive_ballistic_md(xl, kWarmSteps + 1, false);
    auto twin = bomd_rhf();
    const TrajectoryProbe t = drive_ballistic_md(twin, kWarmSteps + 1, false);

    ASSERT_TRUE(xl.converged());
    ASSERT_TRUE(twin.converged());
    EXPECT_NEAR(x.warm_energy, t.warm_energy, kFixedPointTol);
}

/// Turning XL on must not perturb a plain single point: no propagation has
/// happened, so the energy is the ordinary Born–Oppenheimer one. Hard-coded
/// golden (see `kPinnedWaterRhf`), which also fixes the geometry every other
/// test in this file drives.
TEST(MeanFieldMethodXlBomd, SinglePointEnergyMatchesThePinnedReference) {
    auto f = make_water(/*moving=*/false);
    auto m = xl_rhf();
    m.run(f);

    ASSERT_TRUE(m.converged());
    EXPECT_NEAR(m.energy(), kPinnedWaterRhf, kPinnedEnergyTol);
    EXPECT_EQ(m.xl_advances(), 0);
}

/// UHF/UKS propagate α and β independently (RHF/RKS propagate the factor-2
/// total), so the unrestricted leg must advance its own channels and still
/// land on the closed-shell fixed point. Measured on this trajectory: UHF cold
/// seed 21 iterations, last-P warm step 13, extrapolating guess 2; the
/// converged energy equals the RHF one, −74.961141829744 Ha.
TEST(MeanFieldMethodXlBomd, UnrestrictedAdvancesAndKeepsTheFixedPoint) {
    cpu::scf::UHF xl{/*n_alpha=*/5, /*n_beta=*/5};
    xl.basis = BasisSet::STO3G;
    xl.tol = kScfTol;
    xl.max_iter = 200;
    xl.xl_bomd = true;
    const TrajectoryProbe x = drive_ballistic_md(xl, kWarmSteps + 1, false);

    cpu::scf::UHF twin{/*n_alpha=*/5, /*n_beta=*/5};
    twin.basis = BasisSet::STO3G;
    twin.tol = kScfTol;
    twin.max_iter = 200;
    const TrajectoryProbe t = drive_ballistic_md(twin, kWarmSteps + 1, false);

    ASSERT_TRUE(xl.converged());
    ASSERT_TRUE(twin.converged());
    EXPECT_EQ(xl.xl_advances(), kWarmSteps + 1);
    EXPECT_LE(x.warm_iterations, t.warm_iterations);
    EXPECT_NEAR(x.warm_energy, t.warm_energy, kFixedPointTol);
}

// ─── ac-001: one ladder, one order ─────────────────────────────────────
//
// `explicit guess > xl aux() > ASPC predict > last-P reuse`. Each rung is
// pinned by a one-shot energy: same value ⇒ the same density seeded the SCF,
// different value ⇒ a different one (file header, technique 2).

/// The rung XL adds. After a warm trajectory the aux density and the last
/// converged density are different matrices, so if the ladder prefers the aux
/// the two legs cannot produce the same one-shot energy. Both legs walk a
/// bit-identical geometry sequence, so nothing but the guess can separate
/// them — the `ASSERT_NEAR` on the converged energies is that precondition,
/// not the point of the test.
TEST(MeanFieldMethodGuessLadder, XlAuxOutranksLastDensityReuse) {
    auto xl = xl_rhf();
    const TrajectoryProbe x = drive_ballistic_md(xl, kWarmSteps + 1, true);
    auto twin = bomd_rhf();
    const TrajectoryProbe t = drive_ballistic_md(twin, kWarmSteps + 1, true);

    ASSERT_NEAR(x.warm_energy, t.warm_energy, kFixedPointTol)
        << "the two legs must reach the same converged state at the same "
           "geometry before their one-shot energies can be compared";
    EXPECT_GT(std::fabs(x.probe_energy - t.probe_energy), kGuessSeparation)
        << "the XL method seeded its step with the same density the plain "
           "last-P method used — the propagated aux never reached the SCF";
}

/// The top rung. A caller who hands `run()` a density has a reason (a chosen
/// SCF branch, a continuous finite-difference path); a warm internal state may
/// not silently replace it. Both legs are the same method with the same
/// settings at the same geometry with the same explicit guess — only the
/// internal warm state differs, so equal one-shot energies mean the explicit
/// guess won on both.
TEST(MeanFieldMethodGuessLadder, ExplicitGuessOutranksXlAux) {
    const SCFInitialGuess guess = displaced_geometry_guess();

    auto warm = xl_rhf();
    drive_ballistic_md(warm, kWarmSteps, false);
    ASSERT_EQ(warm.xl_advances(), kWarmSteps);
    const double e_warm = one_shot_with_guess(warm, guess);

    auto fresh = xl_rhf();
    const double e_fresh = one_shot_with_guess(fresh, guess);

    EXPECT_NEAR(e_warm, e_fresh, kSameGuessTol)
        << "a warm XL aux displaced the caller's explicit guess";
}

/// The same rung one level down, and the one the current tree gets wrong: the
/// ASPC front-end in `drive_scf_` overwrites `effective_guess` whenever it has
/// a prediction, *including* when the caller passed one (measured 2026-07-27:
/// the two legs below differ by 2.8e-05 Ha). Spec esmm-cpmd-05 rewrites that
/// block as one ordered resolution with the explicit guess on top, and a
/// prediction that is not used must not be "corrected" either — so the warm
/// leg has to take the same path as the fresh one.
TEST(MeanFieldMethodGuessLadder, ExplicitGuessOutranksAspcPrediction) {
    const SCFInitialGuess guess = displaced_geometry_guess();

    auto warm = bomd_rhf();
    warm.predictor(3);
    drive_ballistic_md(warm, kWarmSteps, false);
    const double e_warm = one_shot_with_guess(warm, guess);

    auto fresh = bomd_rhf();
    fresh.predictor(3);
    const double e_fresh = one_shot_with_guess(fresh, guess);

    EXPECT_NEAR(e_warm, e_fresh, kSameGuessTol)
        << "a warm ASPC prediction displaced the caller's explicit guess";
}
