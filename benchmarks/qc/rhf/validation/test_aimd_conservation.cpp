/// @file tests/cpu/scf/test_aimd_conservation.cpp
/// @brief Long-trajectory AIMD validation: NVE energy drift + ASPC
///        state-leakage check.
///
/// Together these pin the two AIMD-specific regressions that the short
/// 5-step BOMD tests in `test_scf_evaluator.cpp` cannot catch:
///
///   1. **Energy drift** — accumulated integrator error over a longer
///      trajectory. A pass here means VelocityVerlet + SCF
///      remain symplectic-equivalent over many steps for the chosen dt.
///
///   2. **State leakage via ASPC** — the predictor-corrector reuses
///      density across MD steps. If anything in `SCF` were to
///      cache geometry-dependent state and forget to invalidate it on
///      the next step (a bug seen during early Driver/Component work),
///      ASPC propagation would silently drift away from a fresh
///      single-point at the same geometry.
///
/// Both tests use the H2 STO-3G fixture for sub-second wall time.

#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/frame/materialize.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_types.h"
#include "atomiverse/molrs.hpp"

using namespace ATV;

namespace {

/// @brief H2 fixture matching `BOMDTest::make_h2`. Local copy keeps the
///        two test files independent.
cpu::Frame make_h2(double bond_length, double v_stretch) {
    cpu::Frame f(2);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = bond_length;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_z)[1] = 0.0;
    f.span(ATV::kAtom_vx)[0] = -v_stretch / 2;
    f.span(ATV::kAtom_vx)[1] = v_stretch / 2;
    f.span(ATV::kAtom_vy)[0] = 0.0;
    f.span(ATV::kAtom_vy)[1] = 0.0;
    f.span(ATV::kAtom_vz)[0] = 0.0;
    f.span(ATV::kAtom_vz)[1] = 0.0;
    f.span(ATV::kAtom_fx)[0] = 0.0;
    f.span(ATV::kAtom_fx)[1] = 0.0;
    f.span(ATV::kAtom_fy)[0] = 0.0;
    f.span(ATV::kAtom_fy)[1] = 0.0;
    f.span(ATV::kAtom_fz)[0] = 0.0;
    f.span(ATV::kAtom_fz)[1] = 0.0;
    f.span(ATV::kAtom_mass)[0] = 1.008;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.span(ATV::kAtom_type)[0] = 1;
    f.span(ATV::kAtom_type)[1] = 1;
    f.meta.get(ATV::kEnergy) = 0;
    f.meta.get(ATV::kKinetic) = 0;
    f.simbox = ATV::cpu::Box::cubic(100.0);
    return f;
}

/// @brief Build an H2 molrs.Frame mirroring `make_h2`, ready for the
///        `cpu::md::MD::run(molrs::FrameRef&, ...)` ingress path.
::rust::Box<molrs::FrameRef> make_h2_molrs(double bond_length,
                                           double v_stretch) {
    auto fr = molrs::frame_new();
    const std::vector<double> x = {0.0, bond_length};
    const std::vector<double> y = {0.0, 0.0};
    const std::vector<double> z = {0.0, 0.0};
    const std::vector<double> vx = {-v_stretch / 2, v_stretch / 2};
    const std::vector<double> vy = {0.0, 0.0};
    const std::vector<double> vz = {0.0, 0.0};
    const std::vector<double> mass = {1.008, 1.008};
    const std::vector<std::int32_t> type = {1, 1};
    auto f64 = [](const std::vector<double>& v) {
        return ::rust::Slice<const double>(v.data(), v.size());
    };
    molrs::frame_set_column_f64(*fr, "atoms", "x", f64(x));
    molrs::frame_set_column_f64(*fr, "atoms", "y", f64(y));
    molrs::frame_set_column_f64(*fr, "atoms", "z", f64(z));
    molrs::frame_set_column_f64(*fr, "atoms", "vx", f64(vx));
    molrs::frame_set_column_f64(*fr, "atoms", "vy", f64(vy));
    molrs::frame_set_column_f64(*fr, "atoms", "vz", f64(vz));
    molrs::frame_set_column_f64(*fr, "atoms", "mass", f64(mass));
    molrs::frame_set_column_i32(
        *fr, "atoms", "type",
        ::rust::Slice<const std::int32_t>(type.data(), type.size()));
    const double h[9] = {100.0, 0, 0, 0, 100.0, 0, 0, 0, 100.0};
    molrs::frame_set_simbox(*fr, ::rust::Slice<const double>(h, 9));
    return fr;
}

} // namespace

// ── NVE energy-drift test ───────────────────────────────────────────
//
// We sample E_total = PE + KE every 10 steps over a 200-step run and
// assert that the **peak-to-peak envelope** stays within a bound. For a
// truly conservative integrator the envelope reflects only initial
// kinetic+potential excitation (which then sloshes back and forth
// between PE and KE without ever leaving the bound). Long-term
// integrator drift would push the envelope wider than the initial
// excitation, which is what this test is here to catch.
//
// Comparing instantaneous E at two arbitrary trajectory steps would
// pick up the phase of the bond oscillation rather than true drift —
// a 0.001 Å/fs initial velocity already moves PE by tens of mEv as
// the bond passes through equilibrium. Envelope-based sampling skips
// that confound.

TEST(AIMDConservation, RHFNVEEnergyDrift200Steps) {
    auto frame = make_h2(0.78, /*v_stretch=*/0.001);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const double dt = 0.5;   // fs
    const int n_steps = 200; // 100 fs total
    const int sample_every = 10;

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.on(md::Phase::force).add(cpu::scf::SCF{config});
    sim.build();

    // Sample E_total every `sample_every` steps. Each sim.run(frame, k)
    // call advances exactly k steps, so we run sample_every steps at a
    // time and read frame.{energy,kinetic} after each chunk.
    double e_min = std::numeric_limits<double>::infinity();
    double e_max = -std::numeric_limits<double>::infinity();
    for (int s = 0; s < n_steps; s += sample_every) {
        sim.run(frame, sample_every);
        const double E = static_cast<double>(frame.meta.get(ATV::kEnergy))
                         + static_cast<double>(frame.meta.get(ATV::kKinetic));
        if (E < e_min)
            e_min = E;
        if (E > e_max)
            e_max = E;
        EXPECT_TRUE(std::isfinite(E));
    }

    const double envelope = e_max - e_min;

    // For H2 / STO-3G starting at bond 0.78 Å with v_stretch = 0.001
    // Å/fs, the initial perturbation excites the molecule by O(30 meV)
    // total (mostly PE relative to bond minimum at 0.74 Å). The true
    // conserved quantity is bounded by this envelope plus a small
    // integrator-noise floor. We set the threshold at 0.10 eV — well
    // above the natural envelope but tight enough that a real drift
    // (e.g. losing a velocity half-kick) would push past it.
    EXPECT_LT(envelope, 0.10)
        << "NVE total-energy envelope = " << envelope << " eV across "
        << n_steps << " steps — expected ≤ 0.10. "
        << "VelocityVerlet + SCF may have lost symplecticity, "
        << "or the SCF unwrap is interfering with conservation.";

    // Sanity: PE alone must always sit in the bonded H2 well, never
    // collapse onto the unphysical RHF ionic limit (≈ -4 eV at large
    // separation). Catches the wrap × non-min-image bug that this
    // test file exists to guard against.
    EXPECT_LT(e_max, -25.0)
        << "Energy escaped the bonded H2 well — wrap × min-image bug?";
    EXPECT_GT(e_min, -35.0);
}

// ── ASPC state-leakage test ─────────────────────────────────────────

TEST(AIMDConservation, ASPCMatchesFreshSinglePointMidTrajectory) {
    auto frame = make_h2(0.80, /*v_stretch=*/0.001);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const double dt = 0.5;
    const int warm_up = 25;

    // Tight corrector tolerance — at 1e-9 the ASPC corrector is forced to
    // converge to the same fixed point as a fresh SCF, so any disagreement
    // points at state leakage rather than predictor truncation.
    cpu::scf::SCF aspc_evaluator{config};
    aspc_evaluator.predictor(/*order=*/3, /*corrector_steps=*/8,
                             /*corrector_tol=*/1e-9, /*fallback=*/true);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.on(md::Phase::force).add(std::move(aspc_evaluator));
    sim.build();

    // Warm up the predictor so it is fully populated and actively
    // extrapolating at the snapshot point.
    sim.run(frame, warm_up);

    // Snapshot the current geometry (positions + velocities) and the
    // forces just produced by the propagation.
    cpu::Frame snapshot = frame;
    const double fx0_prop = static_cast<double>(frame.span(ATV::kAtom_fx)[0]);
    const double fx1_prop = static_cast<double>(frame.span(ATV::kAtom_fx)[1]);
    const double E_prop = static_cast<double>(frame.meta.get(ATV::kEnergy));

    // Stand-alone single-point on the same geometry — Evaluate driver,
    // fresh SCF (no predictor history). Forces here are the
    // ground truth for the snapshot geometry.
    cpu::scf::SCF fresh{config};
    cpu::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate).add_ref(fresh);
    eval.build();
    eval.run(snapshot); // cpu::Evaluate runs exactly one pass

    const double fx0_fresh =
        static_cast<double>(snapshot.span(ATV::kAtom_fx)[0]);
    const double fx1_fresh =
        static_cast<double>(snapshot.span(ATV::kAtom_fx)[1]);
    const double E_fresh = static_cast<double>(snapshot.meta.get(ATV::kEnergy));

    // Convert the spec's 1e-7 Ha/Bohr threshold into eV/Å:
    //   1 Ha / Bohr = units::kHartreeInEv / kBohr2A eV/Å ≈ 27.211386 / 0.529177 ≈ 51.422 eV/Å
    //   1e-7 Ha/Bohr  ≈ 5.14e-6 eV/Å
    // We allow a 2× margin (1e-5 eV/Å) for accumulated FP arithmetic
    // differences across the two converging paths. A genuine state
    // leakage would shift forces by O(1e-3 eV/Å) and trip this easily.
    const double force_tol_eVA = 1e-5;
    EXPECT_NEAR(fx0_prop, fx0_fresh, force_tol_eVA)
        << "ASPC propagation force on atom 0 disagrees with fresh single-point "
        << "at the same geometry — possible state leakage in SCF.";
    EXPECT_NEAR(fx1_prop, fx1_fresh, force_tol_eVA)
        << "ASPC propagation force on atom 1 disagrees with fresh "
           "single-point.";

    // Energy is variational (O(δP²)), so it should agree at least one
    // order tighter than the forces.
    EXPECT_NEAR(E_prop, E_fresh, 1e-6)
        << "ASPC propagation energy disagrees with fresh single-point.";
}

// ── NVE energy-drift through the molrs.Frame ingress path ───────────
//
// rust-frame-lowering-03 ac-003: the same H2 / STO-3G NVE drift check,
// but the starting geometry comes in as a molrs.Frame and is brought
// across by `cpu::materialize` into ONE `cpu::Frame` that the whole
// AIMD trajectory then runs on. This pins that the new molrs-Frame
// ingress preserves the AIMD shared-frame invariant (MD + SCF
// on a single cpu::Frame) and energy conservation.

TEST(AIMDConservation, RHFNVEEnergyDriftViaMaterializedFrame) {
    // Geometry arrives as a molrs.Frame, exactly as the Python path
    // would deliver it, then is materialized once.
    auto molrs_frame = make_h2_molrs(0.78, /*v_stretch=*/0.001);
    cpu::Frame frame = cpu::materialize(*molrs_frame);

    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const double dt = 0.5;   // fs
    const int n_steps = 200; // 100 fs total
    const int sample_every = 10;

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.on(md::Phase::force).add(cpu::scf::SCF{config});
    sim.build();

    // The whole trajectory advances ONE materialized cpu::Frame —
    // consecutive sim.run() calls continue the same simulation.
    double e_min = std::numeric_limits<double>::infinity();
    double e_max = -std::numeric_limits<double>::infinity();
    for (int s = 0; s < n_steps; s += sample_every) {
        sim.run(frame, sample_every);
        const double E = static_cast<double>(frame.meta.get(ATV::kEnergy))
                         + static_cast<double>(frame.meta.get(ATV::kKinetic));
        if (E < e_min)
            e_min = E;
        if (E > e_max)
            e_max = E;
        EXPECT_TRUE(std::isfinite(E));
    }

    const double envelope = e_max - e_min;
    // Same 0.10 eV bound as RHFNVEEnergyDrift200Steps — the materialized
    // ingress must not perturb conservation.
    EXPECT_LT(envelope, 0.10)
        << "NVE total-energy envelope = " << envelope
        << " eV across a materialized-frame AIMD run — expected ≤ 0.10.";
    EXPECT_LT(e_max, -25.0)
        << "Energy escaped the bonded H2 well after materialization.";
    EXPECT_GT(e_min, -35.0);
}
