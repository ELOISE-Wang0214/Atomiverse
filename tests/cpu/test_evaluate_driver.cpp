/// @file tests/cpu/test_evaluate_driver.cpp
/// @brief Phase 5 — `cpu::Evaluate` driver behaves correctly under the
///        one-pass evaluate contract:
///
///   1. Energy-only `SCF(.with_gradient(false))` is accepted by
///      `Evaluate` (capability check passes) but rejected by `Dynamics`
///      (no Force capability).
///   2. `Evaluate + SCF(.with_gradient(false))` reproduces a
///      direct `SCFSolver.run()` energy bit-for-bit on H₂O / STO-3G.
///   3. The free function `cpu::single_point(frame, evaluator)` is
///      equivalent to constructing an `Evaluate`, adding the evaluator,
///      building, and running.
///
/// Acceptance: ac-009 in `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#ifdef ATV_HAS_RHF

#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/hf.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/driver/single_point.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solver.h"
#include "atomiverse/cpu/scf/scf_types.h"

using namespace ATV;

namespace {

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = 0.7572;
    f.span(ATV::kAtom_x)[2] = -0.7572;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.0;
    f.span(ATV::kAtom_y)[2] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_z)[1] = -0.5860;
    f.span(ATV::kAtom_z)[2] = -0.5860;
    f.span(ATV::kAtom_type)[0] = 8;
    f.span(ATV::kAtom_type)[1] = 1;
    f.span(ATV::kAtom_type)[2] = 1;
    f.span(ATV::kAtom_mass)[0] = 15.999;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.span(ATV::kAtom_mass)[2] = 1.008;
    f.simbox = ATV::cpu::Box::cubic(50.0);
    return f;
}

} // namespace

TEST(EvaluateDriver, EnergyOnlyEvaluatorAcceptedByEvaluate) {
    auto frame = make_h2o();
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(integrals);

    cpu::scf::SCF scf{cfg};
    scf.with_gradient(false);
    // Kind HF is enough for Evaluate; Capability write-sets are gone.
    EXPECT_NE(dynamic_cast<const cpu::HF*>(&scf), nullptr);

    cpu::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate).add_ref(scf);
    EXPECT_NO_THROW(eval.build());
}

TEST(EvaluateDriver, EnergyOnlyEvaluatorAcceptedByDynamicsAsHFKind) {
    // MD validates energy/force *kinds* (Potential/HF/…), not a Force bit.
    // Energy-only SCF still satisfies the kind gate; forces_requested is a
    // launch-time channel, not a build-time Capability.
    auto frame = make_h2o();
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(integrals);

    cpu::scf::SCF scf{cfg};
    scf.with_gradient(false);

    cpu::md::MD d;
    d.on(ATV::md::Phase::initial_integrate).add(cpu::md::VelocityVerlet{0.5});
    d.on(ATV::md::Phase::force).add_ref(scf);
    EXPECT_NO_THROW(d.build());
}

TEST(EvaluateDriver, MatchesSCFSolverDirectOnH2O) {
    // Path 1 — direct SCFSolver.
    auto frame_direct = make_h2o();
    auto integrals =
        cpu::scf::build_integrals(frame_direct, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(integrals);

    cpu::scf::SCFSolver solver{cfg};
    solver.build(integrals);
    const auto& result_direct = solver.run(integrals);
    const double E_direct_Ha = result_direct.energy;

    // Path 2 — Evaluate + SCF.
    auto frame_eval = make_h2o();
    cpu::scf::SCF scf_eval{cfg};
    scf_eval.with_gradient(false);

    cpu::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate).add_ref(scf_eval);
    eval.build();
    eval.run(frame_eval);

    const double E_via_eval_Ha = scf_eval.last_result().energy;
    EXPECT_NEAR(E_via_eval_Ha, E_direct_Ha, 1e-10)
        << "Evaluate + SCF should match SCFSolver direct to FP noise";

    // frame.meta.get(ATV::kEnergy) is populated in eV; spot-check the conversion.
    const double frame_eV = frame_eval.meta.get(ATV::kEnergy);
    const double frame_Ha = frame_eV / units::kHartreeInEv;
    EXPECT_NEAR(frame_Ha, E_via_eval_Ha, 1e-10)
        << "frame.meta.get(ATV::kEnergy) (eV) should be SCF energy (Ha) × "
           "Hartree2eV";
}

TEST(EvaluateDriver, SinglePointHelperEquivalentToInline) {
#ifdef _OPENMP
    // Pin the OMP reduction order so the two SCF solves below are
    // bit-comparable: under a heavily parallel ctest schedule the
    // thread-arrival variance in the Fock reductions exceeds the 1e-10
    // equivalence budget (observed 1.2e-10 in the 154-test suite).
    // Same determinism precedent as test_cpu_scf_spin_resolution.
    omp_set_num_threads(1);
#endif
    auto frame_inline = make_h2o();
    auto frame_helper = make_h2o();
    auto integrals =
        cpu::scf::build_integrals(frame_inline, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(integrals);

    // Inline construction.
    cpu::scf::SCF scf_inline{cfg};
    scf_inline.with_gradient(false);
    {
        cpu::Evaluate eval;
        eval.on(ATV::evaluate::Phase::evaluate).add_ref(scf_inline);
        eval.build();
        eval.run(frame_inline);
    }

    // Free-function helper.
    cpu::scf::SCF scf_helper{cfg};
    scf_helper.with_gradient(false);
    cpu::single_point(frame_helper, scf_helper);

    EXPECT_NEAR(scf_inline.last_result().energy,
                scf_helper.last_result().energy, 1e-10)
        << "single_point() must be equivalent to inline Evaluate construction "
        << "(1e-10 — parallel-reduction order can vary by ~1e-12)";
    EXPECT_NEAR(frame_inline.meta.get(ATV::kEnergy),
                frame_helper.meta.get(ATV::kEnergy), 1e-10);
}

#endif // ATV_HAS_RHF
