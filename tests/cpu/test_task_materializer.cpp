/// @file tests/cpu/test_task_materializer.cpp
/// @brief Unit tests for `cpu::task::TaskMaterializer` — the rank-agnostic
///        rebuilder that turns a TaskSpec into a runnable pipeline.
///
/// Spec: .claude/specs/mpi-launcher-embedded-python-02-codec.md. The fidelity
/// contract under test: the driver materialized from a TaskSpec, consumed via
/// the single task-kind-agnostic `run(frame)` entry, is THE SAME object a user
/// would hand-build — so a 5-step LJ MD run and an AM1 single point must agree
/// BITWISE (EXPECT_EQ, no tolerance) with directly-built `cpu::md::MD` /
/// `cpu::Evaluate` pipelines on identically constructed frames. Serial,
/// sub-second, no MPI, no I/O.

#include <gtest/gtest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <optional>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/bond_harmonic.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/cpu/semiempirical/am1_evaluator.h"
#include "atomiverse/cpu/task/task_materializer.h"
#include "atomiverse/cpu/task/task_spec.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;
using namespace ATV::cpu::task;

namespace {

// The bitwise-equality claims below compare two runs of the SAME pipeline
// configuration inside one process. Under OpenMP (> 1 thread) the force
// reduction order varies between runs, so bitwise equality only holds
// single-threaded — pin it before any parallel region fires (the
// test_evaluate_driver / test_scf_lindep convention; setenv would be too
// late, libgomp snapshots OMP_NUM_THREADS in its ELF constructor). This
// does not weaken the fidelity contract (materializer ≡ hand-built
// driver); it only removes scheduler noise from the comparison.
#ifdef _OPENMP
[[maybe_unused]] const bool kOmpPinned = [] {
    omp_set_num_threads(1);
    return true;
}();
#endif

// ── LJ fixture ──────────────────────────────────────────────────────────

// Argon-like LJ parameters; kSpacing ≈ 2^(1/6)·σ so the cluster starts near
// the LJ minimum and 5 steps at 1 fs stay tame.
constexpr double kEps = 0.0104;   // [eV]
constexpr double kSigma = 3.4;    // [Å]
constexpr double kCutoff = 6.0;   // [Å]
constexpr double kSkin = 0.3;     // [Å]
constexpr double kDt = 1.0;       // [fs]
constexpr double kSpacing = 3.8;  // [Å]
constexpr double kTinit = 10.0;   // [K]
constexpr unsigned kSeed = 42u;
constexpr int kNSteps = 5;

/// @brief Deterministic 10-atom LJ cluster (two rows of five along x) in a
///        30 Å cubic box — literals only, so two constructions are bitwise
///        identical (the test_bcast_frame reference-frame pattern).
cpu::Frame make_lj_frame() {
    constexpr int n = 10;
    cpu::Frame frame(n);
    for (int i = 0; i < n; ++i) {
        const int col = i % 5;
        const int row = i / 5;
        frame.atoms.x[i] = 8.0 + kSpacing * static_cast<double>(col);
        frame.atoms.y[i] = 10.0 + kSpacing * static_cast<double>(row);
        frame.atoms.z[i] = 15.0;
        frame.atoms.mass[i] = 39.948;
        frame.atoms.type[i] = 18;
    }
    frame.box = Mat3::cubic(30.0);
    return frame;
}

/// @brief LJ cluster plus a two-atom harmonic bond (atoms 0–1, bond type 0).
///        Tag map populated first so `add_bond` can resolve global tags,
///        exactly as in tests/mpi/test_bcast_frame.cpp / test_pe_single_chain.
cpu::Frame make_bonded_lj_frame() {
    cpu::Frame frame = make_lj_frame();
    frame.atoms.map_clear();
    frame.atoms.map_set(0, frame.atoms.nlocal);
    frame.atoms.add_bond(/*ti=*/0, /*tj=*/1, /*type=*/0);
    frame.atoms.build_special_list();
    return frame;
}

/// @brief Bitwise frame agreement: scalars.energy and every atoms.x/y/z.
void expect_frames_bitwise_equal(const cpu::Frame& got,
                                 const cpu::Frame& want) {
    EXPECT_EQ(got.scalars.energy, want.scalars.energy);
    ASSERT_EQ(got.atoms.size(), want.atoms.size());
    const int n = want.atoms.size();
    for (int i = 0; i < n; ++i) {
        SCOPED_TRACE("atom " + std::to_string(i));
        EXPECT_EQ(got.atoms.x[i], want.atoms.x[i]);
        EXPECT_EQ(got.atoms.y[i], want.atoms.y[i]);
        EXPECT_EQ(got.atoms.z[i], want.atoms.z[i]);
    }
}

// ── (a) MD equivalence: LJ + NeighborList + VelocityVerlet + InitVelocity ─

TEST(TaskMaterializer, MdLjRunMatchesHandBuiltDriverBitwise) {
    // TaskSpec mirror of the hand-built registration below, entry for entry.
    std::vector<ComponentEntry> entries;
    entries.push_back(ComponentEntry{
        ComponentParams{NeighborListParams{kCutoff, kSkin}},
        ExplicitBinding{ATV::md::Phase::force, -10}});
    entries.push_back(ComponentEntry{
        ComponentParams{PairLJCutParams{kEps, kSigma, kCutoff}},
        ExplicitBinding{ATV::md::Phase::force, 0}});
    entries.push_back(
        ComponentEntry{ComponentParams{VelocityVerletParams{kDt}},
                       std::nullopt});
    entries.push_back(
        ComponentEntry{ComponentParams{InitVelocityParams{kTinit, kSeed}},
                       std::nullopt});
    const TaskSpec spec(std::move(entries), kNSteps);

    cpu::Frame materialized_frame = make_lj_frame();
    TaskMaterializer mat(spec);
    EXPECT_EQ(mat.kind(), TaskKind::MD);
    mat.run(materialized_frame);

    // Hand-built driver with identical ctor args and identical registration
    // calls (the minimal working LJ wiring of tests/cpu/test_md_no_alloc.cpp
    // and tests/mpi/test_lj_nve.cpp).
    cpu::Frame direct_frame = make_lj_frame();
    cpu::md::MD sim;
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::NeighborList<>{kCutoff, kSkin}, -10);
    sim.on(ATV::md::Phase::force)
        .add(cpu::PairLJCut{kEps, kSigma, kCutoff}, 0);
    sim.add(cpu::md::VelocityVerlet{kDt});
    sim.add(cpu::md::InitVelocity{kTinit, kSeed});
    sim.build();
    sim.run(direct_frame, kNSteps);

    ASSERT_NE(direct_frame.scalars.energy, 0.0)
        << "fixture produced zero LJ energy — the equivalence check would "
           "be vacuous";
    expect_frames_bitwise_equal(materialized_frame, direct_frame);
}

// ── (b) MD variant: + BondHarmonic (explicit) + NoseHoover (default) ─────

TEST(TaskMaterializer, MdBondedThermostattedRunMatchesHandBuiltDriverBitwise) {
    std::vector<ComponentEntry> entries;
    entries.push_back(ComponentEntry{
        ComponentParams{NeighborListParams{kCutoff, kSkin}},
        ExplicitBinding{ATV::md::Phase::force, -10}});
    entries.push_back(ComponentEntry{
        ComponentParams{PairLJCutParams{kEps, kSigma, kCutoff}},
        ExplicitBinding{ATV::md::Phase::force, 0}});
    entries.push_back(ComponentEntry{
        ComponentParams{BondHarmonicParams{{2.0}, {kSpacing}}},
        ExplicitBinding{ATV::md::Phase::force, 1}});
    entries.push_back(
        ComponentEntry{ComponentParams{VelocityVerletParams{kDt}},
                       std::nullopt});
    entries.push_back(
        ComponentEntry{ComponentParams{InitVelocityParams{kTinit, kSeed}},
                       std::nullopt});
    entries.push_back(ComponentEntry{
        ComponentParams{NoseHooverParams{kTinit, 100.0, kDt}},
        std::nullopt});
    const TaskSpec spec(std::move(entries), kNSteps);

    cpu::Frame materialized_frame = make_bonded_lj_frame();
    TaskMaterializer mat(spec);
    mat.run(materialized_frame);

    cpu::Frame direct_frame = make_bonded_lj_frame();
    cpu::md::MD sim;
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::NeighborList<>{kCutoff, kSkin}, -10);
    sim.on(ATV::md::Phase::force)
        .add(cpu::PairLJCut{kEps, kSigma, kCutoff}, 0);
    sim.on(ATV::md::Phase::force)
        .add(cpu::BondHarmonic{{2.0}, {kSpacing}}, 1);
    sim.add(cpu::md::VelocityVerlet{kDt});
    sim.add(cpu::md::InitVelocity{kTinit, kSeed});
    sim.add(cpu::md::NoseHoover{kTinit, 100.0, kDt});
    sim.build();
    sim.run(direct_frame, kNSteps);

    expect_frames_bitwise_equal(materialized_frame, direct_frame);
}

// ── (c) SinglePoint AM1 equivalence on H2 ───────────────────────────────

/// @brief Fixed-geometry H2 (0.74 Å, 20 Å box) — the exact fixture shape of
///        tests/cpu/semiempirical/test_am1_scaffold.cpp
///        `RunConvergesClosedShellH2AgainstReferenceChargeOracle`.
cpu::Frame make_h2_frame() {
    cpu::Frame frame(2);
    frame.atoms.type[0] = 1; // Z = 1 (hydrogen)
    frame.atoms.type[1] = 1;
    frame.atoms.x[0] = 0.0;
    frame.atoms.y[0] = 0.0;
    frame.atoms.z[0] = 0.0;
    frame.atoms.x[1] = 0.0;
    frame.atoms.y[1] = 0.0;
    frame.atoms.z[1] = 0.74;
    frame.box = Mat3::cubic(20.0);
    return frame;
}

TEST(TaskMaterializer, SinglePointAm1H2MatchesHandBuiltEvaluateBitwise) {
    // AM1EvaluatorParams mirroring the NDDOConfig defaults, spelled out so
    // the two pipelines demonstrably share one configuration.
    AM1EvaluatorParams params{};
    params.charge = 0;
    params.multiplicity = 1;
    params.max_iter = 200;
    params.density_tol = 1.0e-8;
    params.energy_tol = 1.0e-10;
    const TaskSpec spec(params);

    cpu::Frame materialized_frame = make_h2_frame();
    TaskMaterializer mat(spec);
    EXPECT_EQ(mat.kind(), TaskKind::SinglePoint);
    mat.run(materialized_frame);

    // Directly-built single point: Evaluate + AM1Evaluator bound at
    // evaluate::Phase::evaluate (AM1Evaluator deliberately has no
    // default_bindings()), same config field for field.
    cpu::semiempirical::AM1Config cfg{};
    cfg.charge = 0;
    cfg.multiplicity = 1;
    cfg.max_iter = 200;
    cfg.density_tol = 1.0e-8;
    cfg.energy_tol = 1.0e-10;

    cpu::Frame direct_frame = make_h2_frame();
    cpu::Evaluate eval;
    eval.on(ATV::evaluate::Phase::evaluate)
        .add(cpu::semiempirical::AM1Evaluator{cfg});
    eval.build();
    eval.run(direct_frame);

    ASSERT_NE(direct_frame.scalars.energy, 0.0)
        << "AM1 H2 single point produced zero energy — the equivalence "
           "check would be vacuous";
    EXPECT_EQ(materialized_frame.scalars.energy, direct_frame.scalars.energy);
    ASSERT_EQ(materialized_frame.atoms.size(), direct_frame.atoms.size());
    for (int i = 0; i < direct_frame.atoms.size(); ++i) {
        SCOPED_TRACE("atom " + std::to_string(i));
        EXPECT_EQ(materialized_frame.atoms.charge[i],
                  direct_frame.atoms.charge[i]);
    }
}

// ── (d) kind() accessor ─────────────────────────────────────────────────

TEST(TaskMaterializer, KindReportsTaskKindOfTheSpec) {
    // Minimal VALID MD composition: MD::build() (fired by the materializer
    // ctor) demands an Integrator plus a Force-capable component.
    std::vector<ComponentEntry> entries;
    entries.push_back(ComponentEntry{
        ComponentParams{NeighborListParams{kCutoff, kSkin}},
        ExplicitBinding{ATV::md::Phase::force, -10}});
    entries.push_back(ComponentEntry{
        ComponentParams{PairLJCutParams{kEps, kSigma, kCutoff}},
        ExplicitBinding{ATV::md::Phase::force, 0}});
    entries.push_back(
        ComponentEntry{ComponentParams{VelocityVerletParams{kDt}},
                       std::nullopt});
    const TaskSpec md_spec(std::move(entries), 1);
    EXPECT_EQ(TaskMaterializer(md_spec).kind(), TaskKind::MD);

    const TaskSpec sp_spec(AM1EvaluatorParams{});
    EXPECT_EQ(TaskMaterializer(sp_spec).kind(), TaskKind::SinglePoint);
}

} // namespace
