/// @file bench_aimd_step.cpp
/// @brief Per-step wall-time benchmarks for AIMD on H2 / RHF / STO-3G.
///
/// Measures three timings on the same H2 fixture so the user can read off
/// the ASPC payoff:
///
///   BM_AIMD_NoPredict       — one MD step, fresh SCF every step
///   BM_AIMD_ASPC3           — one MD step, ASPC predictor (order 3)
///   BM_SinglePoint_Gradient — stand-alone single-point + analytic gradient
///
/// The Phase 5 spec (.claude/specs/phase5-aimd.md) sets the AIMD-payoff
/// criterion at "≤ 1.2× single gradient" — i.e. the ASPC-warmstarted MD
/// step should cost at most 20% more than a fresh single-point because
/// each subsequent SCF should converge in ~1–2 iterations after a
/// well-extrapolated guess. The benchmark prints both per-step times and
/// the ratio so a human can spot regressions.
///
/// This is a manual benchmark target; CI does not gate on the ratio
/// because it depends on host CPU, basis set, and OpenMP thread count.

#ifdef ATV_HAS_RHF

#include <benchmark/benchmark.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/evaluate.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

namespace {

cpu::Frame make_h2() {
    constexpr double L = 100.0;
    constexpr double cx = L * 0.5;
    cpu::Frame f(2);
    f.atoms.x = {cx - 0.37, cx + 0.37};
    f.atoms.y = {cx, cx};
    f.atoms.z = {cx, cx};
    f.atoms.vx = {-0.001, 0.001};
    f.atoms.vy = {0.0, 0.0};
    f.atoms.vz = {0.0, 0.0};
    f.atoms.mass = {1.008, 1.008};
    f.atoms.type = {1, 1};
    f.box = Mat3::cubic(L);
    return f;
}

cpu::scf::SCFConfig make_config(const cpu::Frame& frame) {
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    return cpu::scf::rhf_config(integrals);
}

} // namespace

// One MD step (initial + force + final + KE) with a fresh SCF every step.
static void BM_AIMD_NoPredict(benchmark::State& state) {
    auto frame = make_h2();
    auto config = make_config(frame);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();

    // Warm-up: one step so component setup costs are out of the timed loop.
    sim.run(frame, 1);

    for (auto _ : state) {
        sim.run(frame, 1);
        benchmark::DoNotOptimize(frame.scalars.energy);
    }
}
BENCHMARK(BM_AIMD_NoPredict)->Unit(benchmark::kMicrosecond);

// One MD step with ASPC density warm-start (order 3 is the spec default).
static void BM_AIMD_ASPC3(benchmark::State& state) {
    auto frame = make_h2();
    auto config = make_config(frame);

    cpu::scf::SCFEvaluator scf{config};
    scf.predictor(/*order=*/3);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(md::Phase::force).add(std::move(scf));
    sim.build();

    // Pre-load the predictor with several steps' worth of density so
    // it is actively extrapolating during the timed loop.
    sim.run(frame, 5);

    for (auto _ : state) {
        sim.run(frame, 1);
        benchmark::DoNotOptimize(frame.scalars.energy);
    }
}
BENCHMARK(BM_AIMD_ASPC3)->Unit(benchmark::kMicrosecond);

// Stand-alone single-point + analytic gradient (Evaluate driver, fresh
// SCFEvaluator each iteration). This is the "single-gradient" baseline
// that the AIMD-vs-MD-overhead ratio is measured against.
static void BM_SinglePoint_Gradient(benchmark::State& state) {
    auto frame = make_h2();
    auto config = make_config(frame);

    for (auto _ : state) {
        // A fresh Evaluate + SCFEvaluator per iteration keeps the cost
        // representative of "one cold single-point" — the comparison the
        // spec makes when judging ASPC efficacy.
        cpu::Evaluate eval;
        cpu::scf::SCFEvaluator scf{config};
        eval.on(ATV::evaluate::Phase::evaluate).add(std::move(scf));
        eval.build();
        eval.run(frame);
        benchmark::DoNotOptimize(frame.scalars.energy);
    }
}
BENCHMARK(BM_SinglePoint_Gradient)->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();

#else
int main() { return 0; }
#endif
