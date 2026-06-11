/// @file benchmarks/md/pe/validation/bench_timing_capture.cpp
/// @brief Timer-instrumentation capture gate: the per-category rollup must
/// account for >= 95% of the step-loop wall time (moved out of the unit-test
/// suite: UTs assert code correctness only; wall-clock-fraction gates live in
/// benchmarks). Non-zero exit on violation.
///
/// Single-rank, bonded-only PE-like chain (no MPI halo, no neighbor list).

#include <cmath>
#include <cstdio>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/timer.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/driver/step_timing.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/potentials/angle_harmonic.h"
#include "atomiverse/cpu/potentials/bond_harmonic.h"
#include "atomiverse/cpu/potentials/dihedral_multi_harmonic.h"

using namespace ATV;
using cpu::md::StepTiming;

namespace {

inline constexpr double kKcalToEv = 1.0 / 23.060547830619026;

cpu::Frame make_chain(int n, double box_len) {
    cpu::Frame frame(n);
    const double theta0 = 109.5 * kPI / 180.0;
    const double half_apex = (kPI - theta0) / 2.0;
    const double dx = 1.53 * std::cos(half_apex);
    const double dy = 1.53 * std::sin(half_apex);
    const double x0 = (box_len - (n - 1) * dx) / 2.0;
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = x0 + i * dx;
        frame.atoms.y[i] = box_len / 2.0 + (i % 2 ? dy : 0.0);
        frame.atoms.z[i] = box_len / 2.0;
        frame.atoms.mass[i] = 14.027;
        frame.atoms.type[i] = 0;
    }
    frame.atoms.map_clear();
    frame.atoms.map_set(0, frame.atoms.nlocal);
    for (int k = 0; k < n - 1; ++k)
        frame.atoms.add_bond(k, k + 1, 0);
    for (int k = 0; k < n - 2; ++k)
        frame.atoms.add_angle(k, k + 1, k + 2, 0);
    for (int k = 0; k < n - 3; ++k)
        frame.atoms.add_dihedral(k, k + 1, k + 2, k + 3, 0);
    frame.atoms.build_special_list();
    frame.box = Mat3::cubic(box_len);
    return frame;
}

} // namespace

int main() {
    auto frame = make_chain(500, 200.0);
    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(ATV::md::Phase::force)
        .add(cpu::BondHarmonic{{268.0 * kKcalToEv}, {1.53}})
        .add(cpu::AngleHarmonic{{60.0 * kKcalToEv}, {109.5 * kPI / 180.0}})
        .add(cpu::DihedralMultiHarmonic{{{1.73 * kKcalToEv, -4.49 * kKcalToEv,
                                     0.776 * kKcalToEv, 6.99 * kKcalToEv,
                                     0.0}}});
    sim.enable_timing(true);
    sim.build();

    constexpr int kSteps = 300;
    sim.run(frame, kSteps);

    const TimerRegistry* reg = sim.timing();
    if (!reg || !reg->enabled()
        || reg->steps() != static_cast<std::uint64_t>(kSteps)
        || reg->loop_seconds() <= 0.0) {
        std::fprintf(stderr, "FAIL: timing registry missing/empty\n");
        return 2;
    }
    const auto roll = reg->rollup_by_category();
    if (roll.size() < static_cast<std::size_t>(StepTiming::Other) + 1
        || roll[StepTiming::Bonded].total_s <= 0.0
        || roll[StepTiming::Other].total_s <= 0.0) {
        std::fprintf(stderr, "FAIL: rollup categories missing\n");
        return 2;
    }
    double sum = 0.0;
    for (const auto& s : roll)
        sum += s.total_s;
    const double frac = sum / reg->loop_seconds();
    std::printf("[timing-capture] component sum / loop wall = %.1f%%\n",
                100.0 * frac);
    if (frac > 1.0 + 1e-6) {
        std::fprintf(stderr, "FAIL: component sum exceeds loop wall\n");
        return 1;
    }
    if (frac < 0.95) {
        std::fprintf(stderr,
                     "FAIL: instrumentation captured only %.1f%% (limit 95%%)\n",
                     100.0 * frac);
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
