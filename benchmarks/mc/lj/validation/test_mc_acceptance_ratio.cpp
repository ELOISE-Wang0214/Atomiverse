/// @file tests/cpu/mc/test_mc_acceptance_ratio.cpp
/// @brief Physics gate (mc-sampling-03-validate ac-002): in a dense LJ fluid
///        the acceptance ratio decreases monotonically as the step size δ
///        grows, and a mid-range δ lands in a broad band. Tuning heuristic —
///        trend + loose band only, never a point value.
///
/// A dense fluid (not an isolated pair) is used so every trial move actually
/// perturbs the energy; isolated atoms drift into empty space where moves are
/// free and the ratio saturates near 1.
///
/// E2E: many sweeps per δ over 256 atoms — registered behind ATV_RUN_E2E_TESTS.

#include <gtest/gtest.h>

#include <cstdint>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/mc/mc.h"
#include "atomiverse/cpu/mc/translate_move.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

#include "mc_test_util.h"

using namespace ATV;

TEST(MCAcceptanceRatio, DecreasesMonotonicallyWithDelta) {
    const double rho = 0.85, Tstar = 1.2, rc = 2.5; // reduced LJ units (σ=ε=1)
    const double T = Tstar / ATV::units::kB_ev_per_k;
    const double deltas[] = {0.05, 0.1, 0.2, 0.4};
    const int nd = sizeof(deltas) / sizeof(deltas[0]);

    double prev = 2.0;
    double mid = -1.0;
    for (int i = 0; i < nd; ++i) {
        int n = 0;
        cpu::Frame frame = mctest::fcc(/*nc=*/4, rho, n); // 256 atoms

        cpu::mc::MC mc{T, 777u};
        mc.add_potential(cpu::PairLJCut{1.0, 1.0, rc});
        mc.set_neighbor(cpu::md::NeighborList<>{rc, 0.3, 1, 0, true});
        mc.set_mover(cpu::mc::TranslateMove{deltas[i]});
        mc.build();

        mc.run(frame, 50 * n); // equilibrate the lattice melt / transient
        // Re-measure over a clean window: acceptance_ratio() is lifetime, so
        // keep the equilibration short relative to the measurement.
        const long long a0 = mc.accepts(), t0 = mc.attempts();
        mc.run(frame, 150 * n);
        const double ratio = static_cast<double>(mc.accepts() - a0)
                             / static_cast<double>(mc.attempts() - t0);

        EXPECT_LE(ratio, prev + 1e-9)
            << "acceptance not monotone at delta=" << deltas[i]
            << " (ratio=" << ratio << ", prev=" << prev << ")";
        prev = ratio;
        if (deltas[i] == 0.2)
            mid = ratio;
    }

    EXPECT_GT(mid, 0.2) << "mid-delta acceptance too low: " << mid;
    EXPECT_LT(mid, 0.75) << "mid-delta acceptance too high: " << mid;
}
