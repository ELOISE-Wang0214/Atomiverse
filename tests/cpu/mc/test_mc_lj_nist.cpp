/// @file tests/cpu/mc/test_mc_lj_nist.cpp
/// @brief Scientific gate (mc-sampling-03-validate ac-003): LJ fluid ⟨U*⟩/N vs
///        the NIST SRSW NVT-MC reference at two state points. Energy is the
///        pure-cut (unshifted) truncated energy at rc=3σ PLUS the analytic
///        long-range correction — the form the SRSW data is defined against
///        (PairLJCut is cut-AND-shifted, so its raw energy is re-evaluated
///        unshifted on snapshots). Tolerance from the run's own block SE with a
///        0.05 absolute floor; block SE + tolerance + caveats are logged so the
///        pass is auditable and never cosmetically tightened.
///
/// type: scientific — verdict is owned by /mol:bench, not /mol:impl. E2E:
/// registered behind ATV_RUN_E2E_TESTS. Reduced LJ units (ε = σ = 1).

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/mc/mc.h"
#include "atomiverse/cpu/mc/translate_move.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

#include "mc_test_util.h"

using namespace ATV;

namespace {
constexpr double kRc = 3.0; // rc = 3σ (σ = 1)

// U_LRC/N = (8/3)π ρ* [ (1/3) rc^-9 − rc^-3 ].
double u_lrc_per_n(double rho) {
    return (8.0 / 3.0) * M_PI * rho
           * ((1.0 / 3.0) * std::pow(kRc, -9.0) - std::pow(kRc, -3.0));
}

struct Result { double u_per_n, se, tol, u_ref; };

Result run_state(double Tstar, double rho, double u_ref) {
    // nc=5 ⇒ N=500 (matches the SRSW reference N and, crucially, keeps the box
    // large enough that rc+skin < L/2 even at ρ*=0.90 — nc=4 puts L/2 below the
    // neighbor cutoff at high density and corrupts the minimum image).
    int n = 0;
    cpu::Frame frame = mctest::fcc(/*nc=*/5, rho, n); // 500 atoms
    const double T = Tstar / ATV::units::kB_ev_per_k; // ε=1 eV ⇒ kT[eV]=T*

    cpu::mc::MC mc{T, 424242u};
    mc.set_pair(cpu::PairLJCut{1.0, 1.0, kRc});
    mc.set_neighbor(cpu::md::NeighborList<>{kRc, 0.3, 1, 0, true});
    mc.set_mover(cpu::mc::TranslateMove{0.12});
    mc.build();

    // Equilibration budget: a dense near-coexistence point (ρ*≈0.9) must MELT
    // the FCC seed before averaging — otherwise residual order pins ⟨U⟩ too low.
    // ρ*=0.78 converges in a few hundred passes (validates the machinery);
    // ρ*=0.90 needs far more. The /mol:bench evaluator owns the final budget +
    // verdict (this scientific criterion is never flipped by /mol:impl); the
    // default here is a compromise for a gated run.
    const int passes_equil = 3000, passes_samp = 1000;
    mc.run(frame, passes_equil * n); // equilibrate

    std::vector<double> u_series; // per-pass ⟨U_cut⟩/N samples
    u_series.reserve(passes_samp);
    for (int p = 0; p < passes_samp; ++p) {
        mc.run(frame, n); // one pass = N attempted moves
        u_series.push_back(mctest::truncated_lj_energy(frame, 1.0, 1.0, kRc)
                           / static_cast<double>(n));
    }

    auto [ucut, se] = mctest::block_mean_se(u_series, 20);
    const double u_per_n = ucut + u_lrc_per_n(rho);
    const double tol = std::max(3.0 * se, 0.05);
    return {u_per_n, se, tol, u_ref};
}
} // namespace

TEST(MCLJNist, MatchesSRSWReference) {
    struct SP { double Tstar, rho, uref; };
    const SP points[] = {{0.85, 0.78, -5.5386}, {0.90, 0.90, -6.1773}};

    for (const SP& sp : points) {
        Result r = run_state(sp.Tstar, sp.rho, sp.uref);
        std::printf("[NIST] T*=%.2f rho*=%.2f  <U*>/N=%.4f (ref %.4f)  "
                    "blockSE=%.4f tol=%.4f  |dev|=%.4f\n",
                    sp.Tstar, sp.rho, r.u_per_n, r.u_ref, r.se, r.tol,
                    std::abs(r.u_per_n - r.u_ref));
        std::printf("[NIST] caveats: energy = pure-cut (rc=3σ) re-eval + LRC "
                    "(sampler runs cut-and-shifted); N=500 matches the ref; "
                    "residual = short-run statistical error + cut-vs-cut-shifted "
                    "ensemble difference.\n");
        EXPECT_LE(std::abs(r.u_per_n - r.u_ref), r.tol)
            << "T*=" << sp.Tstar << " rho*=" << sp.rho;
    }
}
