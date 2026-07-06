/// @file tests/cpu/mc/test_mc_pair_distribution.cpp
/// @brief Physics gate (mc-sampling-03-validate ac-001): the sampler draws
///        configurations ∝ exp(−βE). Reference-free — a two-particle run's
///        separation distribution must match the analytic P(r) ∝ r²·exp(−βV(r))
///        over a window fully inside the minimum-image cell, and the sampled
///        ⟨V⟩ must equal the 1-D quadrature ⟨V⟩ within block standard error.
///
/// E2E: drives MC::run() over many sweeps — registered behind ATV_RUN_E2E_TESTS,
/// not the ms-scale default gate.

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
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
constexpr double kEps = 0.0104, kSig = 3.4, kCut = 8.0, kSkin = 2.0;

double pair_distance(const cpu::Frame& f) {
    double dx = static_cast<double>(f.atoms.x[0] - f.atoms.x[1]);
    double dy = static_cast<double>(f.atoms.y[0] - f.atoms.y[1]);
    double dz = static_cast<double>(f.atoms.z[0] - f.atoms.z[1]);
    f.box.min_image(dx, dy, dz);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
} // namespace

TEST(MCPairDistribution, MatchesBoltzmannRadialDensity) {
    // Small box L=12 ⇒ minimum-image half-width L/2=6: the isolated pair stays
    // within the analysis window instead of drifting apart into free space,
    // and the window [r_lo,r_hi] sits inside L/2 so the radial phase space is
    // exactly ∝ r². cut=5.5 (<L/2) defines the shifted potential.
    const double L = 12.0, cut = 5.5, skin = 0.4;
    const double T = 100.0; // K
    const double beta = 1.0 / (ATV::units::kB_ev_per_k * T);
    const double r_lo = 3.0, r_hi = 5.3;

    cpu::Frame frame(2);
    const double r_min = std::pow(2.0, 1.0 / 6.0) * kSig;
    frame.atoms.x[0] = 0.0; frame.atoms.y[0] = 0.0; frame.atoms.z[0] = 0.0;
    frame.atoms.x[1] = r_min; frame.atoms.y[1] = 0.0; frame.atoms.z[1] = 0.0;
    frame.atoms.mass[0] = frame.atoms.mass[1] = 39.948;
    frame.box = Mat3::cubic(L);
    frame.scalars.energy = 0.0;

    cpu::mc::MC mc{T, 20260703u};
    mc.set_pair(cpu::PairLJCut{kEps, kSig, cut});
    mc.set_neighbor(cpu::md::NeighborList<>{cut, skin, 1, 0, true});
    mc.set_mover(cpu::mc::TranslateMove{0.35});
    mc.build();

    mc.run(frame, 40000); // equilibrate

    const int n_samp = 300000, stride = 2;
    const double kCutV = cut; // shifted-potential cutoff for the analytic V
    std::vector<double> v_win; // V samples inside [r_lo,r_hi]
    std::vector<int> hist(20, 0);
    const double bw = (r_hi - r_lo) / hist.size();
    long long n_win = 0;
    for (int k = 0; k < n_samp; ++k) {
        mc.run(frame, stride);
        const double r = pair_distance(frame);
        if (r >= r_lo && r < r_hi) {
            v_win.push_back(mctest::lj_shifted(r, kEps, kSig, kCutV));
            hist[static_cast<int>((r - r_lo) / bw)]++;
            ++n_win;
        }
    }
    ASSERT_GT(n_win, 10000) << "too few samples in the analysis window";

    // (1) PRIMARY — ⟨V⟩ over the window vs the analytic quadrature, within
    //     block SE. This is the direct proof the sampler draws ∝ exp(−βE).
    auto [vmean, vse] = mctest::block_mean_se(v_win, 20);
    const double vquad =
        mctest::mean_V_quadrature(beta, kEps, kSig, kCutV, r_lo, r_hi);
    EXPECT_NEAR(vmean, vquad, 5.0 * vse + 1e-4)
        << "sampled <V>=" << vmean << " quad=" << vquad << " SE=" << vse;

    // (2) SANITY — histogram shape vs analytic r²·e^{−βV}: reduced χ² ~ O(1).
    const double norm =
        mctest::integrate([&](double r) {
            return mctest::radial_density(r, beta, kEps, kSig, kCutV);
        }, r_lo, r_hi, 20000);
    double chi2 = 0.0;
    for (std::size_t b = 0; b < hist.size(); ++b) {
        const double a = r_lo + b * bw;
        const double exp_frac =
            mctest::integrate([&](double r) {
                return mctest::radial_density(r, beta, kEps, kSig, kCutV);
            }, a, a + bw, 400) / norm;
        const double expc = exp_frac * static_cast<double>(n_win);
        if (expc > 5.0) {
            const double d = hist[b] - expc;
            chi2 += d * d / expc;
        }
    }
    const double red_chi2 = chi2 / (hist.size() - 1);
    EXPECT_LT(red_chi2, 5.0) << "reduced chi^2=" << red_chi2;
}
