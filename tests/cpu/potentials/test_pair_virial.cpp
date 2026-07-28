/// @file tests/cpu/potentials/test_pair_virial.cpp
/// @brief The pair virial each style writes to `kStress` is
///        \(\sum_{\text{pairs}} \Delta r \otimes f\), and its trace is
///        \(-3V\,\mathrm{d}E/\mathrm{d}V\).
///
/// Two independent statements per style, because they fail independently:
///
/// 1. **Against a direct pair sum.** The six Voigt components must equal
///    \(\sum f_{\text{pair}}\,\Delta_a \Delta_b\) recomputed here over the same
///    pair list. This is what catches a kernel that gets energy and forces
///    right and the virial wrong — the shape of the defect that shipped: the
///    `PairLJCut` xsimd ket loop guarded its six batch reduction chains with a
///    loop-invariant `if (do_virial)`, and GCC 14.3 at the project's own
///    Release flags (`-O3 -march=native`, AVX-512) dropped most of the
///    accumulation. Energy and forces were exact to 1e-17; the virial trace was
///    60 % short.
///
/// 2. **Against \(-3V\,\mathrm{d}E/\mathrm{d}V\).** Scale every coordinate and
///    the cell by the same \(\lambda\) and the virial trace is fixed by the
///    energy alone. This is what catches a style that writes no virial at all —
///    the shape `PairBorn` shipped with.
///
/// The cluster is small and sits well inside a large cell with a cutoff larger
/// than its diameter, so no pair crosses the cutoff under the finite difference
/// and the identity is exact, not asymptotic.
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/pair_list.h"
#include "atomiverse/cpu/potentials/pair_born.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

namespace {

constexpr int kN = 24;      ///< ≥ 8 neighbours for the low rows ⇒ the xsimd
                            ///< ket loop runs (W = 8 on AVX-512) plus a tail.
constexpr double kBox = 40.0;  ///< [Å] — cutoff stays < L/2.
constexpr double kCutoff = 14.0; ///< [Å] — larger than the cluster diameter.

/// @brief Deterministic pseudo-random cluster, all pairs inside the cutoff.
///
/// Spread over ~7 Å so the LJ pair force spans repulsive and attractive
/// branches; nothing sits near the cutoff, so an affine scale never changes
/// the pair set.
cpu::Frame cluster(double lambda) {
    cpu::Frame f(kN);
    for (int i = 0; i < kN; ++i) {
        const double a = 1.7 * i + 0.3;
        f.span(kAtom_x)[i] =
            lambda * (3.5 * std::sin(a) + 1.1 * std::cos(2.9 * a));
        f.span(kAtom_y)[i] =
            lambda * (3.5 * std::cos(1.3 * a) + 0.9 * std::sin(3.7 * a));
        f.span(kAtom_z)[i] =
            lambda * (3.5 * std::sin(2.1 * a + 0.5) + 1.3 * std::cos(a));
        f.span(kAtom_mass)[i] = 1.0;
        f.span(kAtom_type)[i] = 0;
        f.span(kAtom_tag)[i] = i + 1;
    }
    f.box = cpu::Box::cubic(kBox * lambda);
    return f;
}

/// @brief Dense half list: every ordered pair i < j, exactly once.
cpu::md::PairList all_pairs_half() {
    cpu::md::PairList p;
    p.inum = kN;
    p.n_rows = kN;
    p.ilist.resize(kN);
    p.numneigh.resize(kN);
    p.firstneigh.resize(kN);
    int off = 0;
    for (int i = 0; i < kN; ++i) {
        p.ilist[static_cast<std::size_t>(i)] = i;
        p.firstneigh[static_cast<std::size_t>(i)] = off;
        p.numneigh[static_cast<std::size_t>(i)] = kN - 1 - i;
        for (int j = i + 1; j < kN; ++j)
            p.neighbors.push_back(j);
        off += kN - 1 - i;
    }
    p.total_pairs = off;
    p.has_shifts = false;
    p.version = 1;
    return p;
}

struct Voigt {
    double v[6] = {0, 0, 0, 0, 0, 0};
    double trace() const { return v[0] + v[1] + v[2]; }
};

/// @brief `sum f_pair * d_a * d_b` over the same pair list, from `f_pair(r2)`.
template <class FPair>
Voigt reference_virial(const cpu::Frame& f, FPair f_pair) {
    Voigt w;
    for (int i = 0; i < kN; ++i)
        for (int j = i + 1; j < kN; ++j) {
            double dx = f.span(kAtom_x)[i] - f.span(kAtom_x)[j];
            double dy = f.span(kAtom_y)[i] - f.span(kAtom_y)[j];
            double dz = f.span(kAtom_z)[i] - f.span(kAtom_z)[j];
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= kCutoff * kCutoff)
                continue;
            const double fp = f_pair(r2);
            w.v[0] += fp * dx * dx;
            w.v[1] += fp * dy * dy;
            w.v[2] += fp * dz * dz;
            w.v[3] += fp * dx * dy;
            w.v[4] += fp * dx * dz;
            w.v[5] += fp * dy * dz;
        }
    return w;
}

/// @brief Run one pair style over the dense half list; return energy + stress.
template <class Style>
void run_style(Style& style, cpu::Frame& f, double& energy, Voigt& stress) {
    cpu::md::PairList pairs = all_pairs_half();
    cpu::Context ctx;
    ctx.md.pair_list = &pairs;
    style.build({ctx, f});
    style.launch({ctx, f, 0, 0});
    energy = f.meta.get(kEnergy);
    for (int k = 0; k < 6; ++k)
        stress.v[k] = f.meta.get(kStress)[k];
}

constexpr double kEps = 0.0104; // [eV]
constexpr double kSig = 3.4;    // [Å]

// Born–Mayer–Huggins, Tosi–Fumi Na+/Cl- unlike-pair row (see pair_born.h).
constexpr double kA = 0.21096;   // [eV]
constexpr double kRho = 0.317;   // [Å]
constexpr double kSigma = 2.755; // [Å]
constexpr double kC = 4.3691;    // [eV Å⁶]
constexpr double kD = 5.4176;    // [eV Å⁸]

} // namespace

// ── PairLJCut ───────────────────────────────────────────────────────────

TEST(PairLJCut, SimdKetLoopVirialMatchesPairSum) {
    cpu::Frame f = cluster(1.0);
    cpu::PairLJCut lj{kEps, kSig, kCutoff};
    double e = 0;
    Voigt got;
    run_style(lj, f, e, got);

    const double sig6 = std::pow(kSig, 6.0);
    const Voigt want = reference_virial(f, [&](double r2) {
        const double r6inv = sig6 / (r2 * r2 * r2);
        return 24.0 * kEps / r2 * (2.0 * r6inv * r6inv - r6inv);
    });

    // The pair sum is exact arithmetic over the identical pair set; only
    // summation order differs, so the tolerance is rounding, not physics.
    for (int k = 0; k < 6; ++k)
        EXPECT_NEAR(got.v[k], want.v[k], 1e-11 * std::max(1.0, std::fabs(want.v[k])))
            << "Voigt component " << k;
    // Guard the guard: a virial of zero must not pass by accident.
    EXPECT_GT(std::fabs(want.trace()), 1.0);
}

TEST(PairLJCut, VirialTraceIsMinusThreeVdEdV) {
    constexpr double h = 1e-6;
    double e_p = 0, e_m = 0, e_0 = 0;
    Voigt s_p, s_m, s_0;
    cpu::Frame f_p = cluster(1.0 + h);
    cpu::Frame f_m = cluster(1.0 - h);
    cpu::Frame f_0 = cluster(1.0);
    cpu::PairLJCut a{kEps, kSig, kCutoff}, b{kEps, kSig, kCutoff},
        c{kEps, kSig, kCutoff};
    run_style(a, f_p, e_p, s_p);
    run_style(b, f_m, e_m, s_m);
    run_style(c, f_0, e_0, s_0);

    const double dEdV =
        (e_p - e_m) / (f_p.box.volume() - f_m.box.volume());
    EXPECT_NEAR(s_0.trace(), -3.0 * f_0.box.volume() * dEdV,
                1e-6 * std::fabs(s_0.trace()));
}

// ── PairBorn ────────────────────────────────────────────────────────────

TEST(PairBorn, VirialMatchesPairSum) {
    cpu::Frame f = cluster(1.0);
    cpu::PairBorn born{kA, kRho, kSigma, kC, kD, kCutoff};
    double e = 0;
    Voigt got;
    run_style(born, f, e, got);

    const Voigt want = reference_virial(f, [&](double r2) {
        const double r = std::sqrt(r2);
        const double r6 = r2 * r2 * r2;
        const double r8 = r6 * r2;
        const double dUdr = -kA / kRho * std::exp((kSigma - r) / kRho)
                            + 6.0 * kC / (r6 * r) - 8.0 * kD / (r8 * r);
        return -dUdr / r;
    });

    for (int k = 0; k < 6; ++k)
        EXPECT_NEAR(got.v[k], want.v[k],
                    1e-11 * std::max(1.0, std::fabs(want.v[k])))
            << "Voigt component " << k;
    EXPECT_GT(std::fabs(want.trace()), 1e-3);
}

TEST(PairBorn, VirialTraceIsMinusThreeVdEdV) {
    constexpr double h = 1e-6;
    double e_p = 0, e_m = 0, e_0 = 0;
    Voigt s_p, s_m, s_0;
    cpu::Frame f_p = cluster(1.0 + h);
    cpu::Frame f_m = cluster(1.0 - h);
    cpu::Frame f_0 = cluster(1.0);
    cpu::PairBorn a{kA, kRho, kSigma, kC, kD, kCutoff},
        b{kA, kRho, kSigma, kC, kD, kCutoff},
        c{kA, kRho, kSigma, kC, kD, kCutoff};
    run_style(a, f_p, e_p, s_p);
    run_style(b, f_m, e_m, s_m);
    run_style(c, f_0, e_0, s_0);

    const double dEdV = (e_p - e_m) / (f_p.box.volume() - f_m.box.volume());
    EXPECT_NEAR(s_0.trace(), -3.0 * f_0.box.volume() * dEdV,
                1e-6 * std::fabs(s_0.trace()));
}

// ── Multi-type PairBorn takes the same route ───────────────────────────

TEST(PairBorn, MultiTypeTableWritesTheSameVirial) {
    cpu::Frame f_single = cluster(1.0);
    cpu::Frame f_multi = cluster(1.0);
    cpu::PairBorn single{kA, kRho, kSigma, kC, kD, kCutoff};
    cpu::PairBorn multi{
        std::vector<cpu::BornTypeCoeff>{{0, 0, kA, kRho, kSigma, kC, kD}},
        kCutoff};
    double e_s = 0, e_m = 0;
    Voigt s_s, s_m;
    run_style(single, f_single, e_s, s_s);
    run_style(multi, f_multi, e_m, s_m);
    EXPECT_NEAR(e_s, e_m, 1e-12 * std::max(1.0, std::fabs(e_s)));
    for (int k = 0; k < 6; ++k)
        EXPECT_NEAR(s_s.v[k], s_m.v[k],
                    1e-11 * std::max(1.0, std::fabs(s_s.v[k])));
}
