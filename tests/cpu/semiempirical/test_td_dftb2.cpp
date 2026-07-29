/// @file tests/cpu/semiempirical/test_td_dftb2.cpp
/// Unit gates for DFTB2::td_casida — the linear-response TD-DFTB (Casida)
/// assembly of Niehaus et al., Phys. Rev. B 63, 085108 (2001).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"
#include "atomiverse/cpu/semiempirical/td_dftb2.h"

using namespace ATV;

namespace {

cpu::Frame make_h2(double r_ang = 0.74) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = static_cast<StorageReal>(r_ang);
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    for (int i = 0; i < 2; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    return f;
}

cpu::Frame make_water() {
    cpu::Frame f(3);
    const double xs[3] = {0.0, 1.05, -0.26};
    const double ys[3] = {0.0, 0.0, 1.0};
    const int zs[3] = {8, 1, 1};
    for (int i = 0; i < 3; ++i) {
        f.span(kAtom_type)[i] = zs[i];
        f.span(kAtom_x)[i] = static_cast<StorageReal>(xs[i]);
        f.span(kAtom_y)[i] = static_cast<StorageReal>(ys[i]);
        f.span(kAtom_z)[i] = 0;
        f.span(kAtom_mass)[i] = (i == 0) ? 15.999 : 1.008;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

cpu::semiempirical::DFTB2 tight_dftb2() {
    cpu::semiempirical::DFTB2 d;
    d.scc_tol = 1e-12;
    d.max_scc = 120;
    return d;
}

} // namespace

/// Ω must come out of the assembly symmetric: K is a contraction of the same
/// transition charges over a symmetric γ (or an on-site W), and the prefactors
/// √Δε_ia … √Δε_jb are symmetric under (ia)↔(jb).
TEST(TdDftb2, CasidaMatrixIsSymmetric) {
    auto f = make_water();
    const auto r = tight_dftb2().td_casida(f, 1);
    ASSERT_GT(r.n_pair, 1);
    double max_asym = 0.0;
    for (int p = 0; p < r.n_pair; ++p)
        for (int q = 0; q < r.n_pair; ++q)
            max_asym = std::max(
                max_asym,
                std::abs(r.casida_ev2[static_cast<std::size_t>(p) * r.n_pair + q]
                         - r.casida_ev2[static_cast<std::size_t>(q) * r.n_pair
                                        + p]));
    EXPECT_LE(max_asym, 1e-12) << "Casida matrix asymmetry [eV^2]";
}

/// Decisive structural check: with the coupling scaled to zero the Casida
/// matrix is diagonal in (ε_a − ε_i)², so every root must be a bare
/// single-particle gap.
TEST(TdDftb2, ZeroCouplingReproducesKsGaps) {
    auto f = make_water();
    const auto d = tight_dftb2();
    const auto r =
        d.td_casida(f, 6, cpu::semiempirical::TdDftbChannel::singlet, 0, 0,
                    /*coupling_scale=*/0.0);
    ASSERT_EQ(r.energy_ev.size(), r.ks_gap_ev.size());
    for (std::size_t s = 0; s < r.energy_ev.size(); ++s)
        EXPECT_NEAR(r.energy_ev[s], r.ks_gap_ev[s], 1e-10)
            << "state " << s << " [eV]";
}

/// Oscillator strengths are non-negative for singlets; triplets are dark.
TEST(TdDftb2, OscillatorStrengthSingletNonNegativeTripletZero) {
    auto f = make_water();
    const auto d = tight_dftb2();
    const auto s =
        d.td_casida(f, 4, cpu::semiempirical::TdDftbChannel::singlet);
    ASSERT_EQ(s.oscillator_strength.size(), 4u);
    ASSERT_EQ(s.transition_dipole.size(), 12u);
    double f_sum = 0;
    for (double fi : s.oscillator_strength) {
        EXPECT_GE(fi, -1e-12);
        f_sum += fi;
    }
    // At least one bright-ish state in the window (water has UV intensity).
    EXPECT_GT(f_sum, 1e-6);

    const auto t =
        d.td_casida(f, 2, cpu::semiempirical::TdDftbChannel::triplet);
    for (double fi : t.oscillator_strength)
        EXPECT_NEAR(fi, 0.0, 1e-14);
}

/// The same limit in the triplet channel — the kernel, not the channel, is
/// what the λ → 0 limit removes.
TEST(TdDftb2, ZeroCouplingReproducesKsGapsTriplet) {
    auto f = make_water();
    const auto d = tight_dftb2();
    const auto r =
        d.td_casida(f, 4, cpu::semiempirical::TdDftbChannel::triplet, 0, 0,
                    /*coupling_scale=*/0.0);
    for (std::size_t s = 0; s < r.energy_ev.size(); ++s)
        EXPECT_NEAR(r.energy_ev[s], r.ks_gap_ev[s], 1e-10)
            << "state " << s << " [eV]";
}

/// H₂ at the DFTB2 bond length has exactly one σ_g → σ_u channel; its singlet
/// excitation must be real, positive and blue-shifted above the bare gap
/// (γ_AA > γ_AB makes the singlet kernel positive).
TEST(TdDftb2, H2LowestSingletIsPositiveAboveKsGap) {
    auto f = make_h2(0.74);
    const auto r = tight_dftb2().td_casida(f, 1);
    ASSERT_EQ(r.n_pair, 1);
    ASSERT_EQ(r.energy_ev.size(), 1u);
    EXPECT_TRUE(std::isfinite(r.energy_ev[0]));
    EXPECT_GT(r.energy_ev[0], 0.0);
    EXPECT_GT(r.energy_ev[0], r.ks_gap_ev[0]);
}

/// Hund-consistent ordering on the same (i,a) pair: the singlet kernel is a
/// contraction over γ_AB > 0, the triplet kernel over the negative spin
/// constants W_A, so ω_S > Δε > ω_T.
TEST(TdDftb2, H2SingletAboveTriplet) {
    auto f = make_h2(0.74);
    const auto d = tight_dftb2();
    const auto s = d.td_casida(f, 1, cpu::semiempirical::TdDftbChannel::singlet);
    const auto t = d.td_casida(f, 1, cpu::semiempirical::TdDftbChannel::triplet);
    ASSERT_EQ(s.pair_occ[0], t.pair_occ[0]);
    ASSERT_EQ(s.pair_virt[0], t.pair_virt[0]);
    EXPECT_GT(s.energy_ev[0], t.energy_ev[0]);
    EXPECT_LT(t.energy_ev[0], t.ks_gap_ev[0]);
}

/// Σ_A q^{ia}_A = ⟨φ_i|φ_a⟩ = δ_ia, which is 0 for every occupied→virtual pair.
TEST(TdDftb2, TransitionChargesSumToZero) {
    auto f = make_water();
    const auto r = tight_dftb2().td_casida(f, 1);
    ASSERT_EQ(r.transition_charge.size(),
              static_cast<std::size_t>(r.n_pair) * r.n_atoms);
    for (int p = 0; p < r.n_pair; ++p) {
        double sum = 0.0;
        for (int a = 0; a < r.n_atoms; ++a)
            sum += r.transition_charge[static_cast<std::size_t>(p) * r.n_atoms
                                       + a];
        EXPECT_NEAR(sum, 0.0, 1e-10) << "pair " << p << " [e]";
    }
}

/// Closed form for a symmetric two-level system: c_i = (a,a), c_a = (b,−b) with
/// a = 1/√(2+2S₁₂), b = 1/√(2−2S₁₂) give q₁ = ab = 1/(2√(1−S₁₂²)) = −q₂. The
/// magnitude therefore exceeds ½ for any nonzero overlap — an analytic bound the
/// general Mulliken code path must satisfy without knowing S.
TEST(TdDftb2, H2TransitionChargeMatchesTwoLevelClosedForm) {
    auto f = make_h2(0.74);
    const auto r = tight_dftb2().td_casida(f, 1);
    ASSERT_EQ(r.n_pair, 1);
    ASSERT_EQ(r.n_atoms, 2);
    const double q1 = r.transition_charge[0];
    const double q2 = r.transition_charge[1];
    EXPECT_NEAR(q1, -q2, 1e-12);
    EXPECT_GT(std::abs(q1), 0.5);
}

/// A 1×1 Casida problem is ω²(λ) = Δε² + λ·(ω²(1) − Δε²): the kernel enters
/// exactly once and linearly.
TEST(TdDftb2, H2CouplingScaleIsLinearInOmegaSquared) {
    auto f = make_h2(0.74);
    const auto d = tight_dftb2();
    const auto full = d.td_casida(f, 1);
    const auto half =
        d.td_casida(f, 1, cpu::semiempirical::TdDftbChannel::singlet, 0, 0, 0.5);
    const double gap2 = full.ks_gap_ev[0] * full.ks_gap_ev[0];
    const double w2_full = full.energy_ev[0] * full.energy_ev[0];
    const double w2_half = half.energy_ev[0] * half.energy_ev[0];
    EXPECT_NEAR(w2_half, gap2 + 0.5 * (w2_full - gap2), 1e-10);
}

/// The active window selects the highest occupied / lowest virtual orbitals.
TEST(TdDftb2, ActiveWindowNarrowsThePairSpace) {
    auto f = make_water();
    const auto d = tight_dftb2();
    const auto full = d.td_casida(f, 1);
    const auto homo_lumo = d.td_casida(
        f, 1, cpu::semiempirical::TdDftbChannel::singlet, 1, 1);
    EXPECT_EQ(homo_lumo.n_pair, 1);
    EXPECT_LT(homo_lumo.n_pair, full.n_pair);
    EXPECT_EQ(homo_lumo.pair_occ[0], full.pair_occ.back());
    EXPECT_EQ(homo_lumo.pair_virt[0], full.pair_virt.front());
}

TEST(TdDftb2, ThreeObPathThrows) {
    auto f = make_h2(0.74);
    auto d = tight_dftb2();
    d.sk_set = cpu::semiempirical::DftbSkSet::ThreeOb;
    EXPECT_THROW(d.td_casida(f, 1), std::runtime_error);
}

TEST(TdDftb2, MultiKThrows) {
    auto f = make_h2(0.74);
    auto d = tight_dftb2();
    d.k_mesh = {2, 1, 1};
    EXPECT_THROW(d.td_casida(f, 1), std::runtime_error);
}

TEST(TdDftb2, OpenShellThrows) {
    auto f = make_h2(0.74);
    auto d = tight_dftb2();
    d.multiplicity = 3;
    EXPECT_THROW(d.td_casida(f, 1), std::runtime_error);
}

TEST(TdDftb2, TooManyStatesThrows) {
    auto f = make_h2(0.74);
    EXPECT_THROW(tight_dftb2().td_casida(f, 2), std::invalid_argument);
}
