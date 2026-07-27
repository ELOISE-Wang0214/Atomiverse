/// @file tests/cpu/evaluate/test_source_cobinding.cpp
/// @brief Unit: two energy sources bound into ONE `cpu::Evaluate` pass must
///        both land in `kEnergy` — the accumulation contract, order-independent.
///
/// Mirrors `include/atomiverse/cpu/evaluate/evaluate.h` (`ATV::cpu::Evaluate`,
/// header-only preset) at the one boundary this spec changes: what a source is
/// allowed to do to the `kEnergy` accumulator channel inside a driver-owned
/// pass. Spec: `.claude/specs/esmm-contract-05-accumulate-switch` (U7).
///
/// ## The contract under test
///
/// `esmm-contract-01` made the **driver** the owner of zeroing: `cpu::Evaluate`
/// binds `ZeroForces` @ −1000, `ZeroEnergy` @ −999, `ZeroStress` @ −998 on
/// `evaluate::Phase::prepare`. Given that, every in-plan source must only ever
/// `+=` into `kEnergy` / `kAtom_f*` / `kStress`. Today the five
/// electronic-structure evaluators (GPW, MeanFieldMethod, PW, AM1, RM1) still
/// **overwrite** with `=` — `src/cpu/scf/mean_field_method.cpp:1199,1201` for
/// the Gaussian methods. The consequence is not a rounding difference: whichever
/// co-bound source ran *first* is silently erased, so "electronic structure ×
/// molecular mechanics 自由组合" cannot be expressed at all.
///
/// ## Why an ordering pair, not one case
///
/// A single co-binding case cannot distinguish "the contract holds" from "the
/// last writer happened to be an accumulator". So both orderings of the same
/// two sources are asserted against the **same** sum:
///
///   - stub @ 0, RHF @ 10 — **RED today**: RHF's `=` clobbers the stub's 1.5 eV
///     and the pass reports the bare SCF energy.
///   - RHF @ 0, stub @ 10 — green today, but only by accident (the accumulator
///     ran last). It is the pin that the GREEN fix does not merely move the
///     overwrite around.
///
/// Two controls bracket the claim so a failure cannot be misread:
///
///   - `TwoAccumulatingStubsSum` — no SCF at all. If this ever fails, the
///     harness (driver-owned zeroing + `Evaluate`'s pass) is broken and the
///     co-binding verdicts mean nothing.
///   - `RhfAloneMatchesLiteratureEnergy` — the **sole-source equivalence
///     witness**. After a driver's Zero*, `+=` ≡ `=` for a single author, so
///     this value must be bit-for-bit unaffected by the flip. It is the only
///     assertion in this file carrying a literature number.
///
/// ## Domain values (hard-coded; no live oracle)
///
/// H₂ at r = 0.7122 Å, RHF/STO-3G. NIST CCCBDB (NIST Standard Reference
/// Database 101, <https://cccbdb.nist.gov>), *calculated* HF/STO-3G data for
/// H₂ (¹Σg⁺): r_e = 0.712 Å (0.7122 Å at table precision) and
/// E_total = −1.11751 Hartree. Geometry originally W. J. Hehre, R. F. Stewart,
/// J. A. Pople, *J. Chem. Phys.* **51**, 2657 (1969), DOI 10.1063/1.1672392.
/// The same pin backs `regressions/minimizer/rhf_h2_linesearch` and
/// `regressions/mc/mc_bias_es` (`−30.4088842472 eV = −1.117506 Hartree`).
///
/// ## Determinism
///
/// `set_threads(1)` on every driver (and `OMP_NUM_THREADS=1` /
/// `OPENBLAS_NUM_THREADS=1` pinned in CMake, as for the sibling SCF units): the
/// pinned energy must not drift with reduction order. No RNG, no wall clock, no
/// filesystem, no network, no third-party package at build or test time.
///
/// **Scope.** Unit: one `Evaluate` pass over a 2-atom molecule, milliseconds.
/// The out-of-library ES+MM witness (LJ dimer + nested RHF, both terms
/// literature-pinned) is `regressions/architecture/lj_rhf_coupled` (R1), not
/// this file.

#ifdef ATV_HAS_HF

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/potential.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/box.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

namespace {

// ── Frozen constants ────────────────────────────────────────────────────────

/// @brief Energy the classical stub contributes, per launch [eV].
///
/// Deliberately neither 1.0 nor 0.0: a clobbered stub is then distinguishable
/// from a unit-valued or absent one, and 1.5 eV is 1.5 × 10⁴ times the
/// assertion tolerance, so the RED cannot be a tolerance artefact.
constexpr double kStubEnergy = 1.5;

/// @brief H–H separation [Å] — CCCBDB HF/STO-3G r_e at table precision.
constexpr double kH2Bond = 0.7122;

/// @brief RHF/STO-3G total energy of H₂ at `kH2Bond` [eV].
///
/// CCCBDB quotes −1.11751 Hartree; the in-tree regressions pin
/// −1.117506 Hartree = −30.4088842472 eV (see the file header). Written here at
/// the literature's own resolution.
constexpr double kRhfH2Energy = -30.408884;

/// @brief The same value in Hartree, for the `rhf.energy()` cross-check.
constexpr double kRhfH2EnergyHa = -1.117506;

/// @brief Assertion tolerance on `kEnergy` [eV].
///
/// Set by the **literature**, not by the solver: CCCBDB quotes the total energy
/// to 5 decimals in Hartree (±5 × 10⁻⁶ Ha ≈ ±1.4 × 10⁻⁴ eV), so a tighter pin
/// would be asserting our own arithmetic rather than the cited value. The SCF is
/// converged to 1e-10 Ha (2.7 × 10⁻⁹ eV) and single-threaded, so the measured
/// deviation from the pin is ~2.5 × 10⁻⁷ eV — 400 × inside. Crucially it is
/// 1.5 × 10⁴ times smaller than `kStubEnergy`, which is what makes the
/// co-binding verdict unambiguous.
constexpr double kEnergyTol = 1e-4;

// ── Stub source ─────────────────────────────────────────────────────────────

/// @brief Minimal classical energy source, written the way every accumulating
///        evaluator in the tree writes: `kEnergy += ...`.
///
/// A `cpu::Potential` so `Evaluate::validate_components()` accepts it as an
/// energy source in its own right (the `TwoAccumulatingStubsSum` control has no
/// other source). It touches no force or stress channel — this file is about
/// `kEnergy` only.
class AccumulatingStub final : public cpu::Potential {
  public:
    /// @brief Construct with the energy to contribute per launch.
    /// @param e Energy added to `kEnergy` on every launch [eV].
    /// @param tag Component name, so two instances stay distinguishable.
    AccumulatingStub(double e, std::string tag)
        : e_(e), tag_(std::move(tag)) {}

    std::string name() const override { return tag_; }

    void launch(const cpu::LaunchContext& ctx) override {
        ++launches;
        ctx.frame.meta.get(kEnergy) += e_;
    }

    int launches = 0; ///< Launch counter — guards against a dead binding.

  private:
    double e_;        ///< Contribution per launch [eV].
    std::string tag_; ///< Component name.
};

// ── Geometry ────────────────────────────────────────────────────────────────

/// @brief H₂ at the CCCBDB HF/STO-3G equilibrium separation, centred at origin.
/// @return A 2-atom Frame in a 20 Å cubic box (isolated molecule).
cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = -0.5 * kH2Bond;
    f.span(kAtom_x)[1] = 0.5 * kH2Bond;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_type)[0] = 1; // Z = 1
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

/// @brief A converging RHF/STO-3G method with the gradient switched off.
///
/// `with_gradient(false)` keeps the epilogue on the pure-energy branch
/// (`mean_field_method.cpp:1201`), which is the write this spec flips; the force
/// branch (`:1191-1199`) is pinned by R1, not here.
void configure(cpu::scf::RHF& rhf) {
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-10;
    rhf.with_gradient(false);
}

} // namespace

// ── Control: the harness itself accumulates ─────────────────────────────────

/// Two accumulating sources, no SCF. `Evaluate` zeroes at `prepare` and both
/// launches land in `kEnergy`. Passes today and after the flip; if it ever
/// fails, the driver-owned zeroing / one-pass harness is broken and every
/// verdict below is meaningless.
TEST(SourceCoBinding, TwoAccumulatingStubsSum) {
    cpu::Frame frame = make_h2();

    AccumulatingStub a{kStubEnergy, "stub_a"};
    AccumulatingStub b{2.0 * kStubEnergy, "stub_b"};

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add_ref(a, /*order=*/0);
    eval.on(evaluate::Phase::evaluate).add_ref(b, /*order=*/10);
    eval.set_threads(1);
    eval.build();
    eval.run(frame);

    ASSERT_EQ(a.launches, 1);
    ASSERT_EQ(b.launches, 1);
    EXPECT_NEAR(static_cast<double>(frame.meta.get(kEnergy)),
                3.0 * kStubEnergy, kEnergyTol)
        << "two classical `+=` sources in one Evaluate pass did not sum — the "
           "harness (ZeroEnergy @ -999 + one pass) is broken, not the ES "
           "evaluators";
}

// ── Control: sole-source equivalence witness ────────────────────────────────

/// RHF as the ONLY source under `Evaluate`. After the driver's `ZeroEnergy`,
/// `+=` ≡ `=` for a single author, so this literature-pinned value must be
/// unchanged by the accumulate switch. Green today; it is the guard that the
/// GREEN fix does not perturb any sole-source configuration in the tree.
TEST(SourceCoBinding, RhfAloneMatchesLiteratureEnergy) {
    cpu::Frame frame = make_h2();

    cpu::scf::RHF rhf;
    configure(rhf);

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add_ref(rhf, /*order=*/0);
    eval.set_threads(1);
    eval.build();
    eval.run(frame);

    ASSERT_TRUE(rhf.converged())
        << "the SCF did not converge — a non-converged run would masquerade as "
           "an energy-channel bug in every case below";
    EXPECT_NEAR(rhf.energy(), kRhfH2EnergyHa, 1e-5)
        << "RHF/STO-3G H2 at r = " << kH2Bond
        << " A no longer reproduces the CCCBDB total energy [Ha]";
    EXPECT_NEAR(static_cast<double>(frame.meta.get(kEnergy)), kRhfH2Energy,
                kEnergyTol)
        << "sole-source RHF under Evaluate must write the literature energy to "
           "kEnergy [eV] — the accumulate switch may not move this number";
}

// ── The contract: co-binding, both orderings, same sum ──────────────────────

/// **RED today.** Classical stub at order 0, RHF at order 10 — the source that
/// runs first is the one that gets erased, so today `kEnergy` reports the bare
/// SCF energy and the stub's 1.5 eV is gone. This is the exact shape of ES+MM
/// co-binding (an MM term, then the QM region) and the root blocker the spec
/// names.
TEST(SourceCoBinding, StubBeforeRhfBothLandInEnergy) {
    cpu::Frame frame = make_h2();

    AccumulatingStub stub{kStubEnergy, "stub_before_rhf"};
    cpu::scf::RHF rhf;
    configure(rhf);

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add_ref(stub, /*order=*/0);
    eval.on(evaluate::Phase::evaluate).add_ref(rhf, /*order=*/10);
    eval.set_threads(1);
    eval.build();
    eval.run(frame);

    ASSERT_EQ(stub.launches, 1) << "the stub binding never ran";
    ASSERT_TRUE(rhf.converged());

    const double got = static_cast<double>(frame.meta.get(kEnergy));
    const double want = kStubEnergy + kRhfH2Energy;
    EXPECT_NEAR(got, want, kEnergyTol)
        << "co-bound sources did not sum: kEnergy = " << got << " eV, expected "
        << want << " eV = " << kStubEnergy << " (stub) + " << kRhfH2Energy
        << " (RHF/STO-3G H2). A measured value of ~" << kRhfH2Energy
        << " eV is the clobber signature — MeanFieldMethod's epilogue "
           "OVERWRITES kEnergy (`=`, mean_field_method.cpp:1199/1201) instead "
           "of accumulating, so the source that ran first is erased and ES+MM "
           "co-binding is unexpressible";
}

/// The order-independence pin. Same two sources, orders swapped. Green today —
/// but only because the accumulator happens to run last, which is not a
/// contract. It fails the moment a fix merely relocates the overwrite (e.g. by
/// having the classical source write `=`, or by zeroing inside the ES source).
TEST(SourceCoBinding, StubAfterRhfBothLandInEnergy) {
    cpu::Frame frame = make_h2();

    cpu::scf::RHF rhf;
    configure(rhf);
    AccumulatingStub stub{kStubEnergy, "stub_after_rhf"};

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add_ref(rhf, /*order=*/0);
    eval.on(evaluate::Phase::evaluate).add_ref(stub, /*order=*/10);
    eval.set_threads(1);
    eval.build();
    eval.run(frame);

    ASSERT_EQ(stub.launches, 1);
    ASSERT_TRUE(rhf.converged());

    const double got = static_cast<double>(frame.meta.get(kEnergy));
    const double want = kStubEnergy + kRhfH2Energy;
    EXPECT_NEAR(got, want, kEnergyTol)
        << "kEnergy = " << got << " eV, expected " << want
        << " eV. The sum must not depend on the binding order of the two "
           "sources — a fix that only makes the ES source write last is not "
           "the accumulation contract";
}

/// The two orderings must agree with each other, not merely each with the pin.
/// Stated separately because it is the invariant a user relies on when
/// composing sources they did not write: binding order is a scheduling detail,
/// never a physics knob.
TEST(SourceCoBinding, EnergySumIsIndependentOfBindingOrder) {
    double results[2] = {0.0, 0.0};

    for (int variant = 0; variant < 2; ++variant) {
        cpu::Frame frame = make_h2();
        AccumulatingStub stub{kStubEnergy, "stub_order_" + std::to_string(variant)};
        cpu::scf::RHF rhf;
        configure(rhf);

        cpu::Evaluate eval;
        const int stub_order = (variant == 0) ? 0 : 10;
        const int rhf_order = (variant == 0) ? 10 : 0;
        eval.on(evaluate::Phase::evaluate).add_ref(stub, stub_order);
        eval.on(evaluate::Phase::evaluate).add_ref(rhf, rhf_order);
        eval.set_threads(1);
        eval.build();
        eval.run(frame);

        ASSERT_EQ(stub.launches, 1) << "variant " << variant;
        ASSERT_TRUE(rhf.converged()) << "variant " << variant;
        results[variant] = static_cast<double>(frame.meta.get(kEnergy));
    }

    EXPECT_NEAR(results[0], results[1], kEnergyTol)
        << "stub-then-RHF gave " << results[0] << " eV but RHF-then-stub gave "
        << results[1]
        << " eV — the total energy of one Evaluate pass depends on the order "
           "the user happened to bind two sources in";
    EXPECT_TRUE(std::isfinite(results[0]));
}

#endif // ATV_HAS_HF
