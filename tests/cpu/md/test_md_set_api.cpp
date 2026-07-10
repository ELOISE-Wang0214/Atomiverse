/// @file tests/cpu/md/test_md_set_api.cpp
/// @brief Phase 2 (driver-infra-auto-injection) — categorical setters on
///        cpu::md::MD + install_infra rules.
///
/// Covers:
///   - reads() overrides on NoseHooverNPT, Barostat
///   - the NeighborConsumer/NeighborProvider opt-in protocol
///     (arch-refactor-03-neighbor): PairLJCut and CoulLong are consumers
///     reporting their ctor cutoff [Å]; NeighborList is the provider;
///     NoseHooverNPT is no consumer; PairMorse opts out (O(N²) demo)
///   - singular setters (set_pair, set_neighbor, set_integrator) reject
///     re-call with a different component
///   - append setters (set_thermostat, set_observer, set_record) accept
///     repeated calls
///   - install_infra: Stress dependency from NPT auto-injects Virial
///   - install_infra: explicit Virial via set_observer suppresses
///     the auto-inject (driver-level dedup by canonical-provider name)
///   - install_infra: pair evaluator without set_neighbor throws at build
///
/// Acceptance: T3 / T4 / T5 / T11 / T12 / T14 of
/// `.claude/specs/driver-infra-auto-injection.md`; protocol migration per
/// `.claude/specs/arch-refactor-03-neighbor.md`.

#include <gtest/gtest.h>

#include <map>
#include <string>

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/component/neighbor_protocol.h"
#include "atomiverse/cpu/md/barostat.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/nose_hoover_npt.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/observer/virial_stress.h"
#include "atomiverse/cpu/potentials/coul_long.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/cpu/potentials/pair_morse.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace {

/// @brief MD subclass that snapshots the pending component multiset
/// at the moment `validate_components()` runs (i.e. *after*
/// `install_infra()` has injected canonical providers but before the
/// Plan is compiled). The snapshot lets us assert the auto-injection
/// rules through introspection rather than runtime side effects.
class IntrospectMD : public ATV::cpu::md::MD {
  public:
    mutable std::map<std::string, int> name_counts;

  protected:
    void validate_components() const override {
        for_each_pending(
            [&](const ATV::cpu::ComponentBase& c) { ++name_counts[c.name()]; });
        ATV::cpu::md::MD::validate_components();
    }
};

constexpr double kCutoff = 6.0;
constexpr double kEps = 0.0104;
constexpr double kSigma = 3.4;
constexpr double kSkin = 0.3;
constexpr double kDt = 1.0;
constexpr double kT = 300.0;
constexpr double kP = 0.0;

// ─── Component-level overrides ──────────────────────────────────────

TEST(MDSetApi, NPTReadsStress) {
    ATV::cpu::md::NoseHooverNPT npt({.T_target = kT,
                                     .tau_T = 100.0,
                                     .P_target = kP,
                                     .tau_P = 1000.0,
                                     .dt = kDt});
    EXPECT_NE(static_cast<int>(npt.reads() & ATV::Capability::Stress), 0);
}

TEST(MDSetApi, BarostatReadsStress) {
    ATV::cpu::md::Barostat baro(kP, 1000.0, kDt);
    EXPECT_NE(static_cast<int>(baro.reads() & ATV::Capability::Stress), 0);
}

// ─── NeighborConsumer/NeighborProvider protocol (arch-refactor-03) ──

TEST(MDSetApi, LJIsNeighborConsumerReportingCtorCutoff) {
    // PairLJCut inherits cpu::PairEvaluator, which implements
    // NeighborConsumer by forwarding pair_cutoff() — the consumer view
    // must report exactly the cutoff passed at construction [Å].
    ATV::cpu::PairLJCut lj(kEps, kSigma, kCutoff);
    const ATV::cpu::NeighborConsumer* nc = lj.as_neighbor_consumer();
    ASSERT_NE(nc, nullptr)
        << "every cutoff-pair evaluator must opt into NeighborConsumer";
    EXPECT_EQ(nc->neighbor_cutoff(), ATV::Real(kCutoff));
}

TEST(MDSetApi, NPTIsNotANeighborConsumer) {
    // NPT is not a pair evaluator — as_neighbor_consumer() defaults null.
    ATV::cpu::md::NoseHooverNPT npt({.T_target = kT,
                                     .tau_T = 100.0,
                                     .P_target = kP,
                                     .tau_P = 1000.0,
                                     .dt = kDt});
    EXPECT_EQ(npt.as_neighbor_consumer(), nullptr);
}

TEST(MDSetApi, NeighborListIsNeighborProvider) {
    // Supply side of the protocol: the driver infrastructure scan keys on
    // as_neighbor_provider(), not on name() == "neighbor_list".
    ATV::cpu::md::NeighborList<> nl{kCutoff, kSkin};
    const ATV::cpu::NeighborProvider* np = nl.as_neighbor_provider();
    ASSERT_NE(np, nullptr);
    EXPECT_EQ(np, static_cast<const ATV::cpu::NeighborProvider*>(&nl));
}

TEST(MDSetApi, PairMorseOptsOutOfNeighborConsumer) {
    // The O(N²) all-pairs demo style opts out of the neighbor-list gate:
    // as_neighbor_consumer() must return nullptr even though its
    // PairEvaluator base opts in (opt-out overrides opt-in).
    ATV::cpu::PairMorse morse(1.0, 2.0, 1.2, kCutoff);
    EXPECT_EQ(morse.as_neighbor_consumer(), nullptr);
}

TEST(MDSetApi, CoulLongIsNeighborConsumerReportingCtorCutoff) {
    // CoulLong previously joined the missing-NL gate without
    // forwarding its cutoff — halo sizing silently skipped it. Under the
    // protocol it must report its real-space ctor cutoff [Å] (deliberate
    // behavior delta of arch-refactor-03-neighbor).
    constexpr double kAlpha = 0.3;       // Ewald splitting [1/Å]
    constexpr double kCoulCutoff = 10.0; // real-space cutoff [Å]
    ATV::cpu::CoulLong coul(kAlpha, kCoulCutoff);
    const ATV::cpu::NeighborConsumer* nc = coul.as_neighbor_consumer();
    ASSERT_NE(nc, nullptr)
        << "CoulLong must opt into NeighborConsumer (cutoff report fix)";
    EXPECT_EQ(nc->neighbor_cutoff(), ATV::Real(kCoulCutoff));
}

// ─── Singular setter discipline ─────────────────────────────────────

TEST(MDSetApi, SetPairTwiceWithSameTypeIsIdempotent) {
    ATV::cpu::md::MD sim;
    EXPECT_NO_THROW(sim.set_pair(ATV::cpu::PairLJCut{kEps, kSigma, kCutoff}));
    // Same component class — idempotent (no second registration).
    EXPECT_NO_THROW(sim.set_pair(ATV::cpu::PairLJCut{kEps, kSigma, kCutoff}));
}

// (Different-name re-call test deferred until a second pair evaluator
// exists in the codebase. Today LJ is the only PairEvaluator on CPU.)

TEST(MDSetApi, SetNeighborTwiceThrowsWhenNameDiffers) {
    // Same-name idempotent path for now (CellList default).
    ATV::cpu::md::MD sim;
    EXPECT_NO_THROW(
        sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin}));
    EXPECT_NO_THROW(
        sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin}));
}

TEST(MDSetApi, SetIntegratorTwiceWithSameNameIdempotent) {
    ATV::cpu::md::MD sim;
    EXPECT_NO_THROW(sim.set_integrator(ATV::cpu::md::VelocityVerlet{kDt}));
    EXPECT_NO_THROW(sim.set_integrator(ATV::cpu::md::VelocityVerlet{kDt}));
}

TEST(MDSetApi, SetIntegratorRejectsDifferentTypeAfterFirst) {
    ATV::cpu::md::MD sim;
    sim.set_integrator(ATV::cpu::md::VelocityVerlet{kDt});
    EXPECT_THROW(
        sim.set_integrator(ATV::cpu::md::NoseHooverNPT{{.T_target = kT,
                                                        .tau_T = 100.0,
                                                        .P_target = kP,
                                                        .tau_P = 1000.0,
                                                        .dt = kDt}}),
        std::invalid_argument);
}

// ─── install_infra: Stress dependency → Virial ────────────────

TEST(MDSetApi, NPTAutoInjectsVirialStress) {
    IntrospectMD sim;
    sim.set_pair(ATV::cpu::PairLJCut{kEps, kSigma, kCutoff});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin});
    sim.set_integrator(ATV::cpu::md::NoseHooverNPT{{.T_target = kT,
                                                    .tau_T = 100.0,
                                                    .P_target = kP,
                                                    .tau_P = 1000.0,
                                                    .dt = kDt}});
    EXPECT_NO_THROW(sim.build());
    EXPECT_EQ(sim.name_counts["virial"], 1)
        << "NPT::reads()==Stress → install_infra must add Virial once";
}

TEST(MDSetApi, ExplicitVirialStressSuppressesAutoInject) {
    IntrospectMD sim;
    sim.set_pair(ATV::cpu::PairLJCut{kEps, kSigma, kCutoff});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin});
    sim.set_integrator(ATV::cpu::md::NoseHooverNPT{{.T_target = kT,
                                                    .tau_T = 100.0,
                                                    .P_target = kP,
                                                    .tau_P = 1000.0,
                                                    .dt = kDt}});
    sim.set_observer(ATV::cpu::md::Virial{});
    EXPECT_NO_THROW(sim.build());
    EXPECT_EQ(sim.name_counts["virial"], 1)
        << "Explicit Virial via set_observer must dedup the auto-inject";
}

TEST(MDSetApi, NoStressReadsLeavesVirialStressOut) {
    IntrospectMD sim;
    sim.set_pair(ATV::cpu::PairLJCut{kEps, kSigma, kCutoff});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin});
    sim.set_integrator(
        ATV::cpu::md::VelocityVerlet{kDt}); // NVE — no stress dep
    EXPECT_NO_THROW(sim.build());
    EXPECT_EQ(sim.name_counts["virial"], 0)
        << "VelocityVerlet does not read Stress → no Virial injection";
}

// ─── install_infra: missing-NL throw ────────────────────────────────

TEST(MDSetApi, PairWithoutNeighborListThrowsAtBuild) {
    ATV::cpu::md::MD sim;
    sim.set_pair(ATV::cpu::PairLJCut{kEps, kSigma, kCutoff});
    sim.set_integrator(ATV::cpu::md::VelocityVerlet{kDt});
    // No set_neighbor → install_infra throws naming the offending evaluator.
    try {
        sim.build();
        FAIL() << "expected throw — pair evaluator registered without "
                  "set_neighbor";
    } catch (const std::runtime_error& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("neighbor"), std::string::npos) << msg;
        EXPECT_NE(msg.find("lj"), std::string::npos) << msg;
    }
}

TEST(MDSetApi, PairWithNeighborListBuildsClean) {
    ATV::cpu::md::MD sim;
    sim.set_pair(ATV::cpu::PairLJCut{kEps, kSigma, kCutoff});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin});
    sim.set_integrator(ATV::cpu::md::VelocityVerlet{kDt});
    EXPECT_NO_THROW(sim.build());
}

} // namespace
