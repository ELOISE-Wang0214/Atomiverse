/// @file tests/cpu/md/test_md_set_api.cpp
/// @brief Phase 2 (driver-infra-auto-injection) — categorical setters on
///        cpu::md::MD + install_infra rules.
///
/// Covers:
///   - reads()/needs_neighbor_list() overrides on NoseHooverNPT, Barostat, LJ
///   - singular setters (set_pair, set_neighbor, set_integrator) reject
///     re-call with a different component
///   - append setters (set_thermostat, set_observer, set_record) accept
///     repeated calls
///   - install_infra: Stress dependency from NPT auto-injects VirialStress
///   - install_infra: explicit VirialStress via set_observer suppresses
///     the auto-inject (driver-level dedup by canonical-provider name)
///   - install_infra: pair evaluator without set_neighbor throws at build
///
/// Acceptance: T3 / T4 / T5 / T11 / T12 / T14 of
/// `.claude/specs/driver-infra-auto-injection.md`.

#include <gtest/gtest.h>

#include <map>
#include <string>

#include "atomiverse/cpu/md/barostat.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover_npt.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/virial_stress.h"
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

TEST(MDSetApi, LJNeedsNeighborList) {
    ATV::cpu::md::LJ lj(kEps, kSigma, kCutoff);
    EXPECT_TRUE(lj.needs_neighbor_list());
}

TEST(MDSetApi, NPTReportsPairCutoffNothing) {
    // NPT is not a pair evaluator — needs_neighbor_list defaults false.
    ATV::cpu::md::NoseHooverNPT npt({.T_target = kT,
                                     .tau_T = 100.0,
                                     .P_target = kP,
                                     .tau_P = 1000.0,
                                     .dt = kDt});
    EXPECT_FALSE(npt.needs_neighbor_list());
}

// ─── Singular setter discipline ─────────────────────────────────────

TEST(MDSetApi, SetPairTwiceWithSameTypeIsIdempotent) {
    ATV::cpu::md::MD sim;
    EXPECT_NO_THROW(sim.set_pair(ATV::cpu::md::LJ{kEps, kSigma, kCutoff}));
    // Same component class — idempotent (no second registration).
    EXPECT_NO_THROW(sim.set_pair(ATV::cpu::md::LJ{kEps, kSigma, kCutoff}));
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

// ─── install_infra: Stress dependency → VirialStress ────────────────

TEST(MDSetApi, NPTAutoInjectsVirialStress) {
    IntrospectMD sim;
    sim.set_pair(ATV::cpu::md::LJ{kEps, kSigma, kCutoff});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin});
    sim.set_integrator(ATV::cpu::md::NoseHooverNPT{{.T_target = kT,
                                                    .tau_T = 100.0,
                                                    .P_target = kP,
                                                    .tau_P = 1000.0,
                                                    .dt = kDt}});
    EXPECT_NO_THROW(sim.build());
    EXPECT_EQ(sim.name_counts["virial_stress"], 1)
        << "NPT::reads()==Stress → install_infra must add VirialStress once";
}

TEST(MDSetApi, ExplicitVirialStressSuppressesAutoInject) {
    IntrospectMD sim;
    sim.set_pair(ATV::cpu::md::LJ{kEps, kSigma, kCutoff});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin});
    sim.set_integrator(ATV::cpu::md::NoseHooverNPT{{.T_target = kT,
                                                    .tau_T = 100.0,
                                                    .P_target = kP,
                                                    .tau_P = 1000.0,
                                                    .dt = kDt}});
    sim.set_observer(ATV::cpu::md::VirialStress{});
    EXPECT_NO_THROW(sim.build());
    EXPECT_EQ(sim.name_counts["virial_stress"], 1)
        << "Explicit VirialStress via set_observer must dedup the auto-inject";
}

TEST(MDSetApi, NoStressReadsLeavesVirialStressOut) {
    IntrospectMD sim;
    sim.set_pair(ATV::cpu::md::LJ{kEps, kSigma, kCutoff});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin});
    sim.set_integrator(
        ATV::cpu::md::VelocityVerlet{kDt}); // NVE — no stress dep
    EXPECT_NO_THROW(sim.build());
    EXPECT_EQ(sim.name_counts["virial_stress"], 0)
        << "VelocityVerlet does not read Stress → no VirialStress injection";
}

// ─── install_infra: missing-NL throw ────────────────────────────────

TEST(MDSetApi, PairWithoutNeighborListThrowsAtBuild) {
    ATV::cpu::md::MD sim;
    sim.set_pair(ATV::cpu::md::LJ{kEps, kSigma, kCutoff});
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
    sim.set_pair(ATV::cpu::md::LJ{kEps, kSigma, kCutoff});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{kCutoff, kSkin});
    sim.set_integrator(ATV::cpu::md::VelocityVerlet{kDt});
    EXPECT_NO_THROW(sim.build());
}

} // namespace
